use crate::{
    constants::{AUTOSOME_CONTIG_NAMES, OUT_SUFFIX_STATS},
    distributions::{expected::EstimatedDistr, trim_distr},
    extract::is_primary_read,
    util::{SampleStats, Utf8String},
};
use clap::{arg, Parser};
use ndarray::{Array1, Axis};
use ndarray_ndimage::gaussian_filter1d;
// use plotters::prelude::*;
use rand::{distributions::WeightedIndex, prelude::*};
use rand_xoshiro::Xoshiro256PlusPlus;
use rust_htslib::bam::{self, Read};
use std::{
    io::Write,
    path::{Path, PathBuf},
};

/// Extracts the read depth and insert size distribution from an alignment file
#[derive(Parser)]
pub struct StatsCommandArgs {
    /// Path to alignment file (BAM or CRAM)
    pub alignment: PathBuf,

    /// Number of regions to sample
    #[arg(short = 'n', long, default_value_t = 100)]
    pub num_regions: usize,

    /// Length of each sampled region
    #[arg(short = 'l', long, default_value_t = 100000)]
    pub region_length: u64,

    /// Minimum average mapping quality for sampled regions considered in depth distribution
    #[arg(short = 'q', long, default_value_t = 60.0)]
    pub min_depth_mapq: f64,

    /// Minimum mapping quality for individual reads considered in insert size distribution
    #[arg(short = 'm', long, default_value_t = 60)]
    pub min_insert_mapq: u8,

    /// Read depth distribution percentile cut-off. Depths accumulated above this percentile are discarded
    #[arg(long = "dc", default_value_t = 0.99)]
    pub depth_distr_cutoff: f64,

    /// Read depth distribution gaussian smoothing filter sigma value
    #[arg(long = "ds", default_value_t = 0.0)]
    pub depth_distr_gaussian_sigma: f64,

    /// The number of standard deviations from the center at which the Gaussian filter is cut off
    #[arg(long = "dt", default_value_t = 3)]
    pub depth_distr_gaussian_truncate: usize,

    /// Insert size distribution percentile cut-off. Sizes accumulated above this percentile are discarded
    #[arg(long = "ic", default_value_t = 0.998)]
    pub insert_size_distr_cutoff: f64,

    /// Insert size distribution gaussian smoothing filter sigma value
    #[arg(long = "is", default_value_t = 6.0)]
    pub insert_size_distr_gaussian_sigma: f64,

    /// The number of standard deviations from the center at which the Gaussian filter is cut off
    #[arg(long = "it", default_value_t = 3)]
    pub insert_size_distr_gaussian_truncate: usize,

    /// Contigs to inlcude in estimation (comma-separated list of contig names). Defaults to autosomes
    #[arg(
        short = 'i',
        long,
        value_delimiter = ',',
        default_value = AUTOSOME_CONTIG_NAMES,
        num_args = 1..,
    )]
    pub include_contigs: Vec<String>,

    /// Random seed
    #[arg(short = 's', long, default_value_t = 0)]
    pub seed: u64,

    /// Do not plot distributions
    #[arg(long)]
    pub no_plot: bool,

    /// Number of threads used to read the alignment file (htslib-specific)
    #[arg(short = '@', default_value_t = 0)]
    pub htslib_read_threads: usize,
}

pub fn run_stats_command(args: &StatsCommandArgs, output_prefix: &Path) -> anyhow::Result<()> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(args.seed);

    info!("Opening alignment file");
    let mut aln_rdr = bam::IndexedReader::from_path(&args.alignment)?;
    if args.htslib_read_threads > 0 {
        aln_rdr.set_threads(args.htslib_read_threads)?;
    }

    let header = aln_rdr.header().clone();

    let tids = args
        .include_contigs
        .iter()
        .filter_map(|c| header.tid(c.as_bytes()))
        .collect::<Vec<_>>();
    let tid_lens = tids
        .iter()
        .map(|&tid| header.target_len(tid).unwrap())
        .collect::<Vec<_>>();

    info!(
        "Found {}/{} of included contigs in the alignment file",
        tids.len(),
        args.include_contigs.len()
    );

    let regions = generate_random_regions(
        &tids,
        &tid_lens,
        args.num_regions,
        args.region_length,
        &mut rng,
    )?;

    let mut read_length = 0;
    let mut depth_distr = EstimatedDistr::new();
    let mut insert_distr = EstimatedDistr::new();

    info!("Estimating read depth distribution");
    for (tid, start, end) in &regions {
        aln_rdr.fetch((*tid, *start, *end))?;

        // Get average mapq
        let mapq_sum: u64 = aln_rdr.records().map(|r| r.unwrap().mapq() as u64).sum();
        let avg_mapq: f64 = mapq_sum as f64 / aln_rdr.records().count() as f64;

        if avg_mapq >= args.min_depth_mapq {
            aln_rdr.fetch((*tid, *start, *end))?;

            // Update depth distribution
            let depths = aln_rdr
                .pileup()
                .map(|p| p.unwrap().depth())
                .collect::<Vec<_>>();
            depth_distr.add_many(&depths);
        } else {
            warn!(
                "Region ({}, {}, {}) has average mapq below 60. Skipping...",
                tid, start, end
            );
        }
    }

    info!("Estimating insert size distribution");
    for (tid, start, end) in &regions {
        aln_rdr.fetch((*tid, *start, *end))?;

        for record_result in aln_rdr.records() {
            let record = record_result?;

            let not_clipped = record.cigar().to_string() == format!("{}M", record.seq_len());

            // Update insert size distribution
            if is_primary_read(&record)
                && record.is_proper_pair()
                && not_clipped
                && record.mapq() >= args.min_insert_mapq
            {
                let insert_size = record.insert_size();
                if insert_size == 0 {
                    warn!(
                        "Proper read pair ({}) has zero insert size",
                        record.qname().to_string()
                    );
                }
                insert_distr.add(insert_size.unsigned_abs() as u32);
            }

            // Update read length
            read_length = read_length.max(record.seq().len());
        }
    }

    // Smooth the distributions
    if args.depth_distr_gaussian_sigma > 0.0 {
        debug!("Smoothing depth distribution with Gaussian filter");
        let depth_distr_arr = Array1::from_vec(depth_distr.into());
        let depth_distr_arr = gaussian_filter1d(
            &depth_distr_arr,
            args.depth_distr_gaussian_sigma,
            Axis(0),
            0,
            ndarray_ndimage::BorderMode::Nearest,
            args.depth_distr_gaussian_truncate,
        );
        depth_distr = EstimatedDistr::from_vec(depth_distr_arr.to_vec());
    }

    if args.insert_size_distr_gaussian_sigma > 0.0 {
        debug!("Smoothing insert size distribution with Gaussian filter");
        let insert_distr_arr = Array1::from_vec(insert_distr.into());
        let insert_distr_arr = gaussian_filter1d(
            &insert_distr_arr,
            args.insert_size_distr_gaussian_sigma,
            Axis(0),
            0,
            ndarray_ndimage::BorderMode::Nearest,
            args.insert_size_distr_gaussian_truncate,
        );
        insert_distr = EstimatedDistr::from_vec(insert_distr_arr.to_vec());
    }

    // Trim the distributions
    if args.depth_distr_cutoff < 1.0 {
        debug!("Trimming depth distribution");
        trim_distr(depth_distr.as_mut(), args.depth_distr_cutoff);
    }

    if args.insert_size_distr_cutoff < 1.0 {
        debug!("Trimming insert size distribution");
        trim_distr(insert_distr.as_mut(), args.insert_size_distr_cutoff);
    }

    // Normalize the distributions
    depth_distr.norm();
    insert_distr.norm();

    let stats = SampleStats::new(
        read_length as u64,
        depth_distr.mean(),
        depth_distr.std_dev(),
        insert_distr.mean(),
        insert_distr.std_dev(),
        depth_distr.distr.clone(),
        insert_distr.distr.clone(),
    );

    // Write the stats to a file
    info!("Writing stats file");
    let depth_out_path = output_prefix.with_extension(OUT_SUFFIX_STATS);
    let mut depth_out_file = std::fs::File::create(depth_out_path)?;
    depth_out_file.write_all(serde_json::to_string(&stats)?.as_bytes())?;

    // Plot the distributions
    // if !args.no_plot {
    //     plot_distr(
    //         &depth_distr.distr,
    //         &output_prefix.with_extension(OUT_SUFFIX_DEPTH_DISTR_PLOT),
    //         "Read Depth",
    //     )?;

    //     plot_distr(
    //         &insert_distr.distr,
    //         &output_prefix.with_extension(OUT_SUFFIX_INSERT_DISTR_PLOT),
    //         "Insert Size",
    //     )?;
    // }

    Ok(())
}

fn generate_random_regions(
    tids: &[u32],
    tid_lens: &[u64],
    num_regions: usize,
    region_len: u64,
    rng: &mut Xoshiro256PlusPlus,
) -> anyhow::Result<Vec<(u32, i64, i64)>> {
    let weight_idx = WeightedIndex::new(tid_lens)?;

    let mut regions = Vec::new();
    for _ in 0..num_regions {
        let tid_idx = weight_idx.sample(rng);
        let tid: u32 = tids[tid_idx];
        let max_start = tid_lens[tid_idx] - region_len;
        let start = rng.gen_range(0..max_start) as i64;
        let end = start + region_len as i64;
        regions.push((tid, start, end));
    }

    Ok(regions)
}

// fn plot_distr(distr: &[f64], out_path: &PathBuf, name: &str) -> anyhow::Result<()> {
// let root = BitMapBackend::new(out_path, (800, 600)).into_drawing_area();
// root.fill(&WHITE)?;

// let mut chart = ChartBuilder::on(&root)
//     .caption(format!("{} Distribution", name), ("sans-serif", 24))
//     .margin(10)
//     .x_label_area_size(60)
//     .y_label_area_size(60)
//     .build_cartesian_2d(
//         0..distr.len() as u32,
//         0.0..distr.iter().cloned().fold(0.0, f64::max),
//     )?;

// chart
//     .configure_mesh()
//     .x_desc(name)
//     .y_desc("Probability")
//     .axis_desc_style(("sans-serif", 20))
//     .draw()?;

// chart.draw_series(LineSeries::new(
//     distr.iter().enumerate().map(|(i, &x)| (i as u32, x)),
//     &RED,
// ))?;

// root.present()?;

// Ok(())
// not implemented
// }
