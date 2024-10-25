use anyhow::{anyhow, Result};
use bio::io::fasta;
use rust_htslib::bam::{self, Read};
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use std::{
    collections::HashMap,
    fs::File,
    path::{Path, PathBuf},
};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SampleStats {
    pub read_length: u64,
    pub depth_mean: f64,
    pub depth_sd: f64,
    pub insert_size_mean: f64,
    pub insert_size_sd: f64,
    depth_distr: Vec<f64>,
    insert_distr: Vec<f64>,
}

impl SampleStats {
    pub fn new(
        read_length: u64,
        depth_mean: f64,
        depth_sd: f64,
        insert_size_mean: f64,
        insert_size_sd: f64,
        depth_distr: Vec<f64>,
        insert_distr: Vec<f64>,
    ) -> Self {
        Self {
            read_length,
            depth_mean,
            depth_sd,
            insert_size_mean,
            insert_size_sd,
            depth_distr,
            insert_distr,
        }
    }

    pub fn depth_distr(&self) -> EmpiricalDepthDistribution {
        let halved_distr = halve_distribution(self.depth_distr.clone());
        EmpiricalDepthDistribution::new(halved_distr)
    }

    pub fn insert_distr(&self) -> EmpiricalInsertSizeDistribution {
        EmpiricalInsertSizeDistribution::new(self.insert_distr.clone())
    }
}

use crate::{
    distributions::{
        halve_distribution, DepthDistribution, EmpiricalDepthDistribution,
        EmpiricalInsertSizeDistribution, InsertSizeDistribution, Normal, Poisson,
    },
    genotype::{GenotypeProblemDefinition, LocusGenotypeResult},
};

#[allow(clippy::type_complexity)]
pub fn load_aln_inputs(
    definitions_path: &PathBuf,
    stats_path: &PathBuf,
    use_theoretical: bool,
    depth_mean: Option<f64>,
    insert_mean: Option<f64>,
    insert_sd: Option<f64>,
) -> Result<(
    Vec<GenotypeProblemDefinition>,
    Box<dyn DepthDistribution>,
    Box<dyn InsertSizeDistribution>,
)> {
    // Load inputs
    let (depth_distr, insert_distr): (Box<dyn DepthDistribution>, Box<dyn InsertSizeDistribution>) =
        match (depth_mean, insert_mean, insert_sd) {
            (Some(depth_mean), Some(insert_mean), Some(insert_sd)) => {
                debug!("Loading sample stats from arguments");
                debug!("Using theoretical distributions");
                (
                    // Divide by two because the distribution is for a single haplotype
                    Box::new(Poisson::new(depth_mean / 2.0)),
                    Box::new(Normal::new(insert_mean, insert_sd)),
                )
            }
            (None, None, None) => {
                debug!("Loading sample stats from file");
                let file = File::open(stats_path)?;
                let sample_stats: SampleStats = serde_json::from_reader(file)?;

                if use_theoretical {
                    debug!("Using theoretical distributions");
                    (
                        // Divide by two because we want the distribution to be for a single haplotype
                        Box::new(Poisson::new(sample_stats.depth_mean / 2.0)),
                        Box::new(Normal::new(
                            sample_stats.insert_size_mean,
                            sample_stats.insert_size_sd,
                        )),
                    )
                } else {
                    (
                        Box::new(sample_stats.depth_distr()),
                        Box::new(sample_stats.insert_distr()),
                    )
                }
            }
            _ => {
                return Err(anyhow!(
                    "Options --depth-mean, --insert-mean, and --insert-sd must be set together"
                ))
            }
        };

    info!(
        "Loaded haplotype depth distribution with mean = {}",
        depth_distr.mean()
    );
    info!(
        "Loaded insert size distribution with mean = {} and sd = {}",
        insert_distr.mean(),
        insert_distr.sd()
    );

    info!("Loading problem definitions");
    let file = File::open(definitions_path)?;
    let definitions: Vec<GenotypeProblemDefinition> = serde_json::from_reader(file)?;

    Ok((definitions, depth_distr, insert_distr))
}

pub trait Utf8String {
    fn to_string(&self) -> String;
}

impl Utf8String for &[u8] {
    fn to_string(&self) -> String {
        String::from_utf8_lossy(self).to_string()
    }
}

impl Utf8String for Vec<u8> {
    fn to_string(&self) -> String {
        String::from_utf8_lossy(self).to_string()
    }
}

pub fn get_contigs_from_fasta(path: &PathBuf) -> Vec<String> {
    let reader = fasta::Reader::from_file(path).unwrap();
    reader
        .records()
        .map(|rec| rec.unwrap().id().to_string())
        .collect()
}

pub fn get_contig_tid_and_len(bam: &bam::Reader, contig: &str) -> Result<(u32, usize)> {
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
    let tid = header_view
        .tid(contig.as_bytes())
        .ok_or(anyhow!("Could not find contig in BAM file."))?;
    let len = header_view
        .target_len(tid)
        .ok_or(anyhow!("Could not find contig in BAM file."))?;
    Ok((tid, len as usize))
}

pub fn get_contig_tid_and_len_index(
    bam: &bam::IndexedReader,
    contig: &str,
) -> Result<(u32, usize)> {
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
    let tid = header_view
        .tid(contig.as_bytes())
        .ok_or(anyhow!("Could not find contig in BAM file."))?;
    let len = header_view
        .target_len(tid)
        .ok_or(anyhow!("Could not find contig in BAM file."))?;
    Ok((tid, len as usize))
}

pub fn fetch_region_from_fasta(
    fasta: &mut fasta::IndexedReader<File>,
    contig: &str,
    start: i64,
    end: i64,
) -> Result<Vec<u8>> {
    fasta.fetch(contig, start as u64, end as u64)?;
    let mut seq: Vec<u8> = Vec::new();
    fasta.read(&mut seq)?;
    Ok(seq)
}

pub fn serialize_to_json<T, P>(
    data: &HashMap<String, LocusGenotypeResult>,
    file_path: P,
) -> Result<()>
where
    T: Serialize,
    P: AsRef<Path>,
{
    let file = File::create(file_path)?;
    serde_json::to_writer_pretty(file, data)?;
    Ok(())
}

pub fn deserialize_from_json<T, P>(file_path: P) -> Result<Vec<T>>
where
    T: DeserializeOwned,
    P: AsRef<Path>,
{
    let file = File::open(file_path)?;
    let data = serde_json::from_reader(file)?;
    Ok(data)
}

pub fn genotype_bounds(len_lower_bound: u32, len_upper_bound: u32, motif_len: usize) -> (f64, f64) {
    let num_repeats_lower = (len_lower_bound as f64 / motif_len as f64).floor();
    let num_repeats_upper = (len_upper_bound as f64 / motif_len as f64).ceil();
    (num_repeats_lower, num_repeats_upper)
}
