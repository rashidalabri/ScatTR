use crate::{
    constants::*,
    extract::LocusId,
    genotype::{GenotypeInterval, GenotypeProblem, LocusGenotypeResult, TandemRepeatGenotype},
    sa::SimulatedAnnealingOptions,
    util::load_aln_inputs,
};
use argmin::core::{observers::ObserverMode, CostFunction, Executor, State};
use argmin_observer_slog::SlogLogger;
use clap::{arg, Parser};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::path::PathBuf;
use std::{
    collections::BTreeMap,
    sync::{Arc, Mutex},
};
use std::{fs::File, path::Path};

/// Estimates the allele lengths for each tandem repeat locus
#[derive(Parser)]
pub struct AlignCommandArgs {
    /// Stats JSON file
    #[arg(long)]
    pub stats: Option<PathBuf>,

    /// Defs JSON file
    #[arg(long)]
    pub defs: Option<PathBuf>,

    /// Number of repeats start
    pub genotype_start: u32,

    /// Number of repeats end
    pub genotype_end: u32,

    /// Number of repeats step size
    pub genotype_step: u32,

    /// Number of trials to run for each genotype
    #[arg(short = 'n', long, default_value_t = 1)]
    pub n_trials: u32,

    /// Output file
    #[arg(short = 'a', long)]
    pub output: Option<PathBuf>,

    /// Simulated annealing options
    #[command(flatten)]
    pub sa_opts: SimulatedAnnealingOptions,

    /// Output alignment positions
    #[arg(short = 'o', long)]
    pub output_positions: bool,

    /// Random seed
    #[arg(short = 's', long, default_value_t = 1)]
    pub seed: u64,

    /// Enables logging
    #[arg(long, default_value_t = false)]
    pub log: bool,
}

pub fn run_align_command(args: &AlignCommandArgs, output_prefix: &Path) -> anyhow::Result<()> {
    // Find the paths to the input files
    let definitions_path = if let Some(definitions_path) = &args.defs {
        definitions_path
    } else {
        &output_prefix.with_extension(OUT_SUFFIX_DEFS)
    };

    let stats_path = if let Some(stats_path) = &args.stats {
        stats_path
    } else {
        &output_prefix.with_extension(OUT_SUFFIX_STATS)
    };

    // Load inputs
    let (definitions, depth_distr, insert_distr) = load_aln_inputs(
        args.defs.as_ref().unwrap_or(definitions_path),
        args.stats.as_ref().unwrap_or(stats_path),
        args.sa_opts.use_theoretical,
        args.sa_opts.depth_mean,
        args.sa_opts.insert_mean,
        args.sa_opts.insert_sd,
    )?;

    // Create a new RNG
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(args.seed);

    // Results
    let gt_results: Arc<Mutex<BTreeMap<LocusId, Vec<LocusGenotypeResult>>>> =
        Arc::new(Mutex::new(BTreeMap::new()));

    for aln_problem_def in definitions.iter() {
        info!(
            "[{}] Starting alignment for locus",
            aln_problem_def.locus.id
        );

        if aln_problem_def.reads.is_empty() {
            warn!(
                "[{}] Locus does not have any reads. Skipping...",
                aln_problem_def.locus.id
            );
            continue;
        }

        let genotypes =
            (args.genotype_start..=args.genotype_end).step_by(args.genotype_step as usize);

        let genotype_rngs = genotypes
            .clone()
            .map(|_| Xoshiro256PlusPlus::from_rng(&mut rng).unwrap())
            .collect::<Vec<_>>();

        genotypes
            .zip(genotype_rngs)
            .par_bridge()
            .for_each(|(genotype, mut genotype_rng)| {
                for _ in 0..args.n_trials {
                    let trial_rng = Xoshiro256PlusPlus::from_rng(&mut genotype_rng).unwrap();
                    let gt_problem = GenotypeProblem::new(
                        aln_problem_def,
                        args.sa_opts.base_error_rate,
                        &*depth_distr,
                        &*insert_distr,
                        &args.sa_opts,
                        None,
                        trial_rng,
                    );

                    let aln_problem = gt_problem
                        .build_aln_problem(TandemRepeatGenotype(genotype, genotype))
                        .unwrap();

                    let solver = gt_problem.build_solver(&aln_problem).unwrap();

                    let slog_logger = SlogLogger::term_noblock();

                    let initial_pos = aln_problem.initial_positions();
                    let initial_cost = aln_problem.cost(&initial_pos).unwrap();

                    let executor = Executor::new(aln_problem.clone(), solver)
                        .configure(|state| {
                            state
                                .param(initial_pos)
                                .max_iters(args.sa_opts.sa_stall_fixed)
                        })
                        .add_observer(slog_logger, ObserverMode::Never);

                    let result = executor.run().unwrap();
                    let positions = result.state().get_best_param().unwrap().clone();
                    let cost = result.state().get_best_cost();

                    let num_unmapped = positions.iter().filter(|p| !p.is_mapped()).count();

                    let depths = aln_problem.get_depths(&positions);
                    let depth_mean = depths.iter().sum::<u32>() as f64 / depths.len() as f64;
                    let depth_var = depths
                        .iter()
                        .map(|d| (*d as f64 - depth_mean).powi(2))
                        .sum::<f64>() / depths.len() as f64;
                    let depth_sd = depth_var.sqrt();

                    let insert_sizes = aln_problem.get_insert_sizes(&positions);
                    let insert_mean =
                        insert_sizes.iter().sum::<u32>() as f64 / insert_sizes.len() as f64;
                    let insert_var = insert_sizes
                        .iter()
                        .map(|d| (*d as f64 - insert_mean).powi(2))
                        .sum::<f64>() / insert_sizes.len() as f64;
                    let insert_sd = insert_var.sqrt();

                    info!(
                        "[{}] Genotype: {}, Initial Cost: {}, Best Cost: {}, Num Unmapped: {}, Depth Mean: {}, Depth SD: {}, Insert Mean: {}, Insert SD: {}",
                        aln_problem_def.locus.id, genotype, initial_cost, cost, num_unmapped, depth_mean, depth_sd, insert_mean, insert_sd
                    );

                    // Print idx of unmapped reads
                    info!("[{}] Unmapped reads: ", aln_problem_def.locus.id);
                    for (i, p) in positions.iter().enumerate() {
                        if !p.is_mapped() {
                            print!("{}, ", i);
                        }
                    }
                    println!();

                    // info!("genotype: {:?}, pos: {:?}", genotype, pos);
                    // let depths = aln_problem.get_depths(&positions);
                    // info!("genotype: {}, depths: {:?}", genotype, depths);

                    // for (i, j) in aln_problem.pairs.iter() {
                    //     let pos_i = pos[*i];
                    //     let pos_j = pos[*j];
                    //     // if pos_i == AlignmentPosition::Un && pos_j != AlignmentPosition::Unaligned {
                    //     //     print!("({}: {:?}, {}: {:?}) ", i, pos_i, j, pos_j);
                    //     // }
                    //     // if pos_j == AlignmentPosition::Unaligned
                    //     //     && pos_i != AlignmentPosition::Unaligned
                    //     // {
                    //     print!("({}: {:?}, {}: {:?}) ", i, pos_i, j, pos_j);
                    //     // }
                    // }
                    // println!();

                    // let depth_cost = aln_problem.depth_cost(&pos);
                    // let insert_cost = aln_problem.insert_cost(&pos);
                    // let edit_cost = aln_problem.edit_cost(&pos);
                    // let entropy_cost = aln_problem.entropy_cost(&pos);
                    // let num_unmapped_cost = aln_problem.num_unmapped_cost(&pos);

                    let locus_gt_result = LocusGenotypeResult {
                        reference_region: format!("{}", aln_problem_def.locus),
                        repeat_unit: aln_problem_def.locus.motif.clone(),
                        genotype: TandemRepeatGenotype(genotype, genotype),
                        genotype_conf_int: GenotypeInterval(
                            (genotype, genotype),
                            (genotype, genotype),
                        ),
                        cost,
                        alignment: if args.output_positions {
                            Some(positions)
                        } else {
                            None
                        },
                    };

                    gt_results
                        .lock()
                        .unwrap()
                        .entry(aln_problem_def.locus.id.clone())
                        .or_default()
                        .push(locus_gt_result);
                }
            });
    }

    // Sort gt_results by locus id and genotype
    for locus_results in gt_results.lock().unwrap().values_mut() {
        locus_results.sort_by(|a, b| a.genotype.0.cmp(&b.genotype.0));
    }

    // Output alignment results
    let gt_results_path = args
        .output
        .clone()
        .unwrap_or_else(|| output_prefix.with_extension(OUT_SUFFIX_ALIGN));
    let file = File::create(gt_results_path)?;
    serde_json::to_writer_pretty(file, &gt_results)?;

    Ok(())
}
