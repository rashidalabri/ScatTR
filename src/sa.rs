use clap::Args;
use serde::{Deserialize, Serialize};

use crate::constants::*;

#[derive(Clone, Serialize, Deserialize, Debug, Args)]
pub struct SimulatedAnnealingOptions {
    /// Base error rate
    #[arg(short = 'e', long, default_value_t = DEFAULT_SA_BASE_ERROR_RATE)]
    pub base_error_rate: f64,

    /// Probability of unmapping a read
    #[arg(long, default_value_t = DEFAULT_UNMAPPED_PROB)]
    pub p_unmap: f64,

    /// Use theoretical distributions instead of emperical
    #[arg(long)]
    pub use_theoretical: bool,

    /// Depth distribution mean (if set, the emperical distribution will not be used)
    #[arg(long)]
    pub depth_mean: Option<f64>,

    /// Insert distribution mean (if set, the emperical distribution will not be used)
    #[arg(long)]
    pub insert_mean: Option<f64>,

    /// Insert distribution mean (has to be set if insert_mean is set)
    #[arg(long)]
    pub insert_sd: Option<f64>,

    /// Insert size tolerance
    #[arg(long = "ist", default_value_t = DEFAULT_INSERT_SD_TOLERANCE)]
    pub insert_sd_tolerance: f64,

    /// Fixed read threshold (any read with number of positions below this threshold will only be sampled once)
    #[arg(long = "frt", default_value_t = DEFAULT_FIXED_READ_THRESHOLD)]
    pub fixed_read_threshold: usize,

    /// The initial temperature for simulated annealing of read alignments.
    #[arg(long = "ti", long, default_value_t = DEFAULT_SA_INIT_TEMP)]
    pub sa_init_temp: f64,

    /// The maximum number of iterations for simulated annealing of read alignments.
    #[arg(long = "sf", default_value_t = DEFAULT_SA_STALL_FIXED)]
    pub sa_stall_fixed: u64,

    /// The number of iterations with no accepted solutions after which the SA algorithm stops.
    #[arg(long = "sa", default_value_t = DEFAULT_SA_STALL_ACCEPTED)]
    pub sa_stall_accepted: u64,

    /// The number of iterations after which reannealing is performed. This may help in overcoming local minima.
    #[arg(long = "rf", default_value_t = DEFAULT_SA_REANNEAL_FIXED)]
    pub sa_reanneal_fixed: u64,

    /// The number of iterations that need to pass after the last accepted solution was found for reannealing to be performed.
    #[arg(long = "ra", default_value_t = DEFAULT_SA_REANNEAL_ACCEPTED)]
    pub sa_reanneal_accepted: u64,
}
