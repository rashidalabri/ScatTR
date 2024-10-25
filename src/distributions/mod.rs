// pub mod dynamic;
pub mod expected;

use rand::{prelude::Distribution as _, Rng};
use statrs::{consts, function::factorial::ln_factorial};

pub trait DistributionProb {
    fn ln_prob(&self, x: u64) -> f64;
    fn mean(&self) -> f64;
    fn sd(&self) -> f64;
}

pub trait DistributionSample {
    fn sample(&self, rng: &mut dyn rand::RngCore) -> u64;
}

pub trait DistributionProbSample: DistributionProb + DistributionSample {}

pub trait DepthDistribution: DistributionProb + Send + Sync {}
pub trait InsertSizeDistribution: DistributionProbSample + Send + Sync {}

pub struct EmpiricalDepthDistribution {
    depth_distr: Vec<f64>,
    mean: f64,
    sd: f64,
    fallback_poisson: Poisson,
}

impl EmpiricalDepthDistribution {
    pub fn new(depth_distr: Vec<f64>) -> Self {
        let (mean, sd) = distr_mean_sd(&depth_distr);
        let fallback_poisson = Poisson::new(mean);
        Self {
            depth_distr,
            mean,
            sd,
            fallback_poisson,
        }
    }
}

impl DistributionProb for EmpiricalDepthDistribution {
    fn ln_prob(&self, x: u64) -> f64 {
        if x < self.depth_distr.len() as u64 {
            self.depth_distr[x as usize].ln()
        } else {
            self.fallback_poisson.ln_prob(x)
        }
    }

    fn mean(&self) -> f64 {
        self.mean
    }

    fn sd(&self) -> f64 {
        self.sd
    }
}

unsafe impl Send for EmpiricalDepthDistribution {}
unsafe impl Sync for EmpiricalDepthDistribution {}
impl DepthDistribution for EmpiricalDepthDistribution {}

pub struct EmpiricalInsertSizeDistribution {
    insert_distr: Vec<f64>,
    mean: f64,
    sd: f64,
    fallback_normal: Normal,
}

impl EmpiricalInsertSizeDistribution {
    pub fn new(insert_distr: Vec<f64>) -> Self {
        let (mean, sd) = distr_mean_sd(&insert_distr);
        let fallback_normal = Normal::new(mean, sd);
        Self {
            insert_distr,
            mean,
            sd,
            fallback_normal,
        }
    }
}

impl DistributionProb for EmpiricalInsertSizeDistribution {
    fn ln_prob(&self, x: u64) -> f64 {
        if x < self.insert_distr.len() as u64 {
            self.insert_distr[x as usize].ln()
        } else {
            self.fallback_normal.ln_prob(x)
        }
    }

    fn mean(&self) -> f64 {
        self.mean
    }

    fn sd(&self) -> f64 {
        self.sd
    }
}

impl DistributionSample for EmpiricalInsertSizeDistribution {
    fn sample(&self, rng: &mut dyn rand::RngCore) -> u64 {
        let mut cum_sum = 0.0;
        let sample = rng.gen::<f64>();
        for (i, &p) in self.insert_distr.iter().enumerate() {
            cum_sum += p;
            if cum_sum >= sample {
                return i as u64;
            }
        }
        self.insert_distr.len() as u64 - 1
    }
}

impl DistributionProbSample for EmpiricalInsertSizeDistribution {}
impl InsertSizeDistribution for EmpiricalInsertSizeDistribution {}

unsafe impl Send for EmpiricalInsertSizeDistribution {}
unsafe impl Sync for EmpiricalInsertSizeDistribution {}

pub struct Poisson {
    lambda: f64,
}

impl Poisson {
    pub fn new(lambda: f64) -> Self {
        Self { lambda }
    }
}

impl DistributionProb for Poisson {
    /// Poisson distribution probability mass function
    /// Source: https://github.com/statrs-dev/statrs/blob/master/src/distribution/poisson.rs
    fn ln_prob(&self, x: u64) -> f64 {
        -self.lambda + x as f64 * self.lambda.ln() - ln_factorial(x)
    }

    fn mean(&self) -> f64 {
        self.lambda
    }

    fn sd(&self) -> f64 {
        self.lambda.sqrt()
    }
}

unsafe impl Send for Poisson {}
unsafe impl Sync for Poisson {}
impl DepthDistribution for Poisson {}

pub struct Normal {
    mean: f64,
    std_dev: f64,
}

impl Normal {
    pub fn new(mean: f64, std_dev: f64) -> Self {
        Self { mean, std_dev }
    }
}

impl DistributionProb for Normal {
    /// Normal distribution probability density function
    /// Source: https://github.com/statrs-dev/statrs/blob/master/src/distribution/normal.rs
    fn ln_prob(&self, x: u64) -> f64 {
        let d = (x as f64 - self.mean) / self.std_dev;
        (-0.5 * d * d) - consts::LN_SQRT_2PI - self.std_dev.ln()
    }

    fn mean(&self) -> f64 {
        self.mean
    }

    fn sd(&self) -> f64 {
        self.std_dev
    }
}

impl DistributionSample for Normal {
    fn sample(&self, rng: &mut dyn rand::RngCore) -> u64 {
        let normal = statrs::distribution::Normal::new(self.mean, self.std_dev).unwrap();
        normal.sample(rng) as u64
    }
}

impl DistributionProbSample for Normal {}

unsafe impl Send for Normal {}
unsafe impl Sync for Normal {}
impl InsertSizeDistribution for Normal {}

fn distr_mean_sd(distr: &[f64]) -> (f64, f64) {
    let mean = distr_mean(distr);
    let sd = distr_sd(distr, mean);
    (mean, sd)
}

fn distr_mean(distr: &[f64]) -> f64 {
    let mut sum = 0.0;
    for (i, &p) in distr.iter().enumerate() {
        sum += i as f64 * p;
    }
    sum
}

fn distr_sd(distr: &[f64], mean: f64) -> f64 {
    let mut sum = 0.0;
    for (i, &p) in distr.iter().enumerate() {
        sum += (i as f64 - mean).powi(2) * p;
    }
    sum.sqrt()
}

pub fn normalize_distr(distr: &mut [f64]) {
    let sum: f64 = distr.iter().sum();
    for x in distr.iter_mut() {
        *x /= sum;
    }
}

pub fn trim_distr(distr: &mut Vec<f64>, pct: f64) {
    let cum_sum: Vec<f64> = distr
        .iter()
        .scan(0.0, |acc, &x| {
            *acc += x;
            Some(*acc)
        })
        .collect();
    let threshold = pct * cum_sum.last().unwrap();
    let prune_to_idx = cum_sum
        .iter()
        .position(|&x| x >= threshold)
        .unwrap_or(distr.len());
    distr.resize(prune_to_idx, 0.0);
}

pub fn weighted_mean(values: &[u32], weights: &[f64]) -> f64 {
    let values = values.iter();
    let weights = weights.iter();

    let mut sum_weighted_values = 0.0;
    let mut sum_weights = 0.0;

    for (value, weight) in values.zip(weights) {
        sum_weighted_values += *value as f64 * weight;
        sum_weights += weight;
    }

    if sum_weights == 0.0 {
        return f64::NAN;
    }

    sum_weighted_values / sum_weights
}

pub fn weighted_percentile(observations: &[u32], weights: &[f64], percentile: f64) -> Option<f64> {
    if observations.is_empty() || weights.is_empty() || observations.len() != weights.len() {
        return None;
    }

    let total_weight: f64 = weights.iter().sum();
    if total_weight == 0.0 {
        return None;
    }

    let mut sorted_data: Vec<(u32, f64)> = observations
        .iter()
        .zip(weights.iter())
        .map(|(&obs, &weight)| (obs, weight))
        .collect();
    sorted_data.sort_by(|a, b| a.0.cmp(&b.0));

    let mut cum_weight = 0.0;
    let target_weight = percentile / 100.0 * total_weight;

    for (obs, weight) in sorted_data.iter() {
        cum_weight += weight;
        if cum_weight >= target_weight {
            return Some(*obs as f64);
        }
    }

    Some(sorted_data.last().unwrap().0 as f64)
}

pub fn halve_distribution(d: Vec<f64>) -> Vec<f64> {
    let mut halved_d = Vec::new();
    let len = d.len();

    for i in (0..len).step_by(2) {
        if i + 1 < len {
            halved_d.push(d[i] + d[i + 1]);
        } else {
            // If the length is odd, append the last element as is
            halved_d.push(d[i]);
        }
    }

    normalize_distr(&mut halved_d);

    halved_d
}
