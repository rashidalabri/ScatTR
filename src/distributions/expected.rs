use std::collections::HashMap;

use anyhow::{Context, Ok, Result};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rust_htslib::bam::{self, Read, Record};
use serde::{Deserialize, Serialize};

use super::{normalize_distr, trim_distr};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ExpectedInsertSizeDistr {
    distr: Vec<f64>,
}

impl ExpectedInsertSizeDistr {
    pub fn from_bam(
        aln_rdr: &mut bam::Reader,
        subsample: f64,
        percentile_limit: f64,
        rng: &mut Xoshiro256PlusPlus,
    ) -> Result<Self> {
        // Maps unpaired read names to their positions and lengths
        let mut unpaired_reads: HashMap<Vec<u8>, (i64, usize)> = HashMap::new();
        let mut distr = Vec::new();

        // Since a read pair is only added if the first read is sampled and added
        // to the unpaired_reads map, we need to double the subsample rate
        // so that the expected fraction of read pairs is equal to `subsample`.
        // The derivation assumes large genome size.
        let adjusted_subsample = subsample * 2.0;

        let mut record = Record::new();
        while let Some(result) = aln_rdr.read(&mut record) {
            result.context("Failed to read record from SAM/BAM/CRAM file while estimating insert size distribution")?;

            if record.is_secondary()
                || record.is_supplementary()
                || record.is_unmapped()
                || record.pos() < 0
            {
                continue;
            }

            let qname = record.qname().to_vec();

            let mut remove = false;
            if unpaired_reads.contains_key(&qname) {
                let pos = record.pos();
                let len = record.seq().len();
                let (mate_pos, mate_len) = *unpaired_reads.get(&qname).unwrap();

                remove = true;

                let insert_size = if pos > mate_pos {
                    (pos + len as i64) - mate_pos
                } else {
                    (mate_pos + mate_len as i64) - pos
                };

                assert!(insert_size >= 0, "Encountered negative insert size");

                if distr.len() <= insert_size as usize {
                    distr.resize(insert_size as usize + 1, 0.0);
                }
                distr[insert_size as usize] += 1.0;
            } else if rng.gen_bool(adjusted_subsample) {
                let pos = record.pos();
                let len = record.seq().len();
                unpaired_reads.insert(qname.clone(), (pos, len));
            }

            if remove {
                unpaired_reads.remove(&qname);
            }
        }

        trim_distr(&mut distr, percentile_limit);
        normalize_distr(&mut distr);

        Ok(Self { distr })
    }

    pub fn to_vec(&self) -> Vec<f64> {
        self.distr.clone()
    }
}

// #[derive(Clone, Debug, Serialize, Deserialize)]
// pub struct ExpectedReadDepthDistr {
//     distr: Vec<f64>,
// }

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EstimatedDistr {
    pub distr: Vec<f64>,
}

impl EstimatedDistr {
    pub fn new() -> Self {
        Self { distr: vec![] }
    }

    pub fn from_vec(distr: Vec<f64>) -> Self {
        Self { distr }
    }

    #[inline]
    pub fn add(&mut self, observation: u32) {
        if self.distr.len() <= observation as usize {
            self.distr.resize(observation as usize + 1, 1.0);
        }
        self.distr[observation as usize] += 1.0;
    }

    pub fn add_many(&mut self, observations: &[u32]) {
        for &obs in observations {
            self.add(obs);
        }
    }

    pub fn norm(&mut self) {
        let sum = self.distr.iter().sum::<f64>();
        self.distr.iter_mut().for_each(|x| *x /= sum);
    }

    pub fn mean(&self) -> f64 {
        self.distr
            .iter()
            .enumerate()
            .map(|(depth, prob)| depth as f64 * prob)
            .sum()
    }

    pub fn std_dev(&self) -> f64 {
        self.variance().sqrt()
    }

    fn variance(&self) -> f64 {
        let mean = self.mean();
        self.distr
            .iter()
            .enumerate()
            .map(|(depth, prob)| (depth as f64 - mean).powi(2) * prob)
            .sum()
    }

    pub fn max(&self) -> (usize, f64) {
        self.distr
            .iter()
            .cloned()
            .enumerate()
            .reduce(
                |(i_max, x_max), (i, x)| {
                    if x > x_max {
                        (i, x)
                    } else {
                        (i_max, x_max)
                    }
                },
            )
            .unwrap()
    }
}

impl Default for EstimatedDistr {
    fn default() -> Self {
        Self::new()
    }
}

// implement iterator for EstimatedDistr that does not consume or mutate the distr field
impl<'a> IntoIterator for &'a EstimatedDistr {
    type Item = &'a f64;
    type IntoIter = std::slice::Iter<'a, f64>;

    fn into_iter(self) -> Self::IntoIter {
        self.distr.iter()
    }
}

impl AsRef<Vec<f64>> for EstimatedDistr {
    fn as_ref(&self) -> &Vec<f64> {
        &self.distr
    }
}

impl AsMut<Vec<f64>> for EstimatedDistr {
    fn as_mut(&mut self) -> &mut Vec<f64> {
        &mut self.distr
    }
}

impl From<EstimatedDistr> for Vec<f64> {
    fn from(distr: EstimatedDistr) -> Vec<f64> {
        distr.distr
    }
}

// impl ExpectedReadDepthDistr {
//     pub fn from_bam(
//         aln_rdr: &mut bam::Reader,
//         subsample: f64,
//         percentile_limit: f64,
//         rng: &mut Xoshiro256PlusPlus,
//     ) -> Result<Self> {
//         // Sample the depths of random positions along the alignment file
//         let depths = aln_rdr
//             .pileup()
//             .filter_map(|p| {
//                 if rng.gen_bool(subsample) {
//                     p.ok()
//                 } else {
//                     None
//                 }
//             })
//             .map(|p| p.depth() as usize);

//         // Create a distribution of the sampled depths
//         let mut distr = depths.fold(vec![0.0; 1], |mut hist, depth| {
//             if hist.len() <= depth {
//                 hist.resize(depth + 1, 0.0);
//             }
//             hist[depth] += 1.0;
//             hist
//         });

//         prune_distr(&mut distr, percentile_limit);
//         normalize_distr(&mut distr);

//         Ok(Self { distr })
//     }

//     pub fn to_vec(&self) -> Vec<f64> {
//         self.distr.clone()
//     }
// }
