use anyhow::{anyhow, Result};
use serde::Serialize;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    sync::{Arc, Mutex, RwLock},
};

use argmin::core::{observers::Observe, State, KV};

use crate::positions::AlignmentPosition;

#[derive(Clone)]
pub struct AlignmentPositionTracker {
    // Map of genotype to best cost seen so far
    best_cost: Arc<RwLock<HashMap<(u32, u32), f64>>>,
    // Map of genotype to best alignment positions seen so far
    #[allow(clippy::type_complexity)]
    best_pos: Arc<Mutex<HashMap<(u32, u32), Vec<AlignmentPosition>>>>,
}

impl AlignmentPositionTracker {
    pub fn new() -> Self {
        Self {
            best_cost: Arc::new(RwLock::new(HashMap::new())),
            best_pos: Arc::new(Mutex::new(HashMap::new())),
        }
    }

    pub fn update(&self, genotype: Vec<u32>, cost: f64, pos: &[AlignmentPosition]) {
        let genotype = (genotype[0], genotype[1]);

        let is_new_best = {
            let best_cost = self.best_cost.read().unwrap();
            best_cost.get(&genotype).map_or(true, |&best| cost < best)
        };

        if is_new_best {
            let mut best_cost = self.best_cost.write().unwrap();
            let mut best_pos = self.best_pos.lock().unwrap();
            best_cost.insert(genotype, cost);
            best_pos.insert(genotype, pos.to_vec());
        }
    }

    pub fn get_best(&self, genotype: Vec<u32>) -> Option<(f64, Vec<AlignmentPosition>)> {
        let genotype = (genotype[0], genotype[1]);
        let best_cost = self.best_cost.read().unwrap();
        let best_pos = self.best_pos.lock().unwrap();
        best_cost
            .get(&genotype)
            .map(|&cost| (cost, best_pos[&genotype].clone()))
    }

    pub fn contains(&self, genotype: Vec<u32>) -> bool {
        let genotype = (genotype[0], genotype[1]);
        let best_cost = self.best_cost.read().unwrap();
        best_cost.contains_key(&genotype)
    }
}

impl Default for AlignmentPositionTracker {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone)]
pub struct SharedFileObserver {
    writer: Arc<Mutex<BufWriter<File>>>,
}

impl SharedFileObserver {
    pub fn new(f: File) -> Self {
        let mut writer = BufWriter::new(f);
        writeln!(writer, "iter\tparam\tcost\ttime").unwrap();
        Self {
            writer: Arc::new(Mutex::new(writer)),
        }
    }

    pub fn flush(&self) -> Result<()> {
        let mut writer = self
            .writer
            .lock()
            .map_err(|_| anyhow!("could not lock log file"))?;
        writer.flush()?;
        Ok(())
    }
}

impl<I> Observe<I> for SharedFileObserver
where
    I: State,
    <I as State>::Param: Serialize,
{
    fn observe_iter(&mut self, state: &I, _kv: &KV) -> Result<(), argmin::core::Error> {
        let mut writer = self
            .writer
            .lock()
            .map_err(|_| anyhow!("could not lock log file"))?;

        let param = state
            .get_param()
            .ok_or(anyhow!("could not get param for log"))?;

        let _param = serde_json::to_string(&param)
            .map_err(|_| anyhow!("could not serialize param for log"))?;

        writeln!(
            writer,
            "{},{},{:?}",
            state.get_iter(),
            state.get_cost(),
            state.get_time(),
        )?;
        Ok(())
    }
}
