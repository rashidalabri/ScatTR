use crate::{catalog::TandemRepeatLocus, sequence::Sequence};
use anyhow::{Context, Ok, Result};
use bio::io::fasta;
use serde::{Deserialize, Serialize};
use std::fs::File;

#[derive(Debug, Serialize, Deserialize, Clone, Hash, PartialEq, Eq)]
pub struct TandemRepeatReference {
    pub motif: Sequence,
    pub left_flank_seq: Sequence,
    pub right_flank_seq: Sequence,
}

impl TandemRepeatReference {
    pub fn from_fasta(
        reader: &mut fasta::IndexedReader<File>,
        locus: &TandemRepeatLocus,
        left_flank_len: u32,
        right_flank_len: u32,
    ) -> Result<Self> {
        let motif = locus.motif.clone();

        let left_flank_start = (locus.start - left_flank_len as i64).max(0);
        let left_flank_end = locus.start;
        let left_flank_seq =
            fetch_sequence(reader, &locus.contig, left_flank_start, left_flank_end)
                .context("unable to fetch left flank sequence")?;

        let right_flank_seq = fetch_sequence(
            reader,
            &locus.contig,
            locus.end,
            locus.end + right_flank_len as i64,
        )
        .context("unable to fetch right flank sequence")?;

        Ok(Self {
            motif,
            left_flank_seq,
            right_flank_seq,
        })
    }

    /// Returns the left flank sequence followed by the repeat with
    /// `repeat_len` bases
    pub fn left_flank_template(&self, repeat_len: u32) -> Sequence {
        self.left_flank_seq
            .iter()
            .chain(self.motif.iter().cycle().take(repeat_len as usize))
            .copied()
            .collect()
    }

    /// Returns the repeat with `repeat_len` bases followed by the
    /// right flank sequence
    pub fn right_flank_template(&self, repeat_len: u32) -> Sequence {
        let num_full_units = repeat_len / self.motif_len();
        let partial_unit_len = repeat_len % self.motif_len();
        let partial_unit_start = (self.motif_len() - partial_unit_len) as usize;
        let partial_unit = &self.motif[partial_unit_start..];
        partial_unit
            .iter()
            .chain(
                self.motif
                    .iter()
                    .cycle()
                    .take((num_full_units * self.motif_len()) as usize),
            )
            .chain(self.right_flank_seq.iter())
            .copied()
            .collect()
    }

    /// Returns the repeat with `repeat_len` bases
    pub fn repeat_template(&self, repeat_len: u32) -> Sequence {
        self.motif
            .iter()
            .cycle()
            .take(repeat_len as usize)
            .copied()
            .collect()
    }

    pub fn whole_template(&self, repeat_len: u32) -> Sequence {
        self.left_flank_seq
            .iter()
            .chain(self.repeat_template(repeat_len).iter())
            .chain(self.right_flank_seq.iter())
            .copied()
            .collect()
    }

    pub fn motif_len(&self) -> u32 {
        self.motif.len() as u32
    }

    pub fn right_flank_len(&self) -> u32 {
        self.right_flank_seq.len() as u32
    }

    pub fn left_flank_len(&self) -> u32 {
        self.left_flank_seq.len() as u32
    }

    pub fn len(&self, num_repeats: u32) -> u32 {
        (num_repeats * self.motif_len()) + self.left_flank_len() + self.right_flank_len()
    }
}

pub fn fetch_sequence(
    fasta_reader: &mut fasta::IndexedReader<File>,
    contig: &str,
    start: i64,
    end: i64,
) -> Result<Sequence> {
    fasta_reader.fetch(contig, start as u64, end as u64)?;
    let mut seq: Vec<u8> = Vec::new();
    fasta_reader.read(&mut seq)?;
    Ok(seq.into())
}
