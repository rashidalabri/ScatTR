use bio::alphabets::dna::revcomp;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::convert::AsRef;
use std::fmt;
use std::ops::RangeFrom;
use std::ops::{Index, IndexMut};
use std::str::FromStr;

#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub struct Sequence {
    pub data: Vec<u8>,
}

impl Sequence {
    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn repeat(&self, n: usize) -> Self {
        let repeated_data = self
            .data
            .iter()
            .cloned()
            .cycle()
            .take(self.data.len() * n)
            .collect();
        Sequence {
            data: repeated_data,
        }
    }

    pub fn rotate_left(&mut self, mid: usize) {
        self.data.rotate_left(mid);
    }

    pub fn revcomp(&self) -> Self {
        let revcomp_data = revcomp(&self.data);
        Sequence { data: revcomp_data }
    }

    pub fn iter(&self) -> SequenceIter {
        SequenceIter {
            seq: self,
            index: 0,
        }
    }
}

#[derive(Clone)]
pub struct SequenceIter<'a> {
    seq: &'a Sequence,
    index: usize,
}

impl<'a> Iterator for SequenceIter<'a> {
    type Item = &'a u8;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.seq.data.len() {
            let item = &self.seq.data[self.index];
            self.index += 1;
            Some(item)
        } else {
            None
        }
    }
}

impl FromIterator<u8> for Sequence {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        let mut seq = Sequence { data: Vec::new() };
        for i in iter {
            seq.data.push(i);
        }
        seq
    }
}

impl Index<RangeFrom<usize>> for Sequence {
    type Output = [u8];

    fn index(&self, index: RangeFrom<usize>) -> &Self::Output {
        &self.data[index]
    }
}

impl IndexMut<RangeFrom<usize>> for Sequence {
    fn index_mut(&mut self, index: RangeFrom<usize>) -> &mut [u8] {
        &mut self.data[index]
    }
}

impl Serialize for Sequence {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let ascii_string = String::from_utf8(self.data.clone()).unwrap();
        serializer.serialize_str(&ascii_string)
    }
}

impl<'de> Deserialize<'de> for Sequence {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let ascii_string = String::deserialize(deserializer)?;
        let data: Vec<u8> = ascii_string.into_bytes();
        Ok(Sequence { data })
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let ascii_string = String::from_utf8(self.data.clone()).unwrap();
        write!(f, "{}", ascii_string)
    }
}

// impl IntoIterator for Sequence {
//     type Item = u8;
//     type IntoIter = std::vec::IntoIter<u8>;

//     fn into_iter(self) -> Self::IntoIter {
//         self.data.into_iter()
//     }
// }

impl FromStr for Sequence {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let data = s.as_bytes().to_vec();
        Ok(Sequence { data })
    }
}

impl From<&str> for Sequence {
    fn from(s: &str) -> Self {
        Sequence {
            data: s.as_bytes().to_vec(),
        }
    }
}

impl AsRef<[u8]> for Sequence {
    fn as_ref(&self) -> &[u8] {
        &self.data
    }
}

impl From<Vec<u8>> for Sequence {
    fn from(data: Vec<u8>) -> Self {
        Sequence { data }
    }
}

impl From<&[u8]> for Sequence {
    fn from(slice: &[u8]) -> Self {
        Sequence {
            data: slice.to_vec(),
        }
    }
}

impl Index<usize> for Sequence {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl Index<std::ops::Range<usize>> for Sequence {
    type Output = [u8];

    fn index(&self, range: std::ops::Range<usize>) -> &Self::Output {
        &self.data[range]
    }
}

impl Index<std::ops::RangeFull> for Sequence {
    type Output = [u8];

    fn index(&self, _: std::ops::RangeFull) -> &Self::Output {
        &self.data[..]
    }
}

impl IndexMut<usize> for Sequence {
    fn index_mut(&mut self, index: usize) -> &mut u8 {
        &mut self.data[index]
    }
}

impl IndexMut<std::ops::Range<usize>> for Sequence {
    fn index_mut(&mut self, range: std::ops::Range<usize>) -> &mut [u8] {
        &mut self.data[range]
    }
}

impl IndexMut<std::ops::RangeFull> for Sequence {
    fn index_mut(&mut self, _: std::ops::RangeFull) -> &mut [u8] {
        &mut self.data[..]
    }
}
