use sha3::{Keccak256, Digest};
use std::fmt::Debug;

pub trait HashWrapper {
  fn absorb(&mut self, data: &Vec<u8>) -> Vec<u8>;
  fn squeeze(&self) -> Vec<u8>;
  fn clear(&mut self);
}

impl HashWrapper for Keccak256 {
  fn absorb(&mut self, data: &Vec<u8>) -> Vec<u8> {
    self.update(data);
    // TODO try removing clone here
    self.clone().finalize().to_vec()
  }

  fn squeeze(&self) -> Vec<u8> {
      //TODO try removing clone here
      self.clone().finalize().to_vec()
  }

  fn clear(&mut self) {
    self.reset();
  }
}

pub trait TranscriptTrait<T: HashWrapper> {
  fn new(hasher: T) -> Self;
  fn append(&mut self, data: &Vec<u8>) -> Vec<u8>;
  fn hash(self) -> Vec<u8>;
  fn clear(&mut self);
}

#[derive(Debug)]
pub(crate) struct Transcript<T: HashWrapper> {
  hasher: T
}

impl<T: HashWrapper> TranscriptTrait<T> for Transcript<T> {

  fn new(hasher: T) -> Self{
    Transcript{
      hasher
    }
  }

  // appends data and returns hash
  fn append(&mut self, data: &Vec<u8>) -> Vec<u8> {
    self.hasher.absorb(data)
  }

  fn hash(self) -> Vec<u8> {
    self.hasher.squeeze()
  }

  fn clear(&mut self) {
    self.hasher.clear();
  }
}