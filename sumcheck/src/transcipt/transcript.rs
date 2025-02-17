use sha3::{Keccak256, Digest};
use std::fmt::Debug;

pub trait HashWrapper {
  fn absorb(&mut self, data: &Vec<u8>);
  fn squeeze(&self) -> Vec<u8>;
  fn clear(&mut self);
}

impl HashWrapper for Keccak256 {
  fn absorb(&mut self, data: &Vec<u8>) {
    self.update(data);
    // TODO try removing clone here
    self.clone().finalize().to_vec();
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
  fn absorb(&mut self, data: &Vec<u8>);
  fn squeeze(&mut self) -> Vec<u8>;
  fn clear(&mut self);
}

#[derive(Debug)]
pub struct Transcript<T: HashWrapper> {
  hasher: T
}

impl<T: HashWrapper> TranscriptTrait<T> for Transcript<T> {

  fn new(hasher: T) -> Self{
    Transcript{
      hasher
    }
  }

  // appends data and returns hash
  fn absorb(&mut self, data: &Vec<u8>){
    self.hasher.absorb(data);
  }

  fn squeeze(&mut self) -> Vec<u8> {
    let hash = self.hasher.squeeze();
    self.hasher.absorb(&hash);
    return hash;
  }

  fn clear(&mut self) {
    self.hasher.clear();
  }
}