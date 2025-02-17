use ark_ff::PrimeField;
use std::fmt::Debug;
use std::iter::repeat_n;
use std::marker::Copy; // TODO: implement copy
use std::clone::Clone;
use std::ops::{Add, Mul};

#[derive(Debug, Clone)]
pub struct MultiLinear<F: PrimeField> {
  pub hypercube: Vec<F>, // result from evaluating the boolean hypercube
}

impl <F: PrimeField> Add for MultiLinear<F> {
    type Output = MultiLinear<F>;

    fn add(self, other: MultiLinear<F>) -> MultiLinear<F> {
      if self.hypercube.len() != other.hypercube.len() {
        panic!("poly hypercube lengths are not equal");
      }

      let mut new_poly_hypercube = vec![];
      for i in 0..self.hypercube.len() {
        new_poly_hypercube.push(self.hypercube[i] + other.hypercube[i]);
      }
      MultiLinear::new(new_poly_hypercube)
    }
}

impl <F: PrimeField> Mul for MultiLinear<F> {
    type Output = MultiLinear<F>;

    fn mul(self, other: MultiLinear<F>) -> MultiLinear<F> {
      if self.hypercube.len() != other.hypercube.len() {
        panic!("poly hypercube lengths are not equal");
      }

      let mut new_poly_hypercube = vec![];
      for i in 0..self.hypercube.len() {
        new_poly_hypercube.push(self.hypercube[i] * other.hypercube[i]);
      }
      MultiLinear::new(new_poly_hypercube)
    }
}

impl <F: PrimeField>MultiLinear<F> {
  pub fn new(hypercube: Vec<F>) -> Self{
    MultiLinear{
      hypercube
    }
  }

  pub fn partial_evaluate(&self, value: F, index: usize) -> Self {
    // let new_hypercube = iter:: ;
    let mut result: MultiLinear<F> = MultiLinear { hypercube: vec![] };
    let skips: usize = 1 << index; // 2 ** index
    let mut new_index = 0;
    let hypercube_len = self.hypercube.len();
    let a = 1 << index + 1; // 2 ** (index + 1)
    let b = (1 << index) - 1; // (2 ** index) - 1
    let mut iteration = 0;

    result.hypercube = repeat_n(0, hypercube_len/2)
      .map(|_| {
        if iteration % a > b {
          new_index = iteration + skips + 1;
        }

        let y0 = self.hypercube[new_index];
        let y1 = self.hypercube[skips + iteration];
        let hypercube_point = y0 + (value * (y1.sub(y0)));

        new_index += 1; 
        iteration += 1;
        return hypercube_point       
      })
      .collect();

    return result;
  }

  pub fn evaluate(&self, values: &Vec<Option<F>>) -> MultiLinear<F> {
      if 2_usize.pow(values.len() as u32) != self.hypercube.len() {
        println!("Polynomial is incorrect");
      }
      let mut hypercube = vec![];
      hypercube.extend( &self.hypercube);
      let mut intermediate_result = MultiLinear::new(hypercube);

      for (i, value) in values.iter().enumerate() {
        intermediate_result = match value {
          Some(_value) => intermediate_result.partial_evaluate(*_value, values.len() - i - 1),
          None => intermediate_result
        }
      }

      intermediate_result
  }
}

pub fn blow_up_right<F: PrimeField>(poly: &MultiLinear<F>, blows: u32) -> MultiLinear<F>{
  let mut new_poly_hypercube = get_blow_up_poly(poly, blows);
  new_poly_hypercube = new_poly_hypercube.iter().enumerate().map( |x| { poly.hypercube[x.0 >> blows]}).collect();
  MultiLinear::new(new_poly_hypercube)
}

pub fn blow_up_left<F: PrimeField>(poly: &MultiLinear<F>, blows: u32) -> MultiLinear<F>{
  let mut new_poly_hypercube = get_blow_up_poly(poly, blows);
  let mask = poly.hypercube.len() - 1;
  new_poly_hypercube = new_poly_hypercube.iter().enumerate().map( |x| poly.hypercube[x.0 & mask]).collect();
  MultiLinear::new(new_poly_hypercube)
}

fn get_blow_up_poly<F: PrimeField>(poly: &MultiLinear<F>, blows: u32) -> Vec<F>{
  if poly.hypercube.len() % 2 != 0 {panic!("poly length is not a multiple of 2")};
  let new_variable_len = 1 << (poly.hypercube.len().trailing_zeros() + blows);
  vec![F::zero(); new_variable_len as usize]
}

pub fn scalar_mul<F: PrimeField>(poly: &MultiLinear<F>, value: F) -> MultiLinear<F> {
  MultiLinear::new(poly.hypercube.iter().map(|x| *x * value).collect())
}

#[cfg(test)]
mod tests {
  use super::{blow_up_left, blow_up_right, scalar_mul, MultiLinear};
    use ark_bn254::Fq;

  #[test]
  fn test_partial_evaluate () {
    let first = MultiLinear{
      hypercube: vec![
        Fq::from(0),
        Fq::from(0),
        Fq::from(0),
        Fq::from(0),
        Fq::from(0),
        Fq::from(4),
        Fq::from(0),
        Fq::from(4),
        Fq::from(0),
        Fq::from(0),
        Fq::from(3),
        Fq::from(3),
        Fq::from(5),
        Fq::from(9),
        Fq::from(8),
        Fq::from(12),
      ]
    };    

    assert_eq!(
      first.evaluate( 
        &vec![Some(Fq::from(4)), None, None, None]).hypercube,
        vec![
          Fq::from(0),
          Fq::from(0),
          Fq::from(12),
          Fq::from(12),
          Fq::from(20),
          Fq::from(24),
          Fq::from(32),
          Fq::from(36)
        ],
    );
  }

  #[test]
  fn test_evaluate() {
    let first = MultiLinear{
      hypercube: vec![
        Fq::from(0),
        Fq::from(0),
        Fq::from(3),
        Fq::from(3),
        Fq::from(0),
        Fq::from(2),
        Fq::from(5),
        Fq::from(7),
      ]
    }; 

    assert_eq!(
      first.evaluate(&vec![
        Some(Fq::from(2)),
        Some(Fq::from(3)),
        Some(Fq::from(2))
      ]).hypercube,
      vec![Fq::from(29)]
    );
  }

  #[test]
  fn test_add() {
    // 2a + 3
    let mut poly_a = MultiLinear::new(vec![3, 5].iter().map(|x| Fq::from(*x)).collect());
    poly_a = blow_up_right(&poly_a, 2);

    // 4b + 2a
    let mut poly_b = MultiLinear::new(vec![0, 4, 2, 6].iter().map(|x| Fq::from(*x)).collect());
    poly_b = blow_up_right(&poly_b, 1);

    // 3c + 2
    let mut poly_c = MultiLinear::new(vec![2, 5].iter().map(|x| Fq::from(*x)).collect());
    poly_c = blow_up_left(&poly_c, 2);

    // 4a + 4b+ 3c + 5
    let result = poly_a + poly_b + poly_c;
    assert_eq!(
      result.hypercube,
      vec![5, 8, 9, 12, 9, 12, 13, 16].iter().map(|x| Fq::from(*x)).collect::<Vec<Fq>>()
    );
  }

  #[test]
  fn test_multiply() {
    // 4b + 2a
    let mut poly_a = MultiLinear::new(vec![0, 4, 2, 6].iter().map(|x| Fq::from(*x)).collect());
    poly_a = blow_up_right(&poly_a, 1);

    // 3c + 2
    let mut poly_b = MultiLinear::new(vec![2, 5].iter().map(|x| Fq::from(*x)).collect());
    poly_b = blow_up_left(&poly_b, 2);

    let result = poly_a * poly_b;
    assert_eq!(
      result.hypercube,
      vec![0, 0, 8, 20, 4, 10, 12, 30].iter().map(|x| Fq::from(*x)).collect::<Vec<Fq>>()
    );
  }

  #[test]
  fn test_scalar_mul() {
    let mut poly = MultiLinear::new(vec![Fq::from(2), Fq::from(3)]);
    poly = scalar_mul(&poly, Fq::from(3));
    assert_eq!(
      poly.hypercube,
      vec![Fq::from(6), Fq::from(9)]
    )
  }
}
