use ark_ff::PrimeField;
use std::fmt::Debug;
use std::io::Read;
use std::iter::repeat_n;
use std::marker::Copy; // TODO: implement copy
use std::clone::Clone;


#[derive(Debug, Clone)]
pub struct MultiLinear<F: PrimeField> {
  pub hypercube: Vec<F>, // result from evaluating the boolean hypercube
}

impl <F: PrimeField>MultiLinear<F> {
  pub fn new(hypercube: Vec<F>) -> Self{
    MultiLinear{
      hypercube
    }
  }

  pub fn partial_evaluate(&self, value: F, index: u32) -> Self {
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
          Some(_value) => intermediate_result.partial_evaluate(*_value, (values.len() - i - 1) as u32),
          None => intermediate_result
        }
      }

      intermediate_result
  }
}




#[cfg(test)]
mod tests {
    use super::MultiLinear;
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
}
