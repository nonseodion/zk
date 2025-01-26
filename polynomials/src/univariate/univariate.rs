use std::ops::{Mul, Add};
use std::fmt::Debug;
use ark_ff::PrimeField;

#[derive(Debug)]
pub(crate) struct UnivariatePolynomial  <F: PrimeField> {
    coefficients: Vec<F>
}

impl <F: PrimeField>Mul for UnivariatePolynomial<F>{
    type Output = Self;

    fn mul(self, other: UnivariatePolynomial<F>) -> UnivariatePolynomial<F> {
        let mut result: UnivariatePolynomial<F> = UnivariatePolynomial{ coefficients: vec![]};
        
        for(i, value1) in self.coefficients.iter().enumerate() {
            for (j, value2) in other.coefficients.iter().enumerate() {
                if result.coefficients.len() == i + j {
                    result.coefficients.push(F::from(0));
                }
                result.coefficients[i + j] += value2.mul(value1);
            }
        }

        return result;
    }
}

impl <F: PrimeField>Add for UnivariatePolynomial<F> {
    type Output = Self;

    fn add(self, other: UnivariatePolynomial<F>) -> UnivariatePolynomial<F> {
        let mut result: UnivariatePolynomial<F> = UnivariatePolynomial{ coefficients: vec![]};
        let lhs: UnivariatePolynomial<F>;
        let rhs: UnivariatePolynomial<F>;
        if self.coefficients.len() > other.coefficients.len() {
            lhs = self;
            rhs = other;
        } else {
            lhs = other;
            rhs = self;
        }

        for (i, value) in lhs.coefficients.iter().enumerate() {
            let mut value2 = F::from(0);
            if rhs.coefficients.len() > i {
                value2 = rhs.coefficients[i];
            }
            if result.coefficients.len() == i {
                result.coefficients.push(F::from(0));
            }
            result.coefficients[i] = value.add(value2);
        }

        return result;
    }
}


pub(crate) fn interpolate<F: PrimeField>(points: &Vec<(F, F)>) -> UnivariatePolynomial<F> {

    let mut polynomial = UnivariatePolynomial{ coefficients: vec![]};

    for (i, point) in points.iter().enumerate() {
        let mut numerator1 = UnivariatePolynomial{ coefficients: vec![point.1]};
        let mut denominator = F::from(1);

        for (_i, _point) in points.iter().enumerate() {
            if i == _i { continue;}
            
            let numerator2 = UnivariatePolynomial{coefficients: vec![-F::from(_point.0), F::from(1)]};
            numerator1 = numerator1 * numerator2; 
            denominator *= point.0 - _point.0;
        };

        polynomial = polynomial + (numerator1 * UnivariatePolynomial{coefficients: vec![ F::from(1) / denominator ]});
    }

    return polynomial;
}

#[cfg(test)]
mod tests {
  use ark_bn254::Fq;
  use super::*;

  #[test]
  fn test_interpolate1() {
    let points = vec![
      (Fq::from(0), Fq::from(5)),
      (Fq::from(1), Fq::from(12)),
      (Fq::from(2), Fq::from(33)),
      (Fq::from(3), Fq::from(74)),
    ];

    assert_eq!(
      interpolate(&points).coefficients,
      vec![Fq::from(5), Fq::from(2), Fq::from(4), Fq::from(1)]
    )
  }

  #[test]
  fn test_interpolate2() {
    let points = vec![
      (Fq::from(0), Fq::from(5)),
      (Fq::from(1), Fq::from(4)),
      (Fq::from(2), Fq::from(5)),
    ];

    assert_eq!(
      interpolate(&points).coefficients,
      vec![Fq::from(5), -Fq::from(2), Fq::from(1)]
    )
  }  

}
