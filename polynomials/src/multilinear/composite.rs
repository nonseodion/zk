use core::panic;

use ark_ff::PrimeField;
use crate::multilinear::multilinear::MultiLinear;
use std::ops::{Add, Mul};

// this is limited to work with 
#[derive(Clone, Debug)]
pub enum OP {
  ADD,
  MUL
}

#[derive(Debug, Clone)]
enum OP_ELEMENT<F: PrimeField> {
  Value(F),
  OP(OP),
  Poly(MultiLinear<F>)
}

enum Poly_or_F<F: PrimeField> {
  Value(F),
  Poly(MultiLinear<F>)  
}

fn get_poly<F: PrimeField>(poly: Poly_or_F<F>) -> MultiLinear<F> {
  if let Poly_or_F::Poly(x) = poly {
    return x;
  }
  panic!("Value is not poly");
}

fn get_F<F: PrimeField>(f: Poly_or_F<F>) -> F {
  if let Poly_or_F::Value(x) = f {
    return x;
  }

  panic!("Value is not poly");
}

impl <F: PrimeField> Add for Poly_or_F<F> {
    type Output = Poly_or_F<F>;

    fn add(self, other: Poly_or_F<F>) -> Poly_or_F<F> {
      match self {
        Poly_or_F::Value(value) => Poly_or_F::Value(value + get_F(other)),
        Poly_or_F::Poly(poly) => Poly_or_F::Poly(poly + get_poly(other)),
        _ => panic!("Operator index")
      }
    }
}

impl <F: PrimeField> Mul for Poly_or_F<F> {
    type Output = Poly_or_F<F>;

    fn mul(self, other: Poly_or_F<F>) -> Poly_or_F<F> {
      match self {
        Poly_or_F::Value(value) => Poly_or_F::Value(value * get_F(other)),
        Poly_or_F::Poly(poly) => Poly_or_F::Poly(poly * get_poly(other)),
        _ => panic!("Operator index")
      }
    }
}

fn getPrecedence(op: &OP) -> u8 {
  match op {
    OP::ADD => 0,
    OP::MUL => 1
  }
}

#[derive(Clone, Debug)]
pub struct Composite<F: PrimeField> {
  pub polys: Vec<MultiLinear<F>>,
  pub ops: Vec<OP>
}

impl <F: PrimeField> Composite<F> {
  pub fn new(hypercubes: &Vec<Vec<F>>, ops: Vec<OP>) -> Self {
    if ops.len() + 1 != hypercubes.len() {
      panic!("ops length is not 1 less than hypercubes length");
    }

    let polys = hypercubes.iter().map( |x| MultiLinear::new(x.to_vec())).collect();
    Composite {
      polys,
      ops
    }
  }

  pub fn partialEvaluate(&self, value: F, index: usize) -> Self {
    return Composite {
      polys: self.polys.iter().map(|x| x.partial_evaluate(value, index)).collect(),
      ops: self.ops.to_vec()
    }
  }

  pub fn evaluate(&self, values: &Vec<Option<F>>) -> F {
    let result = self.polys.iter().map( |x| OP_ELEMENT::Value(x.evaluate(values).hypercube[0])).collect();
    if let OP_ELEMENT::Value(x) = shunting_yard_algo(&result, &self.ops).unwrap() {
      return x;
    }
    panic!("Failed to evaluate to a Field Element");
  }

  pub fn reduce (&self) -> MultiLinear<F> {
    let len = self.polys[0].hypercube.len();
    self.polys.iter().for_each(|x| if x.hypercube.len() != len {panic!("Not all the polys have the same length")});

    let result = shunting_yard_algo(&self.polys.iter().map(|x| OP_ELEMENT::Poly(x.clone())).collect(), &self.ops).unwrap();
    if let OP_ELEMENT::Poly(x) = result{
      return x;
    }
    panic!("Failed to evaluate to a multinear");
  }

  // pub fn degree(&self) -> usize {
  //   let mut degree = 1;
  //   let mut increase = 1;
  //   self.ops.iter().for_each(|x| if let OP::MUL = x {degree += 1;});
  //   degree
  // }
}

fn getOp<F: PrimeField> (list: &Vec<OP_ELEMENT<F>>, index: usize) -> OP{
  if let OP_ELEMENT::OP(_operator) = &list[index] {
    return _operator.clone();
  }
  panic!("Operator index is invalid {}", index);
}

fn getValue<F: PrimeField> (list: &Vec<OP_ELEMENT<F>>, index: usize) -> Poly_or_F<F>{
  match &list[index] {
    OP_ELEMENT::Value(value) => return Poly_or_F::Value(value.clone()),
    OP_ELEMENT::Poly(poly) => return Poly_or_F::Poly(poly.clone()),
    _ => panic!("Operator index")
  }
}

// You'll find the shunting yard algorithm here: https://en.wikipedia.org/wiki/Shunting_yard_algorithm#The_algorithm_in_detail.
fn shunting_yard_algo<F: PrimeField>(values: &Vec<OP_ELEMENT<F>>, ops: &Vec<OP>) -> Result<OP_ELEMENT<F>, &'static str> {
  let mut output: Vec<OP_ELEMENT<F>> = vec![];
  let mut operators: Vec<OP> = vec![];
  let mut operator_indexes = vec![];
  if values.len()-1 != ops.len() {
    return Err("operator length is not 1 less than value length");
  }

  // convert experession from infix to postfix e.g 2+2*2 -> [2,2,2,*,+]
  for i in 0..values.len() {
    output.push(values[i].clone());
    if i < ops.len() {
      // when there are operators in operators
      let precedence = getPrecedence(&ops[i]);
      if operators.len() > 0 &&
        precedence <= getPrecedence(&operators[operators.len() - 1]){
          let mut new_operators = operators.clone();
          // push operators to output vector if the next operator has a lesser or equal precedence
          operators.iter().rev().for_each(|x| {
            if precedence <= getPrecedence(x) {
              operator_indexes.push(output.len());
              output.push( OP_ELEMENT::OP(x.clone()));
              new_operators.pop();
            }
          });
          operators = new_operators;
      }

      operators.push(ops[i].clone());
    }
  }

  operators.iter().rev().for_each(|x| {
    operator_indexes.push(output.len());
    output.push( OP_ELEMENT::OP(x.clone()));
  });

  // evaluate postifx expression in ouput vector
  let mut operator_index;
  for i in 0..ops.len() {
    operator_index = operator_indexes[i] - (i * 2);
    let mut right_operand = getValue(&output, operator_index-1);
    let left_operand = getValue(&output, operator_index -2);    
    right_operand =  match getOp(&output, operator_index) {
        OP::ADD => left_operand + right_operand,
        OP::MUL => left_operand * right_operand
    };

    let _right_operand = match right_operand {
      Poly_or_F::Poly(x) => OP_ELEMENT::Poly(x),
      Poly_or_F::Value(x) => OP_ELEMENT::Value(x),
    };

    output[operator_index - 2] = _right_operand;
    output.drain(operator_index-1..operator_index+1);
  }
  
  Ok(output[0].clone())
}

#[cfg(test)]
mod test {
  use std::vec;
  use super::*;
  use ark_bn254::Fq;
use ark_ff::AdditiveGroup;
  use crate::multilinear::composite::shunting_yard_algo;

  // 2*2*3*3*3+3*7*7
  #[test]
  fn test_shunting_yard_algo(){
    let result = shunting_yard_algo(
      &vec![2, 2, 3, 3, 3,  3, 7, 7].iter().map(|x| OP_ELEMENT::Value(Fq::from(*x))).collect(),
      &vec![OP::MUL, OP::MUL, OP::MUL, OP::MUL, OP::ADD, OP::MUL, OP::MUL]
    ).unwrap();

    let mut _result = Fq::ZERO;
    if let OP_ELEMENT::Value(x) = result {
      _result = x;
    };

    assert_eq!(_result, Fq::from(255));
  }

  #[test]
  fn test_evaluate(){
    // (2a + 3b) * (4b + 7ab) + (2ab + 3b + 6a)
    let poly_a = vec![0, 3, 2, 5].iter().map(|x| Fq::from(*x)).collect();
    let poly_b = vec![0, 4, 0, 11].iter().map(|x| Fq::from(*x)).collect();
    let poly_c = vec![0, 3, 6, 11].iter().map(|x| Fq::from(*x)).collect();

    let main_poly = Composite::new(&vec![poly_a, poly_b, poly_c], vec![OP::MUL, OP::ADD]);
    let result = main_poly.evaluate(
      &vec![2, 3].iter().map(|x| Option::Some(Fq::from(*x))).collect()
    );

    assert_eq!(result, Fq::from(735));
  }

  #[test]
  fn test_reduce() {
    let poly_a = vec![0, 3, 2, 5].iter().map(|x| Fq::from(*x)).collect();
    let poly_b = vec![0, 4, 0, 11].iter().map(|x| Fq::from(*x)).collect();
    let poly_c = vec![0, 3, 6, 11].iter().map(|x| Fq::from(*x)).collect();

    let main_poly = Composite::new(&vec![poly_a, poly_b, poly_c], vec![OP::MUL, OP::ADD]);
    let result = main_poly.reduce();

    assert_eq!(
      result.hypercube,
      vec![0, 15, 6, 66].iter().map(|x| Fq::from(*x)).collect::<Vec<Fq>>()
    );
  }
}

