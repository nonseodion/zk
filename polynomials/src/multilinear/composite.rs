use core::panic;

use ark_ff::PrimeField;
use crate::multilinear::multilinear::MultiLinear;

// this is limited to work with 
#[derive(Clone, Debug)]
pub(crate) enum OP {
  ADD,
  MUL
}

#[derive(Debug)]
enum SY_OUTPUT_ELEMENT<F: PrimeField> {
  Value(F),
  OP(OP)
}

fn getPrecedence(op: &OP) -> u8 {
  match op {
    OP::ADD => 0,
    OP::MUL => 1
  }
}

struct Composite<F: PrimeField> {
  polys: Vec<MultiLinear<F>>,
  ops: Vec<OP>
}

impl <F: PrimeField> Composite<F> {
  fn new(hypercubes: &Vec<Vec<F>>, ops: Vec<OP>) -> Self {
    if ops.len() + 1 != hypercubes.len() {
      panic!("ops length is not 1 less than hypercubes length");
    }

    let polys = hypercubes.iter().map( |x| MultiLinear::new(x.to_vec())).collect();
    Composite {
      polys,
      ops
    }
  }

  fn partialEvaluate(&self, value: F, index: u32) -> Self {
    return Composite {
      polys: self.polys.iter().map(|x| x.partial_evaluate(value, index)).collect(),
      ops: self.ops.to_vec()
    }
  }

  fn evaluate(&self, values: &Vec<Option<F>>) -> F {
    let result = self.polys.iter().map( |x| x.evaluate(values).hypercube[0]).collect();
    shunting_yard_algo(&result, &self.ops).unwrap()
  }
}

fn getOp<F: PrimeField> (list: &Vec<SY_OUTPUT_ELEMENT<F>>, index: usize) -> OP{
  if let SY_OUTPUT_ELEMENT::OP(_operator) = &list[index] {
    return _operator.clone();
  }
  panic!("Operator index is invalid {}", index);
}

fn getValue<F: PrimeField> (list: &Vec<SY_OUTPUT_ELEMENT<F>>, index: usize) -> F{
  if let SY_OUTPUT_ELEMENT::Value(value) = &list[index] {
    return value.clone();
  }
  panic!("Value index is invalid {}", index);
}

// You'll find the shunting yard algorithm here: https://en.wikipedia.org/wiki/Shunting_yard_algorithm#The_algorithm_in_detail.
fn shunting_yard_algo<F: PrimeField>(values: &Vec<F>, ops: &Vec<OP>) -> Result<F, &'static str> {
  let mut output: Vec<SY_OUTPUT_ELEMENT<F>> = vec![];
  let mut operators: Vec<OP> = vec![];
  let mut operator_indexes = vec![];
  if values.len()-1 != ops.len() {
    return Err("operator length is not 1 less than value length");
  }

  // convert experession from infix to postfix e.g 2+2*2 -> [2,2,2,*,+]
  for i in 0..values.len() {
    output.push(SY_OUTPUT_ELEMENT::Value(values[i]));
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
              output.push( SY_OUTPUT_ELEMENT::OP(x.clone()));
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
    output.push( SY_OUTPUT_ELEMENT::OP(x.clone()));
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

    output[operator_index - 2] = SY_OUTPUT_ELEMENT::Value(right_operand);
    output.drain(operator_index-1..operator_index+1);
  }
  
  Ok(getValue(&output, 0))
}

#[cfg(test)]
mod test {
  use std::vec;
  use super::*;

use ark_bn254::Fq;

use crate::multilinear::composite::shunting_yard_algo;
  // 2*2*3*3*3+3*7*7
  #[test]
  fn test_shunting_yard_algo(){
    let result = shunting_yard_algo(
      &vec![Fq::from(2), Fq::from(2), Fq::from(3), Fq::from(3), Fq::from(3),  Fq::from(3), Fq::from(7), Fq::from(7)],
      &vec![OP::MUL, OP::MUL, OP::MUL, OP::MUL, OP::ADD, OP::MUL, OP::MUL]
    );

    assert_eq!(result.unwrap(), Fq::from(255));
  }

  #[test]
  fn test_composite(){
    // (2a + 3b) * (4b + 7ab) + (2ab + 3b + 6a)
    let poly_a = vec![0, 3, 2, 5].iter().map(|x| Fq::from(x.clone())).collect();
    let poly_b = vec![0, 4, 0, 11].iter().map(|x| Fq::from(x.clone())).collect();
    let poly_c = vec![0, 3, 6, 11].iter().map(|x| Fq::from(x.clone())).collect();

    let main_poly = Composite::new(&vec![poly_a, poly_b, poly_c], vec![OP::MUL, OP::ADD]);
    let result = main_poly.evaluate(
      &vec![2, 3].iter().map(|x| Option::Some(Fq::from(x.clone()))).collect()
    );

    assert_eq!(result, Fq::from(735));
  }
}

