use std::ops::Mul;

use polynomials::multilinear::multilinear::MultiLinear;
use ark_ff::PrimeField;
use crate::circut::{self, Circuit, OP, Gate};

#[derive(Debug)]
struct FPOLY<F: PrimeField> {
  mul_poly: MultiLinear<F>,
  add_poly: MultiLinear<F>,
  layer: MultiLinear<F>
}

impl <F: PrimeField> FPOLY<F> {
  fn new(mul_poly: Vec<F>, add_poly: Vec<F>, layer: Vec<F>) -> Self {
    
    return FPOLY {
      mul_poly: MultiLinear::new(mul_poly),
      add_poly: MultiLinear::new(add_poly),
      layer: MultiLinear::new(layer)
    }
  }
}

#[derive(Debug)]
struct GKR<F: PrimeField> {
  f_polys: Vec<FPOLY<F>>,
  circuit: Circuit<F>
}

impl <F: PrimeField> GKR<F> {
  fn new(circuit: Circuit<F>) -> Self {
    return GKR {
      circuit,
      f_polys: vec![]
    }
  }

  fn generate_proof(&mut self, inputs: &Vec<F>) {

    let toOne = |x: usize| -> usize { if x == 0 {1} else {x}};
    self.circuit.evaluate(inputs);

    for i in 0..self.circuit.gates.len() {
      let gates = &self.circuit.gates[i];
      let layer = &self.circuit.layers[i+1];
      let max_layer_bits = toOne((layer.len() as f64).log2().ceil() as usize);
      let max_gates_bits = toOne((gates.len() as f64).log2().ceil() as usize);

      let points_len = 1 << max_gates_bits + (max_layer_bits*2);
      let mut add_poly = vec![F::zero(); points_len];
      let mut mul_poly = vec![F::zero(); points_len];

      for (j, gate) in self.circuit.gates[i].iter().enumerate() {
        let index = (j << max_layer_bits * 2) // gate bits
            + (gate.left_input << max_layer_bits) // left_input bits
            + gate.right_input; // right_input bits
        match gate.op {
          OP::ADD => add_poly[index] = F::one(),
          OP::MUL => mul_poly[index] = F::one()
        }
      }

      self.f_polys.push(FPOLY::new(mul_poly, add_poly, layer.clone()))
    }
  }
}


#[cfg(test)]
mod test {
  use super::*;
  use ark_bn254::Fq;

  #[test]
  fn test_generate_proof() {
    let gates = vec![
      // layer 1
      vec![
        Gate::new(0, 1, OP::ADD, 0),
      ],      
      vec![
        Gate::new(0, 1, OP::ADD, 0),
        Gate::new(2, 3, OP::MUL, 1),
      ]
    ];

    let mut circuit: Circuit<Fq> = Circuit::new(
      gates
    );

    let inputs: Vec<Fq> = vec![ 1, 2, 3, 4 ].iter().map(|x| Fq::from(x.clone())).collect();
    let output = [Fq::from(15)];

    let mut gkr = GKR::new(circuit);

    gkr.generate_proof(&inputs);

  }
}