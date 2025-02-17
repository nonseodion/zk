use std::cmp::max;
use std::iter::repeat_n;
use std::ops::Mul;

use polynomials::multilinear::composite::{Composite, OP as COMPOSITE_OP};
use polynomials::multilinear::multilinear::{MultiLinear, blow_up_left, blow_up_right, scalar_mul};
use polynomials::univariate::univariate::{evaluate, interpolate};
use ark_ff::PrimeField;
use crate::circut::{ Circuit, OP as CIRCUIT_OP, Gate};
use sumcheck::transcipt::transcript::{ HashWrapper, TranscriptTrait, Transcript};
use sumcheck::sumcheck::sumcheck::{add_data_to_transcript, generate_partial_proof, verify_partial_proof}; 
// use sumcheck::

#[derive(Debug)]
struct GKR_PROOF<F: PrimeField> {
  claimed_sums: Vec<F>,
  round_polys: Vec<Vec<Vec<F>>>,
  layers: Vec<Vec<F>>
}

fn generate_proof <F: PrimeField, H: HashWrapper, T: TranscriptTrait<H>> (circuit: &mut Circuit<F>, inputs: &Vec<F>, transcript: &mut T) -> GKR_PROOF<F> {
  circuit.evaluate(inputs);
  let mut gkr_proof = GKR_PROOF{
    claimed_sums: vec![],
    round_polys: vec![],
    layers: vec![]
  };

  let mut add_and_muls = vec![];
  get_add_and_muls(&circuit, &mut add_and_muls);

  let mut _w = circuit.layers[0].clone();
  if _w.len() == 1 { _w = vec![_w[0], F::zero()]; }
  let mut w_i = MultiLinear::new(_w);

  let challenges_length = next_pow_of_2(w_i.hypercube.len()) * 2;  
  let mut challenges = vec![F::zero(); challenges_length];
  add_data_to_transcript(&w_i.hypercube, transcript);
  challenges = challenges.iter().map(|_| F::from_be_bytes_mod_order(&transcript.squeeze())).collect();

  for i in 0..circuit.gates.len() {
    
    let (mut add_poly, mut mul_poly) = add_and_muls[i].clone();
    
    let w_i_plus_1 = MultiLinear::new(circuit.layers[i+1].clone());
    let blows = next_pow_of_2(w_i_plus_1.hypercube.len()) as u32;
    // blow ups
    let w_b = blow_up_right(&w_i_plus_1, blows); // blow up for c
    let w_c = blow_up_left(&w_i_plus_1, blows); // blow up for b
    let w_plus = w_b.clone() + w_c.clone();
    let w_mul = w_b.clone() * w_c.clone();

    let alpha = F::from_be_bytes_mod_order(&transcript.squeeze());
    let beta = F::from_be_bytes_mod_order(&transcript.squeeze());

    // let sum = scalar_mul(&w_i.evaluate(&vec![Some(challenges[0])]), alpha) 
    //   + scalar_mul(&w_i.evaluate(&vec![Some(challenges[1])]), beta);

    apply_alpha_beta(alpha, beta, &challenges, &mul_poly.clone(), &mut mul_poly);
    apply_alpha_beta(alpha, beta, &challenges, &add_poly.clone(), &mut add_poly);
    
    let hypercubes = vec![ add_poly, w_plus, mul_poly, w_mul].iter().map(|x| x.hypercube.clone()).collect();
    let f_poly = Composite::new(
      &hypercubes,
      vec![COMPOSITE_OP::MUL, COMPOSITE_OP::ADD, COMPOSITE_OP::MUL]
    );
    let mut round_polys = vec![];
    challenges = vec![];

    // returns challenges and initial claimed sum
    let sum = generate_partial_proof(&f_poly, transcript, &mut round_polys, &mut challenges);

    // let left = f_poly.evaluate(&challenges.iter().map(|x| Some(x.clone())).collect());
    // let points = round_polys.last().unwrap().iter().enumerate().map( |x| (F::from(x.0 as u64), x.1.clone())).collect::<Vec<(F, F)>>();
    // let univariate_poly = interpolate(&points);    
    // let right = evaluate(&univariate_poly, *challenges.last().unwrap());  
    // dbg!(left); dbg!(right);

    gkr_proof.claimed_sums.push(sum);
    gkr_proof.round_polys.push(round_polys);
    gkr_proof.layers.push(w_i.hypercube.clone());

    w_i = w_i_plus_1;
  }

  gkr_proof
}

fn verify_proof<F: PrimeField, H: HashWrapper, T: TranscriptTrait<H>> (circuit: &mut Circuit<F>, inputs: &Vec<F>, transcript: &mut T, gkr_proof: GKR_PROOF<F>) -> bool {

  let mut add_and_muls = vec![];
  get_add_and_muls(&circuit, &mut add_and_muls);
  let layers = gkr_proof.layers;
  let claimed_sums = gkr_proof.claimed_sums;
  let round_polys = gkr_proof.round_polys;

  let mut _w = circuit.layers[0].clone();
  if _w.len() == 1 { _w = vec![_w[0], F::zero()]; }
  let mut w_i = MultiLinear::new(_w);

  let challenges_length = next_pow_of_2(w_i.hypercube.len()) * 2;  
  let mut challenges = vec![F::zero(); challenges_length];
  add_data_to_transcript(&w_i.hypercube, transcript);
  challenges = challenges.iter().map(|_| F::from_be_bytes_mod_order(&transcript.squeeze())).collect();  

  for i in 0..circuit.gates.len(){
    let (mut add_poly, mut mul_poly) = add_and_muls[i].clone();
    
    let w_i_plus_1 = MultiLinear::new(circuit.layers[i+1].clone());
    let blows = next_pow_of_2(w_i_plus_1.hypercube.len()) as u32;
    // blow ups
    let w_b = blow_up_right(&w_i_plus_1, blows); // blow up for c
    let w_c = blow_up_left(&w_i_plus_1, blows); // blow up for b
    let w_plus = w_b.clone() + w_c.clone();
    let w_mul = w_b.clone() * w_c.clone();

    let alpha = F::from_be_bytes_mod_order(&transcript.squeeze());
    let beta = F::from_be_bytes_mod_order(&transcript.squeeze());

    apply_alpha_beta(alpha, beta, &challenges, &mul_poly.clone(), &mut mul_poly);
    apply_alpha_beta(alpha, beta, &challenges, &add_poly.clone(), &mut add_poly);

    let hypercubes = vec![ add_poly, w_plus, mul_poly, w_mul].iter().map(|x| x.hypercube.clone()).collect();
    let f_poly = Composite::new(
      &hypercubes,
      vec![COMPOSITE_OP::MUL, COMPOSITE_OP::ADD, COMPOSITE_OP::MUL]
    );

    let sum;
    (sum, challenges) = verify_partial_proof(claimed_sums[i], &round_polys[i], transcript);
    
    let evaluated_sum = f_poly.evaluate(&challenges.iter().map(|x| Some(x.clone())).collect());
    if sum != evaluated_sum {
      return false;
    }
  }

  // let points = round_polys[1].iter().enumerate().map( |x| (F::from(x.0 as u64), x.1.clone())).collect::<Vec<(F, F)>>();
  // let univariate_poly = interpolate(&points);
  // let result = evaluate(&univariate_poly, challenges[1]);

  return true;  
}

fn get_add_and_muls<F: PrimeField> (circuit: &Circuit<F>, add_and_muls: &mut Vec<(MultiLinear<F>, MultiLinear<F>)> ) {
  for i in 0..circuit.gates.len() {
    let gates_length = circuit.gates[i].len();
    let layer_length;
    if circuit.gates.len() <= i + 1 {
      layer_length = circuit.gates[i].iter().map(|x| max(x.left_input, x.right_input)).max().unwrap();
    } else {
      layer_length = circuit.gates[i+1].len();
    }
    let max_layer_bits = next_pow_of_2(layer_length);
    let max_gates_bits = next_pow_of_2(gates_length);

    let points_len = 1 << max_gates_bits + (max_layer_bits*2);
    let mut add_poly = vec![F::zero(); points_len];
    let mut mul_poly = vec![F::zero(); points_len];

    for (j, gate) in circuit.gates[i].iter().enumerate() {
      let index = (j << max_layer_bits * 2) // gate bits
          + (gate.left_input << max_layer_bits) // left_input bits
          + gate.right_input; // right_input bits
      match gate.op {
        CIRCUIT_OP::ADD => add_poly[index] = F::one(),
        CIRCUIT_OP::MUL => mul_poly[index] = F::one()
      }
    }

    add_and_muls.push((MultiLinear::new(add_poly), MultiLinear::new(mul_poly)));

    // f_polys.push(FPOLY::new(mul_poly, add_poly, layer.clone()))
  }  
}

fn next_pow_of_2 (no: usize) -> usize {
  let toOne = |x: usize| -> usize { if x == 0 {1} else {x}};
  toOne((no as f64).log2().ceil() as usize)
}

fn apply_alpha_beta <F: PrimeField> (alpha: F, beta: F, challenges: &Vec<F>, former_op_poly: &MultiLinear<F>, new_poly: &mut MultiLinear<F>) {
  let no_of_challenges = challenges.len()/2;
  let mut polys = vec![];

  for  skip in [0, no_of_challenges] {
    let no_of_variables = (former_op_poly.hypercube.len() as f64).log2() as usize;
    let mut _challenges: Vec<Option<F>> = challenges
        .iter()
        .skip(skip)
        .take(no_of_challenges)
        .map(|x| Some(*x))
        .collect();

    _challenges.extend(&vec![None; no_of_variables - no_of_challenges]);
    polys.push(former_op_poly.evaluate(&_challenges));
  }

  *new_poly = scalar_mul(&polys[0], alpha) + scalar_mul(&polys[1], beta);
}


#[cfg(test)]
mod test {
  use super::*;
  use ark_bn254::Fq;
  use sha3::{Keccak256, Digest};  

  #[test]
  fn test_get_add_and_muls() {
        let gates = vec![
      // layer 1
      vec![
        Gate::new(0, 1, CIRCUIT_OP::ADD, 0),
      ],
      vec![
        Gate::new(0, 1, CIRCUIT_OP::ADD, 0),
        Gate::new(2, 3, CIRCUIT_OP::MUL, 1),
      ]
    ];

    let mut circuit: Circuit<Fq> = Circuit::new(
      gates
    );

    let inputs: Vec<Fq> = vec![ 1, 2, 3, 4 ].iter().map(|x| Fq::from(x.clone())).collect();
    let mut add_and_muls = vec![];
    get_add_and_muls(&circuit, &mut add_and_muls);

    assert_eq!(
      add_and_muls[0].0.hypercube,
      vec![ 0, 1, 0, 0, 0, 0, 0, 0].iter().map(|x| Fq::from(*x as u64)).collect::<Vec<Fq>>()
    );

    assert_eq!(
      add_and_muls[0].1.hypercube,
      vec![ 0, 0, 0, 0, 0, 0, 0, 0].iter().map(|x| Fq::from(*x as u64)).collect::<Vec<Fq>>()
    );

    assert_eq!(
      add_and_muls[1].0.hypercube,
      vec![ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ].iter().map(|x| Fq::from(*x as u64)).collect::<Vec<Fq>>()
    );

    assert_eq!(
      add_and_muls[1].1.hypercube,
      vec![ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,].iter().map(|x| Fq::from(*x as u64)).collect::<Vec<Fq>>()
    );
  }

  // 4b + 2a
  #[test]
  fn test_apply_alpha_beta() {
    let poly = MultiLinear::new(vec![0, 4, 3, 7, 2, 6, 5, 9].iter().map(|x| Fq::from(*x)).collect());
    let mut new_poly: MultiLinear<Fq> = MultiLinear::new(vec![]);
    apply_alpha_beta(Fq::from(2), Fq::from(3), &vec![Fq::from(2), Fq::from(3)], &poly, &mut new_poly);

    assert_eq!(
      new_poly.hypercube,
      vec![26, 46, 41, 61].iter().map(|x| Fq::from(*x)).collect::<Vec<Fq>>()
    )
  }

  #[test]
  fn test_generate_proof() {
    let gates = vec![
      // layer 1
      vec![
        Gate::new(0, 1, CIRCUIT_OP::MUL, 0),
      ],   
      vec![
        Gate::new(0, 1, CIRCUIT_OP::ADD, 0),
        Gate::new(2, 3, CIRCUIT_OP::MUL, 1),        
      ],
      vec![
        Gate::new(0, 1, CIRCUIT_OP::ADD, 0),
        Gate::new(2, 3, CIRCUIT_OP::MUL, 1),
        Gate::new(4, 5, CIRCUIT_OP::MUL, 2),
        Gate::new(6, 7, CIRCUIT_OP::ADD, 3)      
      ]
    ];

    let mut circuit: Circuit<Fq> = Circuit::new(
      gates
    );

    let inputs: Vec<Fq> = vec![ 1, 2, 3, 4, 5, 6, 7, 8 ].iter().map(|x| Fq::from(x.clone())).collect();
    
    let mut hasher = Keccak256::new();
    let mut transcript = Transcript::new(hasher);    
    let gkr_proof = generate_proof(&mut circuit, &inputs, &mut transcript);
    
    hasher = Keccak256::new();
    transcript = Transcript::new(hasher);
    assert_eq!(
      true, 
      verify_proof(&mut circuit, &inputs, &mut transcript, gkr_proof)
    );
  }
}