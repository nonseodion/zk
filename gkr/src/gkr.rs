use std::cmp::max;
use std::hash::Hash;

use polynomials::multilinear::composite::{Composite, OP as COMPOSITE_OP};
use polynomials::multilinear::multilinear::{MultiLinear, blow_up_left, blow_up_right, scalar_mul};
use ark_ff::{ biginteger::BigInteger, PrimeField};
use crate::circut::{ Circuit, OP as CIRCUIT_OP, Gate};
use sumcheck::transcipt::transcript::{ HashWrapper, TranscriptTrait, Transcript};

use sumcheck::sumcheck::sumcheck::{add_data_to_transcript, generate_partial_proof, verify_partial_proof};
use ark_ec::{
    pairing::Pairing,
    AdditiveGroup
};
use multilinear_kzg::multilinear_kzg::multilinear_kzg::{commit, open, verify_proof as verify_kzg_proof};

#[derive(Debug)]
struct KZG_PROOF<P: Pairing> {
  commitment: P::G1,
  b_quotients: Vec<P::G1>,
  c_quotients: Vec<P::G1>,
}

#[derive(Debug)]
struct GKR_PROOF<P: Pairing, F: PrimeField> {
  kzg_proof: KZG_PROOF<P>,
  claimed_sums: Vec<F>,
  round_polys: Vec<Vec<Vec<F>>>,
  evaluations: Vec<(F, F)>,
  output: Vec<F>
}

fn generate_proof <F: PrimeField, P:Pairing, H: HashWrapper, T: TranscriptTrait<H>> (
  circuit: &mut Circuit<F>, 
  inputs: &Vec<F>, 
  transcript: &mut T,
  encrypted_lagrange_bases: &Vec<P::G1>
) -> GKR_PROOF<P, F> {
  let layer_evaluations = circuit.evaluate(inputs);
  let mut kzg_proof = KZG_PROOF{ 
    commitment: P::G1::ZERO,
    b_quotients: vec![],
    c_quotients: vec![]
  };
  let commitment = commit::<F, P>(&MultiLinear::new(inputs), encrypted_lagrange_bases);
  transcript.absorb(&commitment.to_string().as_bytes().to_vec());
  kzg_proof.commitment = commitment;

  let mut gkr_proof = GKR_PROOF{
    kzg_proof,
    claimed_sums: vec![],
    round_polys: vec![],
    evaluations: vec![],
    output: vec![]
  };

  let mut add_and_muls = vec![];
  get_add_and_muls(&circuit, &mut add_and_muls);

  let mut _w = layer_evaluations[0].clone();

  if _w.len() == 1 { _w = vec![_w[0], F::zero()]; }
  let w_i = MultiLinear::new(&_w);

  let challenges_length = next_pow_of_2(w_i.hypercube.len());
  let mut challenges = vec![F::zero(); challenges_length];
  add_data_to_transcript(&w_i.hypercube, transcript);
  challenges = challenges.iter().map(|_| F::from_be_bytes_mod_order(&transcript.squeeze())).collect();

  for i in 0..circuit.gates.len() {
    let (mut add_poly, mut mul_poly) = add_and_muls[i].clone();
    
    let w_i_plus_1 = MultiLinear::new(&layer_evaluations[i+1]);
    let blows = next_pow_of_2(w_i_plus_1.hypercube.len()) as u32;
    // blow ups
    let w_b = blow_up_right(&w_i_plus_1, blows); // blow up for c
    let w_c = blow_up_left(&w_i_plus_1, blows); // blow up for b
    let w_plus = MultiLinear::new(&w_b.hypercube) + MultiLinear::new(&w_c.hypercube);
    let w_mul = MultiLinear::new(&w_b.hypercube) * MultiLinear::new(&w_c.hypercube);

    if i != 0 {    
      let alpha = F::from_be_bytes_mod_order(&transcript.squeeze());
      let beta = F::from_be_bytes_mod_order(&transcript.squeeze());
      add_poly = apply_alpha_beta(alpha, beta, &challenges, &add_poly);
      mul_poly = apply_alpha_beta(alpha, beta, &challenges, &mul_poly);
    } else {
      add_poly = add_poly.evaluate(&challenges.iter().map(|x| Some(*x)).collect());
      mul_poly = mul_poly.evaluate(&challenges.iter().map(|x| Some(*x)).collect());
    }

    let hypercubes = vec![ 
      add_poly, 
      MultiLinear::new(&w_plus.hypercube), 
      mul_poly, 
      MultiLinear::new(&w_mul.hypercube)
    ].iter().map(|x| x.hypercube.clone()).collect();

    let f_poly = Composite::new(
      &hypercubes,
      vec![COMPOSITE_OP::MUL, COMPOSITE_OP::ADD, COMPOSITE_OP::MUL]
    );
    let mut round_polys = vec![];
    challenges = vec![];
    // returns challenges and initial claimed sum
    let sum = generate_partial_proof(&f_poly, transcript, &mut round_polys, &mut challenges);
    let (w_b_eval, w_c_eval);
    let b_challenges = challenges.iter().take(blows as usize).map(|x| *x).collect();
    let c_challenges = challenges.iter().skip(blows as usize).map(|x| *x).collect();

    if i == circuit.gates.len()-1 {
      // last_layer compute kzg proofs
      (w_b_eval, gkr_proof.kzg_proof.b_quotients) = open::<F, P>(&w_i_plus_1, encrypted_lagrange_bases, &b_challenges);
      (w_c_eval, gkr_proof.kzg_proof.c_quotients) = open::<F, P>(&w_i_plus_1, encrypted_lagrange_bases, &c_challenges);
    } else {
      w_b_eval = w_i_plus_1.evaluate(&b_challenges.iter().map(|x| Some(*x)).collect()).hypercube[0];
      w_c_eval = w_i_plus_1.evaluate(&c_challenges.iter().map(|x| Some(*x)).collect()).hypercube[0];      
    }
    
    add_data_to_transcript(&vec![w_b_eval, w_c_eval], transcript);

    gkr_proof.claimed_sums.push(sum);
    gkr_proof.round_polys.push(round_polys);
    gkr_proof.evaluations.push((w_b_eval, w_c_eval));
  }

  gkr_proof.output = layer_evaluations[0].clone();

  gkr_proof
}

fn verify_proof<F: PrimeField, P: Pairing, H: HashWrapper, T: TranscriptTrait<H>> (
  circuit: &mut Circuit<F>, 
  encrypted_taus: &Vec<P::G2>, 
  transcript: &mut T, 
  gkr_proof: GKR_PROOF<P, F>
) -> bool {
  transcript.absorb(&gkr_proof.kzg_proof.commitment.to_string().as_bytes().to_vec());
  let mut add_and_muls = vec![];
  get_add_and_muls(&circuit, &mut add_and_muls);

  let evaluations = gkr_proof.evaluations;
  let claimed_sums = gkr_proof.claimed_sums;
  let round_polys = gkr_proof.round_polys;

  let mut _w = gkr_proof.output;
  if _w.len() == 1 { _w.push(F::zero()) }
  let w_i = MultiLinear::new(&_w);

  let challenges_length = next_pow_of_2(w_i.hypercube.len());  
  let mut challenges = vec![F::zero(); challenges_length];

  add_data_to_transcript(&w_i.hypercube, transcript);
  challenges = challenges.iter().map(|_| F::from_be_bytes_mod_order(&transcript.squeeze())).collect();  

  let last_index = circuit.gates.len()-1;
  for i in 0..circuit.gates.len(){
    // follows order of transcript call to ensure it gets the same challenges as prover
    // so alpha and beta are fetched before verify_partial_proof is called even though they aren't used
    let (mut alpha, mut beta)  = (F::zero(), F::zero());
    if i != 0 {
      alpha = F::from_be_bytes_mod_order(&transcript.squeeze());
      beta = F::from_be_bytes_mod_order(&transcript.squeeze()); 
    }

    let (sum, new_challenges, success) = verify_partial_proof(claimed_sums[i], &round_polys[i], transcript);
    if !success {return false};
    let (mut add_poly, mut mul_poly) = add_and_muls[i].clone();

    let (w_b_eval, w_c_eval, w_plus, w_mul);
    if i < last_index {
      (w_b_eval, w_c_eval) = evaluations[i];
      (w_plus , w_mul) = (w_b_eval + w_c_eval, w_b_eval * w_c_eval);
    } else {
      // last layer 
      let challenges_len = new_challenges.len() / 2;
      let b_challenges = new_challenges.iter().take(challenges_len).map(|x| *x).collect();
      let c_challenges = new_challenges.iter().skip(challenges_len).take(challenges_len).map(|x| *x).collect();
      (w_b_eval, w_c_eval) = evaluations[i];
      // verify results using kzg proofs
      let commitment = gkr_proof.kzg_proof.commitment;
      let (check1, check2) = (
        verify_kzg_proof::<F,P>(w_b_eval, commitment, &gkr_proof.kzg_proof.b_quotients, encrypted_taus, &b_challenges),
        verify_kzg_proof::<F,P>(w_c_eval, commitment, &gkr_proof.kzg_proof.c_quotients, encrypted_taus, &c_challenges)
      );
      if !check1 || !check2 {return false;}
      (w_plus, w_mul) = (w_b_eval + w_c_eval, w_b_eval * w_c_eval);
    }

    add_data_to_transcript(&vec![w_b_eval, w_c_eval], transcript);
    
    
    if i != 0 {
      mul_poly = apply_alpha_beta(alpha, beta, &challenges, &mul_poly);
      add_poly = apply_alpha_beta(alpha, beta, &challenges, &add_poly);
    } else {
      mul_poly = mul_poly.evaluate(&challenges.iter().map(|x| Some(*x)).collect());            
      add_poly = add_poly.evaluate(&challenges.iter().map(|x| Some(*x)).collect());
    }

      mul_poly = scalar_mul(&mul_poly, w_mul);
      add_poly = scalar_mul(&add_poly, w_plus);

    let f_poly = mul_poly + add_poly;
    let evaluated_sum = f_poly.evaluate(&new_challenges.iter().map(|x| Some(*x)).collect()).hypercube[0];
    if sum != evaluated_sum {
      return false;
    }

    challenges = new_challenges;
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

    add_and_muls.push((MultiLinear::new(&add_poly), MultiLinear::new(&mul_poly)));

    // f_polys.push(FPOLY::new(mul_poly, add_poly, layer.clone()))
  }  
}

fn next_pow_of_2 (no: usize) -> usize {
  let toOne = |x: usize| -> usize { if x == 0 {1} else {x}};
  toOne((no as f64).log2().ceil() as usize)
}

fn apply_alpha_beta <F: PrimeField> (alpha: F, beta: F, challenges: &Vec<F>, former_op_poly: &MultiLinear<F>) -> MultiLinear<F> {
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

  scalar_mul(&polys[0], alpha) + scalar_mul(&polys[1], beta)
}


#[cfg(test)]
mod test {
  use super::*;
  use multilinear_kzg::trusted_setup::trusted_setup::generate_encrypted_lagrange_bases;
  use sha3::{Keccak256, Digest}; 
  use ark_bls12_381::{Bls12_381, G2Affine, Fr};
  use ark_ec::AffineRepr;

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

    let mut circuit: Circuit<Fr> = Circuit::new(
      gates
    );

    let inputs: Vec<Fr> = vec![ 1, 2, 3, 4 ].iter().map(|x| Fr::from(*x)).collect();
    let mut add_and_muls = vec![];
    get_add_and_muls(&circuit, &mut add_and_muls);

    assert_eq!(
      add_and_muls[0].0.hypercube,
      vec![ 0, 1, 0, 0, 0, 0, 0, 0].iter().map(|x| Fr::from(*x as u64)).collect::<Vec<Fr>>()
    );

    assert_eq!(
      add_and_muls[0].1.hypercube,
      vec![ 0, 0, 0, 0, 0, 0, 0, 0].iter().map(|x| Fr::from(*x as u64)).collect::<Vec<Fr>>()
    );

    assert_eq!(
      add_and_muls[1].0.hypercube,
      vec![ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ].iter().map(|x| Fr::from(*x as u64)).collect::<Vec<Fr>>()
    );

    assert_eq!(
      add_and_muls[1].1.hypercube,
      vec![ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,].iter().map(|x| Fr::from(*x as u64)).collect::<Vec<Fr>>()
    );
  }

  // 4b + 2a
  #[test]
  fn test_apply_alpha_beta() {
    let poly = MultiLinear::new(&vec![0, 4, 3, 7, 2, 6, 5, 9].iter().map(|x| Fr::from(*x)).collect());
    let new_poly: MultiLinear<Fr> = 
      apply_alpha_beta(Fr::from(2), Fr::from(3), &vec![Fr::from(2), Fr::from(3)], &poly);

    assert_eq!(
      new_poly.hypercube,
      vec![26, 46, 41, 61].iter().map(|x| Fr::from(*x)).collect::<Vec<Fr>>()
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

    let mut circuit: Circuit<Fr> = Circuit::new(
      gates
    );

    let inputs: Vec<Fr> = vec![ 1, 2, 3, 4, 5, 6, 7, 8 ].iter().map(|x| Fr::from(*x)).collect();

    let taus = vec![Fr::from(5), Fr::from(2), Fr::from(3)];
    let encrypted_lagrange_bases = generate_encrypted_lagrange_bases::<Fr, Bls12_381>(&taus);    
    let encrypted_taus = taus
        .iter()
        .map(|x| G2Affine::generator().mul_bigint(x.into_bigint()))
        .collect();    
    
    let mut hasher = Keccak256::new();
    let mut transcript = Transcript::new(hasher);    
    let gkr_proof: GKR_PROOF<Bls12_381, Fr> = generate_proof(&mut circuit, &inputs, &mut transcript, &encrypted_lagrange_bases);
    
    hasher = Keccak256::new();
    transcript = Transcript::new(hasher);

    assert_eq!(
      true, 
      verify_proof(&mut circuit, &encrypted_taus, &mut transcript, gkr_proof)
    );
  }
}