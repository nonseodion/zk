use ark_ff::PrimeField;
use ark_ec::pairing::Pairing;
use ark_ec::{PrimeGroup};

pub fn generate_encrypted_lagrange_bases<F: PrimeField, P: Pairing>(taus: &Vec<F>) -> Vec<P::G1>{
    let mut lagrange_bases = vec![F::ONE];

    for (index, tau) in taus.iter().rev().enumerate() {
        let extra1 = lagrange_bases.to_vec();
        lagrange_bases.extend(extra1);

        let extra2 = vec![*tau; 2_usize.pow(index as u32)];
        let mut tau_bases = vec![F::ONE - tau; 2_usize.pow(index as u32)];
        tau_bases.extend(extra2);

        lagrange_bases = lagrange_bases
            .iter()
            .zip(&tau_bases)
            .map(|(x1, x2)| *x1 * x2)
            .collect();
    }

    // 
    let encrypted_lagrange_bases_g1 = lagrange_bases
        .iter()
        .map(|x| P::G1::generator().mul_bigint(F::from(*x).into_bigint()))
        .collect();
    encrypted_lagrange_bases_g1
}

#[cfg(test)]
mod test {
  use super::*;
  use ark_bls12_381::{Fr, Bls12_381};

  #[test]
  fn test_generate_encrypted_bases() {
      generate_encrypted_lagrange_bases::<Fr, Bls12_381>(&vec![Fr::from(5), Fr::from(2), Fr::from(3)]);
  }
}
