use ark_ec::PrimeGroup;
use ark_ec::{
    hashing::curve_maps::parity,
    pairing::{self, Pairing, PairingOutput},
    AdditiveGroup, AffineRepr,
};
use ark_ff::{BigInteger, PrimeField};
use polynomials::multilinear::multilinear::{blow_up_left, MultiLinear};
use std::ops::Mul;
use crate::trusted_setup::generate_encrypted_lagrange_bases;

fn commit<F: PrimeField, P: Pairing>(
    poly: &MultiLinear<F>,
    encrypted_lagrange_bases_g1: &Vec<P::G1>,
) -> P::G1 {
    poly.hypercube
        .iter()
        .enumerate()
        .map(|(index, value)| encrypted_lagrange_bases_g1[index].mul_bigint(&value.into_bigint()))
        .sum()
}

fn open<F: PrimeField, P: Pairing>(
    poly: &MultiLinear<F>,
    encrypted_lagrange_bases: Vec<P::G1>,
    values: Vec<F>,
) -> (F, Vec<P::G1>) {
    let result = poly
        .evaluate(&values.iter().map(|x| Some(*x)).collect())
        .hypercube[0];
    let new_poly = MultiLinear::new(&poly.hypercube.iter().map(|x| *x - result).collect());

    let mut quotients: Vec<P::G1> = vec![];
    let variable_length = values.len();
    let mut remainder = MultiLinear::new(&new_poly.hypercube);

    for (index, value) in values.iter().enumerate() {
        let skips = variable_length - 1 - index;

        // f(1)
        let mut quotient_hypercube: Vec<F> = remainder
            .hypercube
            .iter()
            .skip(2_usize.pow(skips as u32))
            .map(|x| *x)
            .collect();

        // quotient = f(1) - f(0)
        quotient_hypercube = quotient_hypercube
            .iter()
            .enumerate()
            .map(|(index, value)| *value - remainder.hypercube[index])
            .collect();
        let quotient_poly = MultiLinear::new(&quotient_hypercube);

        let quotient = blow_up_left(&quotient_poly, (index + 1) as u32)
            .hypercube
            .iter()
            .enumerate()
            .map(|(index, value)| encrypted_lagrange_bases[index].mul_bigint(value.into_bigint()))
            .sum();

        quotients.push(quotient);
        remainder = remainder.partial_evaluate(*value, skips);
    }

    return (result, quotients);
}

// the proof is the array of quotients
fn verify_proof<F: PrimeField, P: Pairing>(
    result: F,
    commitment: P::G1,
    proofs: Vec<P::G1>,
    encrypted_taus: Vec<P::G2>,
    values: Vec<F>,
) -> bool {
    let rhs = proofs
        .iter()
        .enumerate()
        .fold(PairingOutput::ZERO, |acc, (index, proof)| {
            acc + P::pairing(
                *proof,
                encrypted_taus[index] - P::G2::generator().mul_bigint(values[index].into_bigint()),
            )
        });

    let lhs = commitment - P::G1::generator().mul_bigint(result.into_bigint());
    let expected_sum = P::pairing(lhs, P::G2::generator());
    return rhs == expected_sum;
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bls12_381::{Bls12_381, Config, G1Affine, G1Projective, G2Affine, Fr};
    use ark_ec::{
        bls12, pairing::PairingOutput, short_weierstrass::Projective, AffineRepr, CurveGroup,
    };
    use ark_ff::Field;

    #[test]
    fn test_generate_proof() {
        // 3ab + 4c
        let poly = MultiLinear::new(
            &vec![0, 4, 0, 4, 0, 4, 3, 7]
                .iter()
                .map(|x| Fr::from(*x))
                .collect(),
        );

        let taus = vec![Fr::from(5), Fr::from(2), Fr::from(3)];
        let encrypted_lagrange_bases_g1 = generate_encrypted_lagrange_bases::<Fr, Bls12_381>(&taus);

        let commitment: G1Projective = commit::<Fr, Bls12_381>(&poly, &encrypted_lagrange_bases_g1);

        // values to open proof
        let values = vec![Fr::from(6), Fr::from(4), Fr::from(0)];
        let (result, quotients) =
            open::<Fr, Bls12_381>(&poly, encrypted_lagrange_bases_g1, values.clone());

        let encrypted_taus = taus
            .iter()
            .map(|x| G2Affine::generator().mul_bigint(x.into_bigint()))
            .collect();

        assert_eq!(
            verify_proof::<Fr, Bls12_381>(result, commitment, quotients, encrypted_taus, values),
            true
        );
    }
}
