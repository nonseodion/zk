use ark_ff::{BigInteger, PrimeField};
use polynomials::multilinear::multilinear::{blow_up_left, MultiLinear};

fn commit<F: PrimeField>(poly: &MultiLinear<F>, encrypted_lagrange_bases_g1: &Vec<F>) -> F {
    poly
        .hypercube
        .iter()
        .enumerate()
        .map(|(index, value)| {
            // dbg!(encrypted_lagrange_bases[index], &value);
            encrypted_lagrange_bases_g1[index] * value
        })
        .sum()
}

fn open<F: PrimeField>(
    poly: &MultiLinear<F>, 
    encrypted_lagrange_bases: Vec<F>,
    values: Vec<F>,
) -> (F, Vec<F>) {
    let result = poly
        .evaluate(&values.iter().map(|x| Some(*x)).collect())
        .hypercube[0];
    let new_poly = MultiLinear::new(&poly.hypercube.iter().map(|x| *x - result).collect());

    let mut quotients: Vec<F> = vec![];
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
            .map(|(index, value)| {
                encrypted_lagrange_bases[index] * value
            })
            .sum();
        quotients.push(quotient);
        remainder = remainder.partial_evaluate(*value, skips);
    }

    return (result, quotients);
}

// the proof is the array of quotients
fn verify_proof<F: PrimeField>(
    result: F,
    commitment: F,
    proofs: Vec<F>,
    encrypted_taus: Vec<F>,
    values: Vec<F>,
) -> bool {
    let rhs = proofs
        .iter()
        .enumerate()
        .fold( F::zero(), |acc, (index, proof)| {
            acc + *proof * (encrypted_taus[index] - values[index])
        });

    let lhs = commitment - result;
    return rhs == lhs;
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::{Fq};

    #[test]
    fn test_generate_proof() {
        // 3ab + 4c
        let poly = MultiLinear::new(
            &vec![0, 4, 0, 4, 0, 4, 3, 7]
            .iter()
            .map(|x| Fq::from(*x))
            .collect(),
        );

        let lagrange_bases = vec![-8, 12, 16, -24, 10, -15, -20, 30]
            .iter()
            .map(|x| {
                Fq::from(*x)
            })
            .collect();

        let commitment = commit::<Fq>(&poly, &lagrange_bases);

        // values to open proof
        let values = vec![Fq::from(6), Fq::from(4), Fq::from(0)];
        let (result, quotients) = open::<Fq>(
            &poly,
            lagrange_bases,
            values.clone(),
        );

        let encrypted_taus = vec![Fq::from(5), Fq::from(2), Fq::from(3)];

        assert_eq!(
            verify_proof::<Fq>(result, commitment, quotients, encrypted_taus, values),
            true
        );

    }
}
