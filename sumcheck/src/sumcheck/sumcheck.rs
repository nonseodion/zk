use ark_ff::{BigInteger, PrimeField};
use polynomials::multilinear::multilinear::{MultiLinear, blow_up_left, blow_up_right};
use polynomials::multilinear::composite::{Composite, OP};
use polynomials::univariate::univariate::{evaluate, interpolate, UnivariatePolynomial};
use crate::transcipt::transcript::{HashWrapper, TranscriptTrait, Transcript};

use std::iter::repeat_n;

// Sumcheck protocol generates a proof to show that the sum of the polynomial 
// over the boolean hypercube is a particular value,
pub fn generate_partial_proof <F: PrimeField, H: HashWrapper, T: TranscriptTrait<H>> (poly: &Composite<F>, transcript: &mut T, round_polys: &mut Vec<Vec<F>>, challenges: &mut Vec<F>) -> F {
    // gets the log base 2 of the polynomial length to know how many points are needed to represent the round poly
    let mut main_poly = poly.clone();
    let degree = 2;
    let rounds = main_poly.polys[0].hypercube.len().trailing_zeros() as usize; // this is the number of variables
    let mut sum= F::zero();
    let mut sums = vec![];

    for i in 0..rounds {
        let mut reduced_poly = main_poly.reduce();
        let extra_points = reduced_poly.hypercube.len()/2;
        let mut index = 0;
        repeat_n(0, extra_points).for_each(|_| {
            let mut values = vec![Some(F::zero()); rounds-i];
            values = values.iter().enumerate().map( |x| {
                if x.0 == 0 {
                    return Some(F::from(2));
                } else {
                    // shift to right and find modulus to get the value at that point.
                    return Some(F::from(index >> (rounds-i - x.0 - 1) & 1));
                }
            }).collect();

            let result = main_poly.evaluate(&values);
            reduced_poly.hypercube.push(result);
            index += 1;
        });


        let mut round_poly = vec![];
        for j in 0..(degree + 1) {
            round_poly.push(reduced_poly.hypercube.iter().skip(j * extra_points).take(extra_points).sum());
        }

        
        sum = round_poly[0] + round_poly[1];
        sums.push(sum);
        let mut data = vec![sum];
        data.extend(&round_poly);
        let challenge = add_data_to_transcript(&data, transcript);
        challenges.push(challenge);
        main_poly = main_poly.partialEvaluate(challenge, rounds - i - 1);
        round_polys.push(round_poly);
    }


    sums[0]
}

pub fn verify_partial_proof<F: PrimeField, H: HashWrapper, T: TranscriptTrait<H>> (sum: F, polys: &Vec<Vec<F>>, transcript: &mut T) -> (F , Vec<F>) {
    let mut challenges = vec![];
    let mut challenge;
    let mut sum = sum;

    for i in 0..polys.len() {
        if sum != polys[i][0] + polys[i][1] {
            panic!("Invalid proof for partial sum check");
        }

        let mut data = vec![sum];
        data.extend(&polys[i]);
        challenge = add_data_to_transcript(&data, transcript);
        challenges.push(challenge);

        let points = polys[i].iter().enumerate().map( |x| (F::from(x.0 as u64), x.1.clone())).collect::<Vec<(F, F)>>();
        let univariate_poly = interpolate(&points);
        sum = evaluate(&univariate_poly, challenge);
    }

    (sum, challenges)
}

    // adds the proof to the Sumcheck struct
    // the proof inclueds the sums and polynomials fields in the Sumcheck struct
fn generate_proof <F: PrimeField, H: HashWrapper, T: TranscriptTrait<H>> (poly: &MultiLinear<F>, transcript: &mut T, round_polys: &mut Vec<MultiLinear<F>>) -> F {
    // initialise script with original polynomial
    add_data_to_transcript(&poly.hypercube, transcript);

    let poly_length = poly.hypercube.len();
    // this does the base 2 log of the length to get the no. of variables = no. of rounds.
    let rounds = poly_length.trailing_zeros() as usize;
    // clone the original polynomial to perform partial evaluation at random values using this polynomial
    let mut last_polynomial = poly.clone();

    // using the points from the original polynomial get the points of a univariate polynomial 
    // in the first variable for each variable
    let mut sums = vec![F::zero(); rounds];
    for j in 0..rounds {
        sums[j as usize] = last_polynomial.hypercube.iter().sum();
        // get univariate polynomial
        let mut univariate_polynomial = MultiLinear::new(vec![F::from(0), F::from(0)]);
        let poly_half_length = last_polynomial.hypercube.len()/2;

        // get the univariate polynomial for the current first variable by summing the points on the 
        // boolean hypercube where the variable is 0 and the points on the boolean hypercube where it is 1
        univariate_polynomial.hypercube[0] = last_polynomial.hypercube.iter().take(poly_half_length).sum();
        univariate_polynomial.hypercube[1] = last_polynomial.hypercube.iter().skip(poly_half_length).sum();

        let mut data = vec![sums[j]];
        data.extend(&univariate_polynomial.hypercube);
        let challenge = add_data_to_transcript(&data, transcript);
        last_polynomial = last_polynomial.partial_evaluate(challenge, rounds - j - 1);

        round_polys.push(univariate_polynomial);
    }

    sums[0]
}


fn verify_proof <F: PrimeField, H: HashWrapper, T: TranscriptTrait<H>> (poly: &MultiLinear<F>, round_polys: &Vec<MultiLinear<F>>, sum: F,  transcript: &mut T) -> bool{

    // initialise script with original polynomial
    add_data_to_transcript(&poly.hypercube, transcript);
    
    let len = round_polys.len();
    let mut challenges: Vec<Option<_>> = vec![];
    let mut claimed_sum: F = sum; // initial claimed sum is the sum of the original polynomial over the boolean hypercube
    let mut challenge;

    // check that each polynomial matches its claimed sum 
    // by comparing the sum of the polynomial with the challenge substitution in the previous univariate polynomial
    for i in 0..len {
        if claimed_sum != round_polys[i].hypercube.iter().sum() {
            return false;
        }
        let mut data = vec![claimed_sum];
        data.extend(&round_polys[i].hypercube);
        challenge = add_data_to_transcript( &data, transcript);
        challenges.push(Some(challenge));

        claimed_sum = round_polys[i].partial_evaluate(F::from(challenge),  0).hypercube[0];
    }

    let partial_evaluation = round_polys[len-1].partial_evaluate(challenges[len-1].unwrap(), 0).hypercube[0];

    let _challenges = &challenges;
    let evaluation = poly.evaluate(_challenges).hypercube[0];

    // oracle check
    if evaluation != partial_evaluation{
        return false;
    }
    return true;
}

pub fn add_data_to_transcript <F: PrimeField, H: HashWrapper, T: TranscriptTrait<H>> (data: &Vec<F>, transcript: &mut T) -> F {
    let mut bytes = vec![];
    data.iter().for_each(|x| {
        bytes.extend(x.into_bigint().to_bytes_be())
    });
    transcript.absorb(&bytes);
    let challenge = F::from_be_bytes_mod_order(&transcript.squeeze());
    return challenge;
}




#[cfg(test)]
mod test {
    // use super::
    use super::*;
    use ark_bn254::Fq;
    use sha3::{Keccak256, Digest};

    #[test]
    fn test_generate_proof(){
        let polynomial = MultiLinear::new(
            vec![
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(3),
                Fq::from(11),
                Fq::from(9),
                Fq::from(2),
                Fq::from(5),
            ]
        );

        let hasher = Keccak256::new();
        let mut transcript = Transcript::new(hasher);
        let mut round_polys = vec![];
        let sum = generate_proof(&polynomial, &mut transcript, &mut round_polys);

        transcript = Transcript::new(Keccak256::new());

        assert_eq!(verify_proof(&polynomial, &round_polys, sum, &mut transcript), true);
    }

    #[test]
    fn test_generate_partial_proof() {
        // 2a + 3
        let mut poly_a = MultiLinear::new(vec![3, 5].iter().map(|x| Fq::from(*x)).collect());
        poly_a = blow_up_right(&poly_a, 2);

        // 4b + 2a
        let mut poly_b = MultiLinear::new(vec![0, 4, 2, 6].iter().map(|x| Fq::from(*x)).collect());
        poly_b = blow_up_right(&poly_b, 1);

        // 3c + 2
        let mut poly_c = MultiLinear::new(vec![2, 5].iter().map(|x| Fq::from(*x)).collect());
        poly_c = blow_up_left(&poly_c, 2); 

        // 3c + 2
        let mut poly_d = MultiLinear::new(vec![2, 5].iter().map(|x| Fq::from(*x)).collect());
        poly_d = blow_up_left(&poly_d, 2);         

        let composite = Composite::new(
            &vec![poly_a.hypercube, poly_b.hypercube, poly_c.hypercube, poly_d.hypercube],
            vec![OP::MUL, OP::ADD, OP::MUL]
        );

        let mut round_polys = vec![];

        let hasher = Keccak256::new();
        let mut transcript = Transcript::new(hasher);
        let mut challenges = vec![];
        let initial_sum = generate_partial_proof(&composite, &mut transcript, &mut round_polys, &mut challenges);

        let hasher = Keccak256::new();
        let mut transcript = Transcript::new(hasher);
        let (sum, challenges) = verify_partial_proof(initial_sum, &round_polys, &mut transcript);

        assert_eq!(
            sum,
            composite.evaluate(&challenges.iter().map(|x| Some(*x)).collect())
        );
    }
}
