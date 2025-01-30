use std::marker::PhantomData;
use std::fmt::Debug;
use ark_ff::{BigInteger, PrimeField};
use polynomials::multilinear::multilinear::MultiLinear;
mod transcript;

use transcript::{HashWrapper, TranscriptTrait};

// Sumcheck protocol generates a proof to show that the sum of the polynomial 
// 6over the boolean hypercube is a particular value
#[derive(Debug)]
pub struct Sumcheck<F: PrimeField, H: HashWrapper, T: TranscriptTrait<H>> {
    polynomials: Vec<MultiLinear<F>>, // univariate_polynomials for proof
    original_polynomial: MultiLinear<F>, 
    sum: F, // claimed sums corresponding to polynomials
    transcript: T, // transcript shared by proof generation and verifier
    phantom: PhantomData<H>
}


impl<F: PrimeField, H: HashWrapper, T: TranscriptTrait<H>> Sumcheck<F, H, T> {
    // creates a Sumcheck Protocol instance
    fn new(polynomial: MultiLinear<F>, transcript: T) -> Self {
        Sumcheck{
            sum: polynomial.hypercube.iter().sum(),
            original_polynomial: polynomial,
            polynomials: vec![],
            transcript,
            phantom: PhantomData
        }
    }

    // adds the proof to the Sumcheck struct
    // the proof inclueds the sums and polynomials fields in the Sumcheck struct
    fn generate_proof(&mut self) {
        // initialise script with original polynomial
        self.add_original_polynomial_to_transcript();

        let poly_length = self.original_polynomial.hypercube.len();
        // this does the base 2 log of the length to get the no. of variables = no. of rounds.
        let rounds = poly_length.trailing_zeros();
        // clone the original polynomial to perform partial evaluation at random values using this polynomial
        let mut last_polynomial = self.original_polynomial.clone();

        for j in 0..rounds {

            let sum = last_polynomial.hypercube.iter().sum();
            // get univariate polynomial
            let mut univariate_polynomial = MultiLinear::new(vec![F::from(0), F::from(0)]);
            let poly_half_length = last_polynomial.hypercube.len()/2;

            for i in 0..poly_half_length {
                let zero_hypercube = last_polynomial.hypercube[i];
                let one_hypercube = last_polynomial.hypercube[i+poly_half_length];
                univariate_polynomial.hypercube[0] += zero_hypercube;
                univariate_polynomial.hypercube[1] += one_hypercube; // same as constant
            }
            self.polynomials.push(univariate_polynomial);

            let challenge = self.add_sum_and_univariate_polynomial_to_transcript(sum, j as usize);
            last_polynomial = last_polynomial.partial_evaluate(challenge, rounds - j - 1);
            // add univariate polynomial to polynomial list
        }

        // clear transcript to allow future calls start from an empty state.
        self.transcript.clear();
    }

    fn verify_proof(&mut self) -> bool{

        if self.original_polynomial.hypercube.iter().sum::<F>() != self.sum {
            return false;
        }

        // initialise script with original polynomial
        self.add_original_polynomial_to_transcript();
        
        let len = self.polynomials.len();
        let mut challenges = vec![];
        let mut claimed_sum: F = self.sum;
        let mut challenge;

        // 
        for i in 0..len {
            if claimed_sum != self.polynomials[i].hypercube.iter().sum() {
                return false;
            }

            challenge = self.add_sum_and_univariate_polynomial_to_transcript(claimed_sum, i);
            challenges.push(Some(challenge));            

            claimed_sum = self.polynomials[i].partial_evaluate(F::from(challenge),  0).hypercube[0];
        }

        let partial_evaluation = self.polynomials[len-1].partial_evaluate(challenges[len-1].unwrap(), 0).hypercube[0];
        let evaluation = self.original_polynomial.evaluate(challenges).hypercube[0];

        // oracle check
        if evaluation != partial_evaluation{
            return false;
        }
        
        // clear transcript to remove appended data
        self.transcript.clear();
        return true;
    }

    fn add_original_polynomial_to_transcript (&mut self) {
        let mut original_polynomial_bytes = vec![];
        self.original_polynomial.hypercube.iter().for_each(|x| {
            original_polynomial_bytes.extend(x.into_bigint().to_bytes_be())
        });
        self.transcript.append(&original_polynomial_bytes);
    }

    fn add_sum_and_univariate_polynomial_to_transcript (&mut self, sum: F, index: usize) -> F {
            let polynomial = &self.polynomials[index];
            // add sum and polynomial to previous transcript data and hash
            let mut data = sum.into_bigint().to_bytes_be();
            data.extend(polynomial.hypercube[0].into_bigint().to_bytes_be());
            data.extend(polynomial.hypercube[1].into_bigint().to_bytes_be());

            let challenge = F::from_be_bytes_mod_order(&mut self.transcript.append( &data ) );

            return challenge;
    }
}


#[cfg(test)]
mod test {
    use crate::transcript::Transcript;

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
        let transcript = Transcript::new(hasher);

        let mut protocol = Sumcheck::new(polynomial, transcript);
        protocol.generate_proof();

        assert_eq!(protocol.verify_proof(), true);
    }
}
