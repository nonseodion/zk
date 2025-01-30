use ark_ff::PrimeField;
use polynomials::univariate::univariate::{UnivariatePolynomial, interpolate, evaluate};
use rand::{thread_rng, Rng};

// This file implements the Shamir Secret Sharing protocol.
// The protocol allows a secret to be shared into different parts using a polynomial.
// The secret is placed at the zero index of the polynomial.
// random x and y values are generated and used as the points of the polynomial.
// These points are used to get the original polynomial. 
// Then the orginal polynomial is used to calculate for random values of x to get
// their corresponding y values.
// These new x and y values are the passwords.

// the secret to share
// the minimum number of points neded to get the secret.
fn generate_passwords<F: PrimeField>(secret: F, quorum: u64) -> Vec<(F, F)> {
    let mut points: Vec<(F, F)> = vec![]; // holds the passwords/genrated points
    let mut rng = thread_rng(); // random number generator

    points.push((F::from(0), secret)); 
    // generate quorum random numbers and store in r
    for i in 1..quorum {
        let random_number: u64 = rng.gen();
        points.push((F::from(i), F::from(random_number)));
    }

    let polynomial = interpolate(&points);

    for i in 0..quorum {
        let x = F::from(rng.gen::<u64>());
        let y = evaluate(&polynomial, x);
        points[i as usize] = (x, y);
    }

    return points;
}


#[cfg(test)]
mod test {
    use ark_bn254::Fq;
    use super::*;

    use crate::generate_passwords;

    #[test]
    fn test_generate_passwords1() {
        let secret = Fq::from(17487101313200509019_u64);
        let passwords = generate_passwords(secret, 2);
        let original_polynomial = interpolate(&passwords);
        assert_eq!(evaluate(&original_polynomial, Fq::from(0)), secret);
    }
}