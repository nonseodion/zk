use ark_ff::PrimeField;
   
enum OP{
  ADD,
  MUL
}

struct Gate {
  left_input: usize,
  right_input: usize,
  op: OP,
  output: usize
}

impl Gate {
  fn new(left_input: usize, right_input: usize, op: OP, output: usize) -> Gate{
    Gate {
      left_input, right_input, op, output
    }
  }
}

struct Circuit<F: PrimeField> {
  layers: Vec<Vec<F>>,
  gates: Vec<Vec<Gate>>
}

impl <F: PrimeField> Circuit<F> {
  fn new(gates: Vec<Vec<Gate>>) -> Self{
    Circuit { layers: vec![], gates }
  }

  fn evaluate(&self, inputs: &Vec<F>) -> Vec<Vec<F>> {
    let mut layer_values = vec![];

    let mut inputs = inputs.clone();
    let mut outputs: Vec<F>;

    for layerGate in &self.gates {
      // update outputs
      outputs = vec![F::from(0); layerGate.len()];

      for gate in layerGate {
        let output = match gate.op{
          OP::ADD => inputs[gate.left_input] + inputs[gate.right_input],
          OP::MUL => inputs[gate.left_input] * inputs[gate.right_input],
        };

        outputs[gate.output] = output;
      }

      layer_values.push(inputs);
      inputs = outputs;
    }

    // push last output
    layer_values.push(inputs);
    return layer_values;
  }
}


#[cfg(test)]

mod test {
  use super::*;
  use ark_bn254::Fq;

  #[test]
  fn test_evaluate1() {
    let gates = vec![
      // layer 1
      vec![
        Gate::new(0, 1, OP::ADD, 0),
        Gate::new(2, 3, OP::MUL, 1),
      ],
      vec![
        Gate::new(0, 1, OP::ADD, 0),
      ]
    ];

    let circuit: Circuit<Fq> = Circuit::new(
      gates
    );

    let inputs = vec![ 1, 2, 3, 4 ].iter().map(|x| Fq::from(x.clone())).collect();
    let output = [Fq::from(15)];
    assert_eq!(
      circuit.evaluate(&inputs)[2], 
      output
    );
  }

  #[test]
  fn test_evaluate2() {
    let gates = vec![
      // layer 1
      vec![
        Gate::new(0, 1, OP::ADD, 0),
        Gate::new(2, 3, OP::MUL, 1),
        Gate::new(4, 5, OP::MUL, 2),
        Gate::new(6, 7, OP::MUL, 3),
      ],
      vec![
        Gate::new(0, 1, OP::ADD, 0),
        Gate::new(2, 3, OP::MUL, 1),
      ]
    ];

    let circuit: Circuit<Fq> = Circuit::new(
      gates
    );

    let inputs = vec![ 1, 2, 3, 4, 5, 6, 7, 8 ].iter().map(|x| Fq::from(x.clone())).collect();
    let output = [Fq::from(15), Fq::from(1680)];
    dbg!(circuit.evaluate(&inputs));
    assert_eq!(
      circuit.evaluate(&inputs)[2], 
      output
    );
  }
}
