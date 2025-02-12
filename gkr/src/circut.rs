use ark_ff::PrimeField;

#[derive(Debug)]
pub(crate) enum OP{
  ADD,
  MUL
}

#[derive(Debug)]
pub(crate) struct Gate {
  pub(crate) left_input: usize,
  pub(crate) right_input: usize,
  pub(crate) op: OP,
  output: usize
}

impl Gate {
  pub(crate) fn new(left_input: usize, right_input: usize, op: OP, output: usize) -> Gate{
    Gate {
      left_input, right_input, op, output
    }
  }
}

#[derive(Debug)]
pub(crate) struct Circuit<F: PrimeField> {
  pub(crate) layers: Vec<Vec<F>>,
  pub(crate) gates: Vec<Vec<Gate>>
}

impl <F: PrimeField> Circuit<F> {
  pub(crate) fn new(gates: Vec<Vec<Gate>>) -> Self{
    Circuit { layers: vec![], gates }
  }

  pub(crate) fn evaluate(&mut self, inputs: &Vec<F>) -> Vec<Vec<F>> {
    let layers_len = self.gates.len() + 1;
    let mut layer_values = vec![vec![]; layers_len];

    let mut inputs = inputs.clone();
    let mut outputs: Vec<F>;
    let mut index = 0;

    for layer_gate in self.gates.iter().rev() {
      // update outputs
      outputs = vec![F::from(0); layer_gate.len()];

      for gate in layer_gate {
        let output = match gate.op{
          OP::ADD => inputs[gate.left_input] + inputs[gate.right_input],
          OP::MUL => inputs[gate.left_input] * inputs[gate.right_input],
        };
        outputs[gate.output] = output;
      }
      
      layer_values[layers_len - index - 1] = inputs;
      index += 1;
      inputs = outputs;
    }

    // push last output
    layer_values[0] = inputs;
    self.layers = layer_values.clone();
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
      ],      
      vec![
        Gate::new(0, 1, OP::ADD, 0),
        Gate::new(2, 3, OP::MUL, 1),
      ]
    ];

    let mut circuit: Circuit<Fq> = Circuit::new(
      gates
    );

    let inputs = vec![ 1, 2, 3, 4 ].iter().map(|x| Fq::from(x.clone())).collect();
    let output = [Fq::from(15)];
    assert_eq!(
      circuit.evaluate(&inputs)[0], 
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
      ],      
      vec![
        Gate::new(0, 1, OP::ADD, 0),
        Gate::new(2, 3, OP::MUL, 1),
        Gate::new(4, 5, OP::MUL, 2),
        Gate::new(6, 7, OP::MUL, 3),
      ]
    ];

    let mut circuit: Circuit<Fq> = Circuit::new(
      gates
    );

    let inputs = vec![ 1, 2, 3, 4, 5, 6, 7, 8 ].iter().map(|x| Fq::from(x.clone())).collect();
    let output = [Fq::from(15), Fq::from(1680)];
    assert_eq!(
      circuit.evaluate(&inputs)[0], 
      output
    );
  }
}
