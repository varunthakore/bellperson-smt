use pasta_curves::{
    group::ff::PrimeField,
};
// use ff::{Field, PrimeField};
use neptune::{
    circuit::poseidon_hash,
    poseidon::{PoseidonConstants},
    Arity
};

use bellperson::{Circuit, ConstraintSystem, SynthesisError, gadgets::num::AllocatedNum};

struct VerifyCircuit<F: PrimeField, A: Arity<F>> {
    pub path: Vec<F>,
    pub root: F,
    pub val: F,
    pub idx: Vec<bool>,
    pub hash_params: PoseidonConstants<F, A>
}

impl<F: PrimeField, A: Arity<F>> Circuit<F> for VerifyCircuit<F, A> {

    fn synthesize<CS: ConstraintSystem<F>>(self, cs: &mut CS) -> Result<(), SynthesisError> {
        
    }
    
}