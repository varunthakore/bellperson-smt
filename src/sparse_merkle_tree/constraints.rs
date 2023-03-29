use pasta_curves::{
    group::ff::PrimeField,
};
// use ff::{Field, PrimeField};
use neptune::{
    circuit::poseidon_hash,
    poseidon::{PoseidonConstants},
    Arity
};
use generic_array::typenum::{U2, U1};

use bellperson::{Circuit, ConstraintSystem, SynthesisError, gadgets::{num::AllocatedNum, boolean::{AllocatedBit, Boolean}}};

struct PathVerifyCircuit<F: PrimeField> {
    pub siblings: Vec<F>, // from root to leaf
    pub root: F,
    pub val: F,
    pub idx: Vec<bool>, // from root to leaf
    pub leaf_hash_params: PoseidonConstants<F, U1>,
    pub node_hash_params: PoseidonConstants<F, U2>
}

impl<F: PrimeField> Circuit<F> for PathVerifyCircuit<F> {

    fn synthesize<CS: ConstraintSystem<F>>(self, cs: &mut CS) -> Result<(), SynthesisError> {
        let root_var = AllocatedNum::alloc_input(cs.namespace(|| "root"), || Ok(self.root))?;
        let val_var = AllocatedNum::alloc_input(cs.namespace(|| "leaf value"), || Ok(self.val))?;
        let siblings_var = self.siblings
            .into_iter()
            .enumerate()
            .map(|(i, s)| AllocatedNum::alloc(cs.namespace(|| format!("sibling {}", i)),|| Ok(s)))
            .collect::<Result<Vec<AllocatedNum<F>>, SynthesisError>>()?
        ;

        let mut idx_var = self.idx
            .into_iter()
            .enumerate()
            .map(|(i, b)| AllocatedBit::alloc(cs.namespace(|| format!("idx {}", i)),Some(b)))
            .collect::<Result<Vec<AllocatedBit>, SynthesisError>>()?
        ;

        let mut cur_hash_var = poseidon_hash(&mut *cs, vec![val_var], &self.leaf_hash_params)?;

        idx_var.reverse(); // Going from leaf to root

        for (i, sibling) in siblings_var.clone().into_iter().rev().enumerate() {
    

            let (lc, rc) = AllocatedNum::conditionally_reverse(&mut *cs, &cur_hash_var, &sibling, &Boolean::from(idx_var[i].clone()))?;

        }

        Ok(())


    }
    
}