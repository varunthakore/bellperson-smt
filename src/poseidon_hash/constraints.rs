use pasta_curves::{
    group::ff::PrimeField,
};
// use ff::{Field, PrimeField};
use neptune::{
    // circuit::poseidon_hash,
    poseidon::{PoseidonConstants},
    Arity
};

use neptune::circuit2::poseidon_hash_allocated;

use bellperson::{Circuit, ConstraintSystem, SynthesisError, gadgets::num::AllocatedNum};


struct PoseidonCircuit<F: PrimeField, A: Arity<F>> {
    pub xl: F,
    pub xr: F,
    pub node: F,
    pub params: PoseidonConstants<F, A>
}

impl<F: PrimeField, A: Arity<F>> Circuit<F> for PoseidonCircuit<F, A> {
    fn synthesize<CS: ConstraintSystem<F>>(self, cs: &mut CS) -> Result<(), SynthesisError> {
        let xl = AllocatedNum::alloc(cs.namespace(|| "preimage xl"), || Ok(self.xl))?;
        let xr = AllocatedNum::alloc(cs.namespace(|| "preimage xr"), || Ok(self.xr))?;
        let node = AllocatedNum::alloc_input(cs.namespace(|| "hash node"), || Ok(self.node))?;

        let calc_node = poseidon_hash_allocated(&mut *cs, vec![xl, xr], &self.params)?;

        cs.enforce(
            || "node = calc_node", 
            |lc| lc + node.get_variable(), 
            |lc| lc + CS::one(), 
            |lc| lc + calc_node.get_variable()
        );

        Ok(())
        
    }
}

#[cfg(test)]
mod tests {
    use bellperson::Circuit;
    // use neptune::circuit2::poseidon_hash_allocated;
    use pasta_curves::pallas::Base as Fp;
    
    use neptune::poseidon::{PoseidonConstants, Poseidon};
    // use neptune::circuit::poseidon_hash;
    use generic_array::typenum::U2;
    use bellperson::{util_cs::test_cs::TestConstraintSystem ,
        // ConstraintSystem, gadgets::num::AllocatedNum
    };

    use crate::poseidon_hash::constraints::PoseidonCircuit;

    #[test]
    fn test_poseidon_circuit() {
        // Need to think of a way to test these output values
        let leaf_present  = Fp::one();
		
        let node_hash_params = PoseidonConstants::<Fp, U2>::new();
        let mut node_hasher = Poseidon::new_with_preimage(&[leaf_present.double(), leaf_present],&node_hash_params);
        let hash_value = node_hasher.hash();
        // println!("the hash value is {:?}", hash_value);

        let circ = PoseidonCircuit {
            xl: leaf_present.double(),
            xr: leaf_present,
            node: hash_value,
            params: node_hash_params
        };

        let mut cs = TestConstraintSystem::<Fp>::new();

        println!("the number of constraints are {}", cs.num_constraints());

        
        assert!(!circ.synthesize(& mut cs).is_err());

        println!("the number of constraints are {}", cs.num_constraints());

        assert!(cs.is_satisfied());

        println!("constraint satisfied is {:?}", cs.is_satisfied());

        // println!("the inputs are {:?}", cs.get_inputs());

        // let xl = AllocatedNum::alloc(cs.namespace(|| "preimage xl"), || Ok(leaf_present.double())).unwrap();
        // let xr = AllocatedNum::alloc(cs.namespace(|| "preimage xr"), || Ok(leaf_present)).unwrap();
        // let calc_node = poseidon_hash_allocated(&mut cs, vec![xl, xr], &node_hash_params).unwrap();
        // println!("calc_node hash is {:?}",calc_node.get_value().unwrap());

    }
}