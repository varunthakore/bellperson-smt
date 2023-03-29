use pasta_curves::{
    group::ff::{PrimeField, PrimeFieldBits},
};
// use ff::{Field, PrimeField};
use neptune::{
    // circuit::poseidon_hash,
    circuit2::poseidon_hash_allocated,
    // poseidon::{PoseidonConstants},
    sponge::circuit::SpongeCircuit,
};
// use generic_array::typenum::{U2, U1};
use bellperson::{Circuit, ConstraintSystem, SynthesisError, gadgets::{num::AllocatedNum, boolean::{AllocatedBit, Boolean}}};
use crate::sparse_merkle_tree::smt::Path;
struct PathVerifyCircuit<'a, F: PrimeField + PrimeFieldBits, const N: usize> {
    // pub siblings: Vec<F>, // from root to leaf
    pub root: F,
    pub val: F,
    pub idx: Vec<bool>, // from root to leaf
    // pub leaf_hash_params: PoseidonConstants<F, U1>,
    // pub node_hash_params: PoseidonConstants<F, U2>
    pub path: Path<'a, F, N>
}

impl<'a, F: PrimeField + PrimeFieldBits, const N: usize> Circuit<F> for PathVerifyCircuit<'a, F, N> {

    fn synthesize<CS: ConstraintSystem<F>>(self, cs: &mut CS) -> Result<(), SynthesisError> {
        let root_var = AllocatedNum::alloc_input(cs.namespace(|| "root"), || Ok(self.root))?;
        let val_var = AllocatedNum::alloc_input(cs.namespace(|| "leaf value"), || Ok(self.val))?;
        let siblings_var = self.path.siblings
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

        let mut cur_hash_var = poseidon_hash_allocated(&mut *cs, vec![val_var], &self.path.leaf_hash_params)?;

        idx_var.reverse(); // Going from leaf to root

        for (i, sibling) in siblings_var.clone().into_iter().rev().enumerate() {

            let (lc, rc) = AllocatedNum::conditionally_reverse(&mut *cs, &cur_hash_var, &sibling, &Boolean::from(idx_var[i].clone()))?;
            cur_hash_var = poseidon_hash_allocated(&mut *cs, vec![lc, rc], &self.path.node_hash_params)?;

        }

        cs.enforce(
            || "root = calc_root", 
            |lc| lc + root_var.get_variable(), 
            |lc| lc + CS::one(), 
            |lc| lc + cur_hash_var.get_variable()
        );

        Ok(())


    }
    
}

#[cfg(test)]
mod test {
    use std::collections::HashMap;

    use bellperson::Circuit;
    // use neptune::circuit2::poseidon_hash_allocated;
    // use pasta_curves::pallas::Base as Fp;
    
    // use neptune::poseidon::{PoseidonConstants, Poseidon};
    // use neptune::circuit::poseidon_hash;
    // use generic_array::typenum::U2;
    use bellperson::{util_cs::test_cs::TestConstraintSystem ,
        // ConstraintSystem, gadgets::num::AllocatedNum
    };

    use rand::thread_rng;
    use pasta_curves::{Fp, group::ff::Field};

    use crate::sparse_merkle_tree::{smt::{
        SparseMerkleTree,
        // Path,
        idx_to_bits
    }, constraints::PathVerifyCircuit};

    #[test]
    pub fn valid_path() {
        let mut db = HashMap::<String, (Fp, Fp)>::new();
        const HEIGHT: usize = 256;
        let empty_leaf_val = Fp::zero();    

        let mut tree: SparseMerkleTree<Fp, HEIGHT> = SparseMerkleTree::new(empty_leaf_val, HEIGHT, &mut db);

        let present_leaf_val = Fp::one();

        let idx = Fp::random(&mut thread_rng());
        let idx_in_bits = idx_to_bits(HEIGHT, idx);
        tree.update(idx_in_bits.clone(), present_leaf_val, &mut db);

        let path = tree.get_siblings_path(idx_in_bits.clone(), &db);

        assert!(path.verify(idx_in_bits.clone(), present_leaf_val, tree.root));

        let circ = PathVerifyCircuit {
            root: tree.root,
            val: present_leaf_val,
            idx: idx_in_bits,
            path: path
        };

        let mut cs = TestConstraintSystem::<Fp>::new();

        println!("the number of constraints are {}", cs.num_constraints());
        
        assert!(!circ.synthesize(& mut cs).is_err());

        // println!("the number of constraints are {}", cs.num_constraints());

        // assert!(cs.is_satisfied());

        // println!("constraint satisfied is {:?}", cs.is_satisfied());
    }
}