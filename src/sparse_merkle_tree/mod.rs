// pub mod constraints;

// pub mod db;
// use db::HashValueDb;

// use std::hash::Hash;

use std::collections::HashMap;
// use fnv::FnvHashMap;

use num_bigint::BigUint;

use pasta_curves::group::ff::{PrimeField, PrimeFieldBits};
// use ff::PrimeField;
// use ark_ff::PrimeField;


use std::marker::PhantomData;

use pasta_curves::pallas::Base;
use neptune::poseidon::{PoseidonConstants, Poseidon};
use generic_array::typenum::{U2, U1};
use rand::thread_rng;
use bitvec::array::BitArray;


pub struct SparseMerkleTree<F: PrimeField + PrimeFieldBits, const N: usize> {
    // pub depth: usize,
    pub root: F,
    pub leaf_hash_params: PoseidonConstants<F, U1>,
    pub node_hash_params: PoseidonConstants<F, U2>
}


impl<F: PrimeField + PrimeFieldBits, const N: usize> SparseMerkleTree<F, N> {
    
    /// Create a new tree. `empty_leaf_val` is the default value for leaf of empty tree. It could be zero.
    /// Requires a database to hold leaves and nodes. The db should implement the `HashValueDb` trait
    pub fn new(
        empty_leaf_val: F,
        depth: usize,
        hash_db: &mut HashMap<F, (F,F)>,
    ) -> SparseMerkleTree<F, N> {
        assert!(depth > 0);
        let leaf_hash_params = PoseidonConstants::<F, U1>::new();
        let node_hash_params = PoseidonConstants::<F, U2>::new();
        let mut cur_hash = Poseidon::new_with_preimage(&[empty_leaf_val],&leaf_hash_params).hash();
        for _ in 0..depth {
            let val = (cur_hash.clone(), cur_hash.clone());
            cur_hash = Poseidon::new_with_preimage(&[cur_hash.clone(), cur_hash.clone()], &node_hash_params).hash();
            hash_db.insert(cur_hash.clone(), val);
        }
        Self {
            root: cur_hash,
            leaf_hash_params: leaf_hash_params,
            node_hash_params: node_hash_params
        }
    }

    pub fn update(
        &mut self,
        mut idx_in_bits: Vec<bool>,
        val: F,
        hash_db: &mut HashMap<F, (F,F)>,
    ) {

        let mut siblings = self.get_siblings_path(idx_in_bits.clone(), hash_db).unwrap().path;

        // let mut path = self.get_path(idx);
        // Reverse since path was from root to leaf but i am going leaf to root
        idx_in_bits.reverse();
        let mut cur_hash = Poseidon::new_with_preimage(&[val],&self.leaf_hash_params).hash();

        // Iterate over the bits
        for d in idx_in_bits {
            let sibling = siblings.pop().unwrap();
            let (l, r) = if d == false {
                // leaf falls on the left side
                (cur_hash, sibling)
            } else {
                // leaf falls on the right side
                (sibling, cur_hash)
            };
            let val = (l, r);
            cur_hash = Poseidon::new_with_preimage(&[l, r], &self.node_hash_params).hash();
            hash_db.put(cur_hash.clone(), val)?;
        }

        self.root = cur_hash;

        Ok(())

    }



     /// Get siblings given leaf index index
     pub fn get_siblings_path(
        &self,
        idx_in_bits: Vec<bool>,
        hash_db: &HashMap<F, (F,F)>,
    ) -> Path<F, N> {
        // let path = self.get_path(idx);
        let mut cur_node = self.root;
        let mut siblings = Vec::<F>::new();


        let mut children;
        for d in idx_in_bits {
            children = hash_db.get(&cur_node)?;
            if d == false {
                // leaf falls on the left side
                cur_node = children.0;
                siblings.push(children.1);

            } else {
                // leaf falls on the right side
                cur_node = children.1;
                siblings.push(children.0);

            }
        }
        // Ok((cur_node.clone() , siblings))
        Ok(Path{ path: siblings})
    }

}

pub fn idx_to_bits<F: PrimeField + PrimeFieldBits>(depth: usize, idx: F) -> Vec<bool> {
    let mut path: Vec<bool> = vec![];
    for d in idx.to_le_bits() {
        if path.len() >= depth {
            break;
        }
        
        if d==true {
            path.push(true)
        } else {
            path.push(false)
        }
    }

    while path.len() != depth {
        path.push(false);
    }

    path.reverse();
    path

}


pub struct Path<F, const N: usize>
where
	F: PrimeField + PrimeFieldBits,
{
	// #[allow(clippy::type_complexity)]
	path: Vec<F>, // siblings from root to leaf
	// phantom: PhantomData<H>,
}

impl<F: PrimeField + PrimeFieldBits, const N: usize> Path<F, N> {
    pub fn compute_root(
        &self,
        mut idx_in_bits: Vec<bool>,
        val: F,
    ) -> std::result::Result<F, MerkleTreeError<F>> {
        assert_eq!(self.path.len(), N);
        idx_in_bits.reverse();

        let mut cur_hash = Poseidon::new_with_preimage(&[val],&self.leaf_hash_params).hash();

        for (i, sibling) in self.path.clone().into_iter().rev().enumerate() {
            let (l, r) = if idx_in_bits[i] == false {
                // leaf falls on the left side
                (cur_hash, sibling)
            } else {
                // leaf falls on the right side
                (sibling, cur_hash)
            };
            cur_hash = hasher.hash_two(&l, &r)?;
        }
        Ok(cur_hash)

    }

    pub fn verify(
        &self,
        idx_in_bits: Vec<bool>,
        val: F,
        hasher: H,
        root: F
    ) -> std::result::Result<bool, MerkleTreeError<F>> {
        let computed_root = self.compute_root(idx_in_bits, val, hasher)?;
        Ok(computed_root == root)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{UniformRand, test_rng, Zero, One};
    use ark_bls12_381::Fr;
    
    use crate::poseidon_hash::Poseidon;
    use crate::poseidon_hash::poseidon_parameters::setup;
    use crate::sparse_merkle_tree::db::InMemoryHashValueDb;

    #[test]
    fn test_tree() {
        let rng = &mut test_rng();

        let mut db = InMemoryHashValueDb::<Fr>::new();
        // let tree_depth = 256;
        const HEIGHT: usize = 256;
        let empty_leaf_val = Fr::zero();
        
        let params = setup::<Fr>();
        let hasher = Poseidon::new(params);

        let mut tree: SparseMerkleTree<_, Poseidon<ark_ff::Fp<ark_ff::MontBackend<ark_bls12_381::FrConfig, 4>, 4>>, HEIGHT> = SparseMerkleTree::new(
            empty_leaf_val.clone(),
            hasher.clone(),
            HEIGHT,
            &mut db,
        ).unwrap();
        

        let val = Fr::one();

        let test_cases = 300;

        for _ in 0..test_cases {
            let idx = Fr::rand(rng);
            let idx_in_bits = idx_to_bits(HEIGHT, idx);
            assert!(!tree.update(idx_in_bits.clone(), val, &mut db).is_err());

            let path = tree.get_siblings_path(idx_in_bits.clone(), &db).unwrap();

            

            assert!(path.verify(idx_in_bits.clone(), val, hasher.clone(), tree.root).unwrap());
        }
        
    }
}