use std::collections::HashMap;
use pasta_curves::group::ff::{PrimeField, PrimeFieldBits};
use neptune::poseidon::{PoseidonConstants, Poseidon};
use generic_array::typenum::{U2, U1};


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
        hash_db: &mut HashMap<String, (F,F)>,
    ) -> SparseMerkleTree<F, N> {
        assert!(depth > 0);
        let leaf_hash_params = PoseidonConstants::<F, U1>::new();
        let node_hash_params = PoseidonConstants::<F, U2>::new();
        let mut cur_hash = Poseidon::new_with_preimage(&[empty_leaf_val],&leaf_hash_params).hash();
        for _ in 0..depth {
            let val = (cur_hash.clone(), cur_hash.clone());
            cur_hash = Poseidon::new_with_preimage(&[cur_hash.clone(), cur_hash.clone()], &node_hash_params).hash();
            hash_db.insert(format!("{:?}",cur_hash.clone()), val);
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
        hash_db: &mut HashMap<String, (F,F)>,
    ) {

        let mut siblings = self.get_siblings_path(idx_in_bits.clone(), hash_db).path;

        // let mut path = self.get_path(idx);

        // Reverse since path was from root to leaf but I am going leaf to root
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
            hash_db.insert(format!("{:?}",cur_hash.clone()), val);
        }

        self.root = cur_hash;

    }



     /// Get siblings given leaf index index
     pub fn get_siblings_path(
        &self,
        idx_in_bits: Vec<bool>,
        hash_db: &HashMap<String, (F,F)>,
    ) -> Path<F, N> {
        // let path = self.get_path(idx);
        let mut cur_node = self.root;
        let mut siblings = Vec::<F>::new();


        let mut children;
        for d in idx_in_bits {
            children = hash_db.get(&format!("{:?}",cur_node.clone())).unwrap();
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
        Path{ 
            path: siblings,
            leaf_hash_params: &self.leaf_hash_params,
            node_hash_params: &self.node_hash_params
        }
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


pub struct Path<'a, F, const N: usize>
where
	F: PrimeField + PrimeFieldBits,
{
	// #[allow(clippy::type_complexity)]
	pub path: Vec<F>, // siblings from root to leaf
    pub leaf_hash_params: &'a PoseidonConstants<F, U1>,
    pub node_hash_params: &'a PoseidonConstants<F, U2>
	// phantom: PhantomData<H>,
}

impl<'a, F: PrimeField + PrimeFieldBits, const N: usize> Path<'a, F, N> {
    pub fn compute_root(
        &self,
        mut idx_in_bits: Vec<bool>,
        val: F,
    ) -> F {
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
            cur_hash = Poseidon::new_with_preimage(&[l,r],&self.node_hash_params).hash();
        }
        cur_hash

    }

    pub fn verify(
        &self,
        idx_in_bits: Vec<bool>,
        val: F,
        root: F
    ) -> bool {
        let computed_root = self.compute_root(idx_in_bits, val);
        computed_root == root
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    use rand::thread_rng;
    use pasta_curves::{Fp, group::ff::Field};

    #[test]
    fn test_smt() {
        let mut db = HashMap::<String, (Fp, Fp)>::new();
        const HEIGHT: usize = 256;
        let empty_leaf_val = Fp::zero();    

        let mut tree: SparseMerkleTree<Fp, HEIGHT> = SparseMerkleTree::new(empty_leaf_val, HEIGHT, &mut db);

        let val = Fp::one();

        let test_cases = 300;

        for _ in 0..test_cases {
            let idx = Fp::random(&mut thread_rng());
            let idx_in_bits = idx_to_bits(HEIGHT, idx);
            tree.update(idx_in_bits.clone(), val, &mut db);

            let path = tree.get_siblings_path(idx_in_bits.clone(), &db);

            

            assert!(path.verify(idx_in_bits.clone(), val, tree.root));
        }
        
    }
}