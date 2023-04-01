pub mod gadget;
pub mod matrix;
pub mod mds;

// use pasta_curves::group::ff::PrimeField;
// use neptune::{
//     poseidon::{PoseidonConstants, Poseidon},
//     Arity
// };

// pub fn hash<F: PrimeField, A: Arity<F>>(inputs: &[F]) -> F {
//     let p = PoseidonConstants::<F, A>::new();
//     let mut hasher = Poseidon::new_with_preimage(inputs, &p);
//     hasher.hash()
// }

// #[cfg(test)]
// mod tests {
//     use pasta_curves::pallas::Base as Fp;
    
//     use neptune::poseidon::{PoseidonConstants, Poseidon};
//     use generic_array::typenum::{U1, U2};

//     #[test]
//     fn test_hash() {
//         // Need to think of a way to test these output values
//         let leaf_present  = Fp::one();

//         let leaf_hash_params = PoseidonConstants::<Fp, U1>::new();
//         let mut leaf_hasher = Poseidon::new_with_preimage(&[leaf_present],&leaf_hash_params);
//         println!("the leaf hash is {:?}",leaf_hasher.hash());
		
//         let node_hash_params = PoseidonConstants::<Fp, U2>::new();
//         let mut node_hasher = Poseidon::new_with_preimage(&[leaf_hasher.hash(), leaf_hasher.hash()],&node_hash_params);
//         println!("the node hash is {:?}",node_hasher.hash());
//     }
// }