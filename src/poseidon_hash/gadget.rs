// Copied from https://github.com/lurk-lab/neptune/blob/master/src/circuit2.rs

use pasta_curves::{
    group::ff::PrimeField,
};
use neptune::{
    poseidon::{PoseidonConstants, Arity},
};

use crate::poseidon_hash::matrix::Matrix;
use crate::poseidon_hash::mds::SparseMatrix;
use bellperson::gadgets::boolean::Boolean;
use bellperson::gadgets::num::{self, AllocatedNum};
use bellperson::{ConstraintSystem, LinearCombination, SynthesisError};
use std::marker::PhantomData;

/// Similar to `num::Num`, we use `Elt` to accumulate both values and linear combinations, then eventually
/// extract into a `num::AllocatedNum`, enforcing that the linear combination corresponds to the result.
#[derive(Clone)]
pub enum Elt<Scalar: PrimeField> {
    Allocated(AllocatedNum<Scalar>),
    Num(num::Num<Scalar>),
}

impl<Scalar: PrimeField> From<AllocatedNum<Scalar>> for Elt<Scalar> {
    fn from(allocated: AllocatedNum<Scalar>) -> Self {
        Self::Allocated(allocated)
    }
}

impl<Scalar: PrimeField> Elt<Scalar> {
    pub fn is_allocated(&self) -> bool {
        matches!(self, Self::Allocated(_))
    }

    pub fn is_num(&self) -> bool {
        matches!(self, Self::Num(_))
    }

    pub fn num_from_fr<CS: ConstraintSystem<Scalar>>(fr: Scalar) -> Self {
        let num = num::Num::<Scalar>::zero();
        Self::Num(num.add_bool_with_coeff(CS::one(), &Boolean::Constant(true), fr))
    }

    pub fn ensure_allocated<CS: ConstraintSystem<Scalar>>(
        &self,
        cs: &mut CS,
        enforce: bool,
        hnum: i32
    ) -> Result<AllocatedNum<Scalar>, SynthesisError> {
        match self {
            Self::Allocated(v) => Ok(v.clone()),
            Self::Num(num) => {
                let v = AllocatedNum::alloc(cs.namespace(|| format!("hash num {} : allocate for Elt::Num", hnum)), || {
                    num.get_value().ok_or(SynthesisError::AssignmentMissing)
                })?;

                if enforce {
                    cs.enforce(
                        || format!("hash num {} : enforce num allocation preserves lc", hnum),
                        |_| num.lc(Scalar::one()),
                        |lc| lc + CS::one(),
                        |lc| lc + v.get_variable(),
                    );
                }
                Ok(v)
            }
        }
    }

    pub fn val(&self) -> Option<Scalar> {
        match self {
            Self::Allocated(v) => v.get_value(),
            Self::Num(num) => num.get_value(),
        }
    }

    pub fn lc(&self) -> LinearCombination<Scalar> {
        match self {
            Self::Num(num) => num.lc(Scalar::one()),
            Self::Allocated(v) => LinearCombination::<Scalar>::zero() + v.get_variable(),
        }
    }

    /// Add two Elts and return Elt::Num tracking the calculation.
    #[allow(clippy::should_implement_trait)]
    pub fn add(self, other: Elt<Scalar>) -> Result<Elt<Scalar>, SynthesisError> {
        match (self, other) {
            (Elt::Num(a), Elt::Num(b)) => Ok(Elt::Num(a.add(&b))),
            (a, b) => Ok(Elt::Num(a.num().add(&b.num()))),
        }
    }

    pub fn add_ref(self, other: &Elt<Scalar>) -> Result<Elt<Scalar>, SynthesisError> {
        match (self, other) {
            (Elt::Num(a), Elt::Num(b)) => Ok(Elt::Num(a.add(b))),
            (a, b) => Ok(Elt::Num(a.num().add(&b.num()))),
        }
    }

    /// Scale
    pub fn scale<CS: ConstraintSystem<Scalar>>(
        self,
        scalar: Scalar,
    ) -> Result<Elt<Scalar>, SynthesisError> {
        match self {
            Elt::Num(num) => Ok(Elt::Num(num.scale(scalar))),
            Elt::Allocated(a) => Elt::Num(a.into()).scale::<CS>(scalar),
        }
    }

    /// Square
    pub fn square<CS: ConstraintSystem<Scalar>>(
        &self,
        mut cs: CS,
    ) -> Result<AllocatedNum<Scalar>, SynthesisError> {
        match self {
            Elt::Num(num) => {
                let mut tmp = num.get_value().ok_or(SynthesisError::AssignmentMissing)?;
                tmp = tmp * tmp;
                let allocated =
                    AllocatedNum::alloc(&mut cs.namespace(|| "squared num"), || Ok(tmp))?;
                cs.enforce(
                    || "squaring constraint",
                    |_| num.lc(Scalar::one()),
                    |_| num.lc(Scalar::one()),
                    |lc| lc + allocated.get_variable(),
                );
                Ok(allocated)
            }
            Elt::Allocated(a) => a.square(cs),
        }
    }

    pub fn num(&self) -> num::Num<Scalar> {
        match self {
            Elt::Num(num) => num.clone(),
            Elt::Allocated(a) => a.clone().into(),
        }
    }
}


/// Circuit for Poseidon hash.
pub struct PoseidonCircuit2<'a, Scalar, A>
where
    Scalar: PrimeField,
    A: Arity<Scalar>,
{
    constants_offset: usize,
    width: usize,
    pub(crate) elements: Vec<Elt<Scalar>>,
    pub(crate) pos: usize,
    current_round: usize,
    constants: &'a PoseidonConstants<Scalar, A>,
    _w: PhantomData<A>,
}

/// PoseidonCircuit2 implementation.
impl<'a, Scalar, A> PoseidonCircuit2<'a, Scalar, A>
where
    Scalar: PrimeField,
    A: Arity<Scalar>,
{
    /// Create a new Poseidon hasher for `preimage`.
    pub fn new(elements: Vec<Elt<Scalar>>, constants: &'a PoseidonConstants<Scalar, A>) -> Self {
        let width = constants.width();

        PoseidonCircuit2 {
            constants_offset: 0,
            width,
            elements,
            pos: 1,
            current_round: 0,
            constants,
            _w: PhantomData::<A>,
        }
    }

    pub fn new_empty<CS: ConstraintSystem<Scalar>>(
        constants: &'a PoseidonConstants<Scalar, A>,
    ) -> Self {
        let elements = Self::initial_elements::<CS>();
        Self::new(elements, constants)
    }

    pub fn hash<CS: ConstraintSystem<Scalar>>(
        &mut self,
        cs: &mut CS,
        hnum: i32
    ) -> Result<Elt<Scalar>, SynthesisError> {
        self.full_round(cs.namespace(|| format!("hash num {} : first round", hnum)), true, false)?;

        for i in 1..self.constants.full_rounds / 2 {
            self.full_round(
                cs.namespace(|| format!("hash num {} : initial full round {}", hnum, i)),
                false,
                false,
            )?;
        }

        for i in 0..self.constants.partial_rounds {
            self.partial_round(cs.namespace(|| format!("hash num {} : partial round {}", hnum, i)))?;
        }

        for i in 0..(self.constants.full_rounds / 2) - 1 {
            self.full_round(
                cs.namespace(|| format!("hash num {} : final full round {}", hnum, i)),
                false,
                false,
            )?;
        }
        self.full_round(cs.namespace(|| format!("hash num {} : terminal full round", hnum)), false, true)?;

        let elt = self.elements[1].clone();
        self.reset_offsets();

        Ok(elt)
    }

    // pub fn apply_padding<CS: ConstraintSystem<Scalar>>(&mut self) {
    //     if let HashType::ConstantLength(l) = self.constants.hash_type {
    //         let final_pos = 1 + (l % self.constants.arity());

    //         assert_eq!(
    //             self.pos, final_pos,
    //             "preimage length does not match constant length required for hash"
    //         );
    //     };
    //     match self.constants.hash_type {
    //         HashType::ConstantLength(_) | HashType::Encryption => {
    //             for elt in self.elements[self.pos..].iter_mut() {
    //                 *elt = Elt::num_from_fr::<CS>(Scalar::zero());
    //             }
    //             self.pos = self.elements.len();
    //         }
    //         HashType::VariableLength => todo!(),
    //         HashType::Sponge => (),
    //         _ => (),
    //     }
    // }

    fn hash_to_allocated<CS: ConstraintSystem<Scalar>>(
        &mut self,
        mut cs: CS,
        hnum: i32
    ) -> Result<AllocatedNum<Scalar>, SynthesisError> {
        let elt = self.hash(&mut cs, hnum).unwrap();
        elt.ensure_allocated(&mut cs, true, hnum)
    }

    // fn hash_to_num<CS: ConstraintSystem<Scalar>>(
    //     &mut self,
    //     mut cs: CS,
    // ) -> Result<num::Num<Scalar>, SynthesisError> {
    //     self.hash(&mut cs).map(|elt| elt.num())
    // }

    fn full_round<CS: ConstraintSystem<Scalar>>(
        &mut self,
        mut cs: CS,
        first_round: bool,
        last_round: bool,
    ) -> Result<(), SynthesisError> {
        let mut constants_offset = self.constants_offset;

        let pre_round_keys = if first_round {
            (0..self.width)
                .map(|i| self.constants.compressed_round_constants[constants_offset + i])
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        constants_offset += pre_round_keys.len();

        let post_round_keys = if first_round || !last_round {
            (0..self.width)
                .map(|i| self.constants.compressed_round_constants[constants_offset + i])
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        constants_offset += post_round_keys.len();

        // Apply the quintic S-Box to all elements
        for i in 0..self.elements.len() {
            let pre_round_key = if first_round {
                let rk = pre_round_keys[i];
                Some(rk)
            } else {
                None
            };

            let post_round_key = if first_round || !last_round {
                let rk = post_round_keys[i];
                Some(rk)
            } else {
                None
            };

            if first_round {
                {
                    self.elements[i] = quintic_s_box_pre_add(
                        cs.namespace(|| format!("quintic s-box {}", i)),
                        &self.elements[i],
                        pre_round_key,
                        post_round_key,
                    )?;
                }
            } else {
                self.elements[i] = quintic_s_box(
                    cs.namespace(|| format!("quintic s-box {}", i)),
                    &self.elements[i],
                    post_round_key,
                )?;
            }
        }
        self.constants_offset = constants_offset;

        // Multiply the elements by the constant MDS matrix
        self.product_mds::<CS>()?;
        Ok(())
    }

    fn partial_round<CS: ConstraintSystem<Scalar>>(
        &mut self,
        mut cs: CS,
    ) -> Result<(), SynthesisError> {
        let round_key = self.constants.compressed_round_constants[self.constants_offset];
        self.constants_offset += 1;
        // Apply the quintic S-Box to the first element.
        self.elements[0] = quintic_s_box(
            cs.namespace(|| "solitary quintic s-box"),
            &self.elements[0],
            Some(round_key),
        )?;

        // Multiply the elements by the constant MDS matrix
        self.product_mds::<CS>()?;
        Ok(())
    }

    fn product_mds_m<CS: ConstraintSystem<Scalar>>(&mut self) -> Result<(), SynthesisError> {
        self.product_mds_with_matrix::<CS>(&self.constants.mds_matrices.m)
    }

    /// Set the provided elements with the result of the product between the elements and the appropriate
    /// MDS matrix.
    #[allow(clippy::collapsible_else_if)]
    fn product_mds<CS: ConstraintSystem<Scalar>>(&mut self) -> Result<(), SynthesisError> {
        let full_half = self.constants.half_full_rounds;
        let sparse_offset = full_half - 1;
        if self.current_round == sparse_offset {
            self.product_mds_with_matrix::<CS>(&self.constants.pre_sparse_matrix)?;
        } else {
            if (self.current_round > sparse_offset)
                && (self.current_round < full_half + self.constants.partial_rounds)
            {
                let index = self.current_round - sparse_offset - 1;
                let sparse_matrix = &SparseMatrix {
                    v_rest: self.constants.sparse_matrixes[index].v_rest.clone(),
                    w_hat: self.constants.sparse_matrixes[index].w_hat.clone()
                };

                self.product_mds_with_sparse_matrix::<CS>(sparse_matrix)?;
            } else {
                self.product_mds_m::<CS>()?;
            }
        };

        self.current_round += 1;
        Ok(())
    }

    #[allow(clippy::ptr_arg)]
    fn product_mds_with_matrix<CS: ConstraintSystem<Scalar>>(
        &mut self,
        matrix: &Matrix<Scalar>,
    ) -> Result<(), SynthesisError> {
        let mut result: Vec<Elt<Scalar>> = Vec::with_capacity(self.constants.width());

        for j in 0..self.constants.width() {
            let column = (0..self.constants.width())
                .map(|i| matrix[i][j])
                .collect::<Vec<_>>();

            let product = scalar_product::<Scalar, CS>(self.elements.as_slice(), &column)?;

            result.push(product);
        }

        self.elements = result;

        Ok(())
    }

    // Sparse matrix in this context means one of the form, M''.
    fn product_mds_with_sparse_matrix<CS: ConstraintSystem<Scalar>>(
        &mut self,
        matrix: &SparseMatrix<Scalar>,
    ) -> Result<(), SynthesisError> {
        let mut result: Vec<Elt<Scalar>> = Vec::with_capacity(self.constants.width());

        result.push(scalar_product::<Scalar, CS>(
            self.elements.as_slice(),
            &matrix.w_hat,
        )?);

        for j in 1..self.width {
            result.push(
                self.elements[j].clone().add(
                    self.elements[0]
                        .clone() // First row is dense.
                        .scale::<CS>(matrix.v_rest[j - 1])?, // Except for first row/column, diagonals are one.
                )?,
            );
        }

        self.elements = result;

        Ok(())
    }

    fn initial_elements<CS: ConstraintSystem<Scalar>>() -> Vec<Elt<Scalar>> {
        std::iter::repeat(Elt::num_from_fr::<CS>(Scalar::zero()))
            .take(A::to_usize() + 1)
            .collect()
    }
    pub fn reset<CS: ConstraintSystem<Scalar>>(&mut self) {
        self.reset_offsets();
        self.elements = Self::initial_elements::<CS>();
    }

    pub fn reset_offsets(&mut self) {
        self.constants_offset = 0;
        self.current_round = 0;
        self.pos = 1;
    }

    // fn debug(&self) {
    //     let element_frs: Vec<_> = self.elements.iter().map(|n| n.val()).collect::<Vec<_>>();
    //     dbg!(element_frs, self.constants_offset);
    // }
}

/// Create circuit for Poseidon hash, returning an unallocated `Num` at the cost of one constraint.
pub fn poseidon_hash_allocated<CS, Scalar, A>(
    cs: CS,
    preimage: Vec<AllocatedNum<Scalar>>,
    constants: &PoseidonConstants<Scalar, A>,
    hnum: i32
) -> Result<AllocatedNum<Scalar>, SynthesisError>
where
    CS: ConstraintSystem<Scalar>,
    Scalar: PrimeField,
    A: Arity<Scalar>,
{
    let arity = A::to_usize();
    let tag_element = Elt::num_from_fr::<CS>(constants.domain_tag);
    let mut elements = Vec::with_capacity(arity + 1);
    elements.push(tag_element);
    elements.extend(preimage.into_iter().map(Elt::Allocated));

    // if let HashType::ConstantLength(length) = constants.hash_type {
    //     assert!(length <= arity, "illegal length: constants are malformed");
    //     // Add zero-padding.
    //     for i in 0..(arity - length) {
    //         let allocated = AllocatedNum::alloc(cs.namespace(|| format!("padding {}", i)), || {
    //             Ok(Scalar::zero())
    //         })?;
    //         let elt = Elt::Allocated(allocated);
    //         elements.push(elt);
    //     }
    // }

    let mut p = PoseidonCircuit2::new(elements, constants);

    p.hash_to_allocated(cs, hnum)
}


/// Compute l^5 and enforce constraint. If round_key is supplied, add it to result.
fn quintic_s_box<CS: ConstraintSystem<Scalar>, Scalar: PrimeField>(
    mut cs: CS,
    l: &Elt<Scalar>,
    post_round_key: Option<Scalar>,
) -> Result<Elt<Scalar>, SynthesisError> {
    // If round_key was supplied, add it after all exponentiation.
    let l2 = l.square(cs.namespace(|| "l^2"))?;
    let l4 = l2.square(cs.namespace(|| "l^4"))?;
    let l5 = mul_sum(
        cs.namespace(|| "(l4 * l) + rk)"),
        &l4,
        l,
        None,
        post_round_key,
        true,
    );

    Ok(Elt::Allocated(l5?))
}

/// Compute l^5 and enforce constraint. If round_key is supplied, add it to l first.
fn quintic_s_box_pre_add<CS: ConstraintSystem<Scalar>, Scalar: PrimeField>(
    mut cs: CS,
    l: &Elt<Scalar>,
    pre_round_key: Option<Scalar>,
    post_round_key: Option<Scalar>,
) -> Result<Elt<Scalar>, SynthesisError> {
    if let (Some(pre_round_key), Some(post_round_key)) = (pre_round_key, post_round_key) {
        // If round_key was supplied, add it to l before squaring.
        let l2 = square_sum(cs.namespace(|| "(l+rk)^2"), pre_round_key, l, true)?;
        let l4 = l2.square(cs.namespace(|| "l^4"))?;
        let l5 = mul_sum(
            cs.namespace(|| "l4 * (l + rk)"),
            &l4,
            l,
            Some(pre_round_key),
            Some(post_round_key),
            true,
        );

        Ok(Elt::Allocated(l5?))
    } else {
        panic!("pre_round_key and post_round_key must both be provided.");
    }
}


/// Calculates square of sum and enforces that constraint.
pub fn square_sum<CS: ConstraintSystem<Scalar>, Scalar: PrimeField>(
    mut cs: CS,
    to_add: Scalar,
    elt: &Elt<Scalar>,
    enforce: bool,
) -> Result<AllocatedNum<Scalar>, SynthesisError>
where
    CS: ConstraintSystem<Scalar>,
{
    let res = AllocatedNum::alloc(cs.namespace(|| "squared sum"), || {
        let mut tmp = elt.val().ok_or(SynthesisError::AssignmentMissing)?;
        tmp.add_assign(&to_add);
        tmp = tmp.square();

        Ok(tmp)
    })?;

    if enforce {
        cs.enforce(
            || "squared sum constraint",
            |_| elt.lc() + (to_add, CS::one()),
            |_| elt.lc() + (to_add, CS::one()),
            |lc| lc + res.get_variable(),
        );
    }
    Ok(res)
}

/// Calculates (a * (pre_add + b)) + post_add — and enforces that constraint.
#[allow(clippy::collapsible_else_if)]
pub fn mul_sum<CS: ConstraintSystem<Scalar>, Scalar: PrimeField>(
    mut cs: CS,
    a: &AllocatedNum<Scalar>,
    b: &Elt<Scalar>,
    pre_add: Option<Scalar>,
    post_add: Option<Scalar>,
    enforce: bool,
) -> Result<AllocatedNum<Scalar>, SynthesisError>
where
    CS: ConstraintSystem<Scalar>,
{
    let res = AllocatedNum::alloc(cs.namespace(|| "mul_sum"), || {
        let mut tmp = b.val().ok_or(SynthesisError::AssignmentMissing)?;
        if let Some(x) = pre_add {
            tmp.add_assign(&x);
        }
        tmp.mul_assign(&a.get_value().ok_or(SynthesisError::AssignmentMissing)?);
        if let Some(x) = post_add {
            tmp.add_assign(&x);
        }

        Ok(tmp)
    })?;

    if enforce {
        if let Some(x) = post_add {
            let neg = -x;

            if let Some(pre) = pre_add {
                cs.enforce(
                    || "mul sum constraint pre-post-add",
                    |_| b.lc() + (pre, CS::one()),
                    |lc| lc + a.get_variable(),
                    |lc| lc + res.get_variable() + (neg, CS::one()),
                );
            } else {
                cs.enforce(
                    || "mul sum constraint post-add",
                    |_| b.lc(),
                    |lc| lc + a.get_variable(),
                    |lc| lc + res.get_variable() + (neg, CS::one()),
                );
            }
        } else {
            if let Some(pre) = pre_add {
                cs.enforce(
                    || "mul sum constraint pre-add",
                    |_| b.lc() + (pre, CS::one()),
                    |lc| lc + a.get_variable(),
                    |lc| lc + res.get_variable(),
                );
            } else {
                cs.enforce(
                    || "mul sum constraint",
                    |_| b.lc(),
                    |lc| lc + a.get_variable(),
                    |lc| lc + res.get_variable(),
                );
            }
        }
    }
    Ok(res)
}

/// Calculates a * (b + to_add) — and enforces that constraint.
pub fn mul_pre_sum<CS: ConstraintSystem<Scalar>, Scalar: PrimeField>(
    mut cs: CS,
    a: &AllocatedNum<Scalar>,
    b: &AllocatedNum<Scalar>,
    to_add: Scalar,
    enforce: bool,
) -> Result<AllocatedNum<Scalar>, SynthesisError>
where
    CS: ConstraintSystem<Scalar>,
{
    let res = AllocatedNum::alloc(cs.namespace(|| "mul_sum"), || {
        let mut tmp = b.get_value().ok_or(SynthesisError::AssignmentMissing)?;
        tmp.add_assign(&to_add);
        tmp.mul_assign(&a.get_value().ok_or(SynthesisError::AssignmentMissing)?);

        Ok(tmp)
    })?;

    if enforce {
        cs.enforce(
            || "mul sum constraint",
            |lc| lc + b.get_variable() + (to_add, CS::one()),
            |lc| lc + a.get_variable(),
            |lc| lc + res.get_variable(),
        );
    }
    Ok(res)
}

// fn scalar_product_with_add<Scalar: PrimeField, CS: ConstraintSystem<Scalar>>(
//     elts: &[Elt<Scalar>],
//     scalars: &[Scalar],
//     to_add: Scalar,
// ) -> Result<Elt<Scalar>, SynthesisError> {
//     let tmp = scalar_product::<Scalar, CS>(elts, scalars)?;
//     let tmp2 = tmp.add(Elt::<Scalar>::num_from_fr::<CS>(to_add))?;

//     Ok(tmp2)
// }

fn scalar_product<Scalar: PrimeField, CS: ConstraintSystem<Scalar>>(
    elts: &[Elt<Scalar>],
    scalars: &[Scalar],
) -> Result<Elt<Scalar>, SynthesisError> {
    elts.iter()
        .zip(scalars)
        .try_fold(Elt::Num(num::Num::zero()), |acc, (elt, &scalar)| {
            acc.add(elt.clone().scale::<CS>(scalar)?)
        })
}

// struct PoseidonCircuit<F: PrimeField, A: Arity<F>> {
//     pub xl: F,
//     pub xr: F,
//     pub node: F,
//     pub params: PoseidonConstants<F, A>
// }

// impl<F: PrimeField, A: Arity<F>> Circuit<F> for PoseidonCircuit<F, A> {
//     fn synthesize<CS: ConstraintSystem<F>>(self, cs: &mut CS) -> Result<(), SynthesisError> {
//         let xl = AllocatedNum::alloc(cs.namespace(|| "preimage xl"), || Ok(self.xl))?;
//         let xr = AllocatedNum::alloc(cs.namespace(|| "preimage xr"), || Ok(self.xr))?;
//         let node = AllocatedNum::alloc_input(cs.namespace(|| "hash node"), || Ok(self.node))?;

//         let i = 0;

//         let calc_node = poseidon_hash_allocated(&mut *cs, vec![xl, xr], &self.params, i)?;

//         cs.enforce(
//             || "node = calc_node", 
//             |lc| lc + node.get_variable(), 
//             |lc| lc + CS::one(), 
//             |lc| lc + calc_node.get_variable()
//         );

//         Ok(())
        
//     }
// }

// #[cfg(test)]
// mod tests {
//     use bellperson::Circuit;
//     // use neptune::circuit2::poseidon_hash_allocated;
//     use pasta_curves::pallas::Base as Fp;
    
//     use neptune::poseidon::{PoseidonConstants, Poseidon};
//     // use neptune::circuit::poseidon_hash;
//     use generic_array::typenum::U2;
//     use bellperson::{util_cs::test_cs::TestConstraintSystem ,
//         // ConstraintSystem, gadgets::num::AllocatedNum
//     };

//     use crate::poseidon_hash::constraints::PoseidonCircuit;

//     #[test]
//     fn test_poseidon_circuit() {
//         // Need to think of a way to test these output values
//         let leaf_present  = Fp::one();
		
//         let node_hash_params = PoseidonConstants::<Fp, U2>::new();
//         let mut node_hasher = Poseidon::new_with_preimage(&[leaf_present.double(), leaf_present],&node_hash_params);
//         let hash_value = node_hasher.hash();
//         // println!("the hash value is {:?}", hash_value);

//         let circ = PoseidonCircuit {
//             xl: leaf_present.double(),
//             xr: leaf_present,
//             node: hash_value,
//             params: node_hash_params
//         };

//         let mut cs = TestConstraintSystem::<Fp>::new();

//         println!("the number of constraints are {}", cs.num_constraints());

        
//         assert!(!circ.synthesize(& mut cs).is_err());

//         println!("the number of constraints are {}", cs.num_constraints());

//         assert!(cs.is_satisfied());

//         println!("constraint satisfied is {:?}", cs.is_satisfied());

//         // println!("the inputs are {:?}", cs.get_inputs());

//         // let xl = AllocatedNum::alloc(cs.namespace(|| "preimage xl"), || Ok(leaf_present.double())).unwrap();
//         // let xr = AllocatedNum::alloc(cs.namespace(|| "preimage xr"), || Ok(leaf_present)).unwrap();
//         // let calc_node = poseidon_hash_allocated(&mut cs, vec![xl, xr], &node_hash_params).unwrap();
//         // println!("calc_node hash is {:?}",calc_node.get_value().unwrap());

//     }
// }