use crate::poly::{Poly, compare_monomials};
use crate::rational::Rational;
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

pub fn reduced_groebner_basis(generators: &[Poly]) -> Vec<Poly> {
    let mut basis = generators
        .iter()
        .filter(|poly| !poly.is_zero())
        .map(Poly::monic)
        .collect::<Vec<_>>();

    if basis.iter().any(Poly::is_nonzero_constant) {
        return vec![constant_one(&basis)];
    }
    if basis.is_empty() {
        return basis;
    }

    let mut pairs = VecDeque::new();
    for left in 0..basis.len() {
        for right in left + 1..basis.len() {
            pairs.push_back((left, right));
        }
    }

    while let Some((left, right)) = pairs.pop_front() {
        if relatively_prime_leading_monomials(&basis[left], &basis[right]) {
            continue;
        }

        let s_poly = s_polynomial(&basis[left], &basis[right]);
        let remainder = normal_form(&s_poly, &basis);
        if remainder.is_zero() {
            continue;
        }

        let next = remainder.monic();
        if next.is_nonzero_constant() {
            return vec![Poly::constant(&next.vars, Rational::one())];
        }

        let next_index = basis.len();
        for index in 0..next_index {
            pairs.push_back((index, next_index));
        }
        basis.push(next);
    }

    reduce_basis(basis)
}

pub fn normal_form(poly: &Poly, basis: &[Poly]) -> Poly {
    let mut working = poly.clone();
    let mut remainder = Poly::zero(&poly.vars);

    while let Some((leading_exp, leading_coeff)) = working.leading_term() {
        let mut reduced = false;
        for divisor in basis {
            if divisor.is_zero() {
                continue;
            }
            let Some((divisor_exp, divisor_coeff)) = divisor.leading_term() else {
                continue;
            };
            if monomial_divides(&divisor_exp, &leading_exp) {
                let quotient_exp = subtract_exponents(&leading_exp, &divisor_exp);
                let quotient_coeff = leading_coeff.clone() / divisor_coeff;
                working = working.sub(&divisor.mul_monomial(quotient_coeff, &quotient_exp));
                reduced = true;
                break;
            }
        }

        if !reduced {
            let leading = Poly::monomial(&poly.vars, leading_coeff, leading_exp);
            remainder = remainder.add(&leading);
            working = working.sub(&leading);
        }
    }

    remainder
}

pub fn decompose_ideal(generators: &[Poly]) -> Vec<Vec<Poly>> {
    let mut raw_components = Vec::new();
    let mut seen_nodes = BTreeSet::new();
    decompose_recursive(generators, &mut raw_components, &mut seen_nodes);
    remove_redundant_components(raw_components)
}

pub fn ideal_dimension(groebner_basis: &[Poly], var_count: usize) -> Option<usize> {
    if groebner_basis.iter().any(Poly::is_nonzero_constant) {
        return None;
    }

    let leading_supports = groebner_basis
        .iter()
        .filter_map(|poly| poly.leading_term().map(|(exp, _)| support_mask(&exp)))
        .filter(|mask| *mask != 0)
        .collect::<Vec<_>>();

    let mut best = 0usize;
    for mask in 0usize..(1usize << var_count) {
        let is_independent = leading_supports
            .iter()
            .all(|leading_support| leading_support & !mask != 0);
        if is_independent {
            best = best.max(mask.count_ones() as usize);
        }
    }
    Some(best)
}

pub fn ideal_contains(groebner_basis: &[Poly], poly: &Poly) -> bool {
    normal_form(poly, groebner_basis).is_zero()
}

pub fn ideal_contains_all(left_groebner_basis: &[Poly], right_generators: &[Poly]) -> bool {
    right_generators
        .iter()
        .all(|generator| ideal_contains(left_groebner_basis, generator))
}

pub fn groebner_key(groebner_basis: &[Poly]) -> String {
    groebner_basis
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join("; ")
}

pub fn generic_multiplicity(hypersurface: &Poly, component_groebner_basis: &[Poly]) -> usize {
    if hypersurface.is_zero()
        || component_groebner_basis
            .iter()
            .any(Poly::is_nonzero_constant)
    {
        return 0;
    }

    let max_power = hypersurface.total_degree().max(1);
    let mut multiplicity = 0usize;
    for power in 1..=max_power {
        let power_generators = ideal_power_generators(component_groebner_basis, power);
        if power_generators.is_empty() {
            break;
        }
        let power_basis = reduced_groebner_basis(&power_generators);
        if ideal_contains(&power_basis, hypersurface) {
            multiplicity = power;
        } else {
            break;
        }
    }
    multiplicity
}

fn decompose_recursive(
    generators: &[Poly],
    components: &mut Vec<Vec<Poly>>,
    seen_nodes: &mut BTreeSet<String>,
) {
    let groebner_basis = reduced_groebner_basis(generators);
    if groebner_basis.iter().any(Poly::is_nonzero_constant) {
        return;
    }

    let key = groebner_key(&groebner_basis);
    if !seen_nodes.insert(key.clone()) {
        return;
    }

    if let Some(factors) = choose_variety_split(&groebner_basis) {
        let mut branched = false;
        for factor in factors {
            if ideal_contains(&groebner_basis, &factor) {
                continue;
            }

            let mut next_generators = groebner_basis.clone();
            next_generators.push(factor);
            let next_basis = reduced_groebner_basis(&next_generators);
            if groebner_key(&next_basis) == key {
                continue;
            }
            decompose_recursive(&next_basis, components, seen_nodes);
            branched = true;
        }

        if branched {
            return;
        }
    }

    components.push(groebner_basis);
}

fn remove_redundant_components(components: Vec<Vec<Poly>>) -> Vec<Vec<Poly>> {
    let mut unique = BTreeMap::new();
    for component in components {
        unique.entry(groebner_key(&component)).or_insert(component);
    }
    let components = unique.into_values().collect::<Vec<_>>();

    let mut keep = vec![true; components.len()];
    for left in 0..components.len() {
        for right in 0..components.len() {
            if left == right {
                continue;
            }
            let left_contains_right = ideal_contains_all(&components[left], &components[right]);
            let right_contains_left = ideal_contains_all(&components[right], &components[left]);
            if left_contains_right && !right_contains_left {
                keep[left] = false;
                break;
            }
        }
    }

    let mut minimal = components
        .into_iter()
        .enumerate()
        .filter_map(|(index, component)| keep[index].then_some(component))
        .collect::<Vec<_>>();
    minimal.sort_by(|left, right| {
        let left_dim = ideal_dimension(left, left.first().map_or(0, |poly| poly.vars.len()));
        let right_dim = ideal_dimension(right, right.first().map_or(0, |poly| poly.vars.len()));
        right_dim
            .cmp(&left_dim)
            .then_with(|| groebner_key(left).cmp(&groebner_key(right)))
    });
    minimal
}

fn reduce_basis(mut basis: Vec<Poly>) -> Vec<Poly> {
    basis.retain(|poly| !poly.is_zero());
    if basis.iter().any(Poly::is_nonzero_constant) {
        return vec![constant_one(&basis)];
    }
    let mut unique_basis = BTreeMap::new();
    for poly in basis {
        unique_basis.entry(poly.to_string()).or_insert(poly);
    }
    let basis = unique_basis.into_values().collect::<Vec<_>>();

    let mut reduced = Vec::new();
    for index in 0..basis.len() {
        let others = basis
            .iter()
            .enumerate()
            .filter_map(|(other_index, poly)| (other_index != index).then_some(poly.clone()))
            .collect::<Vec<_>>();
        let remainder = normal_form(&basis[index], &others);
        if !remainder.is_zero() {
            reduced.push(remainder.monic());
        }
    }

    let mut unique = BTreeMap::new();
    for poly in reduced {
        unique.entry(poly.to_string()).or_insert(poly);
    }
    let mut reduced = unique.into_values().collect::<Vec<_>>();
    reduced.sort_by(|left, right| {
        let left_exp = left.leading_term().map(|(exp, _)| exp).unwrap_or_default();
        let right_exp = right.leading_term().map(|(exp, _)| exp).unwrap_or_default();
        compare_monomials(&left_exp, &right_exp)
            .reverse()
            .then_with(|| left.to_string().cmp(&right.to_string()))
    });
    reduced
}

fn s_polynomial(left: &Poly, right: &Poly) -> Poly {
    let (left_exp, left_coeff) = left.leading_term().expect("nonzero left polynomial");
    let (right_exp, right_coeff) = right.leading_term().expect("nonzero right polynomial");
    let lcm = lcm_exponents(&left_exp, &right_exp);
    let left_multiplier = subtract_exponents(&lcm, &left_exp);
    let right_multiplier = subtract_exponents(&lcm, &right_exp);

    left.mul_monomial(Rational::one() / left_coeff, &left_multiplier)
        .sub(&right.mul_monomial(Rational::one() / right_coeff, &right_multiplier))
}

fn relatively_prime_leading_monomials(left: &Poly, right: &Poly) -> bool {
    let Some((left_exp, _)) = left.leading_term() else {
        return false;
    };
    let Some((right_exp, _)) = right.leading_term() else {
        return false;
    };
    left_exp
        .iter()
        .zip(right_exp.iter())
        .all(|(left_power, right_power)| *left_power == 0 || *right_power == 0)
}

fn choose_variety_split(groebner_basis: &[Poly]) -> Option<Vec<Poly>> {
    for poly in groebner_basis {
        let factors = variety_factors(poly);
        if factors.len() > 1
            || factors
                .first()
                .is_some_and(|factor| !same_up_to_unit(poly, factor))
        {
            return Some(factors);
        }
    }
    None
}

fn variety_factors(poly: &Poly) -> Vec<Poly> {
    if poly.is_zero() || poly.is_nonzero_constant() {
        return vec![poly.clone()];
    }

    if poly.terms.len() == 1 {
        let exp = poly.terms.keys().next().expect("single term");
        let mut factors = exp
            .iter()
            .enumerate()
            .filter_map(|(index, power)| (*power > 0).then(|| Poly::var(&poly.vars, index)))
            .collect::<Vec<_>>();
        if factors.is_empty() {
            factors.push(poly.clone());
        }
        return unique_factors(factors);
    }

    let common_exp = common_monomial_exponent(poly);
    if common_exp.iter().any(|power| *power > 0) {
        let mut factors = common_exp
            .iter()
            .enumerate()
            .filter_map(|(index, power)| (*power > 0).then(|| Poly::var(&poly.vars, index)))
            .collect::<Vec<_>>();
        let residual_terms = poly.terms.iter().map(|(exp, coeff)| {
            let residual_exp = subtract_exponents(exp, &common_exp);
            (residual_exp, coeff.clone())
        });
        let residual = Poly::from_terms(&poly.vars, residual_terms).monic();
        if !residual.is_nonzero_constant() {
            factors.push(residual);
        }
        return unique_factors(factors);
    }

    if let Some(factors) = binomial_difference_of_squares(poly) {
        return unique_factors(factors);
    }

    if let Some((var_index, coeffs)) = univariate_coefficients(poly) {
        let (factors, changed) = factor_univariate_by_rational_roots(&poly.vars, var_index, coeffs);
        if changed {
            return unique_factors(factors);
        }
    }

    vec![poly.monic()]
}

fn binomial_difference_of_squares(poly: &Poly) -> Option<Vec<Poly>> {
    if poly.terms.len() != 2 {
        return None;
    }

    let terms = poly.terms.iter().collect::<Vec<_>>();
    let (left_exp, left_coeff, right_exp, right_coeff) =
        if terms[0].1.is_negative() && !terms[1].1.is_negative() {
            (
                terms[1].0,
                terms[1].1.clone(),
                terms[0].0,
                terms[0].1.clone(),
            )
        } else {
            (
                terms[0].0,
                terms[0].1.clone(),
                terms[1].0,
                terms[1].1.clone(),
            )
        };

    if left_coeff.is_negative() || !right_coeff.is_negative() {
        return None;
    }
    if !left_exp.iter().all(|power| power % 2 == 0) || !right_exp.iter().all(|power| power % 2 == 0)
    {
        return None;
    }

    let left_sqrt = rational_sqrt(left_coeff)?;
    let right_sqrt = rational_sqrt(-right_coeff)?;
    let left_root_exp = left_exp.iter().map(|power| power / 2).collect::<Vec<_>>();
    let right_root_exp = right_exp.iter().map(|power| power / 2).collect::<Vec<_>>();
    let left_root = Poly::monomial(&poly.vars, left_sqrt, left_root_exp);
    let right_root = Poly::monomial(&poly.vars, right_sqrt, right_root_exp);

    Some(vec![left_root.sub(&right_root), left_root.add(&right_root)])
}

fn factor_univariate_by_rational_roots(
    vars: &[String],
    var_index: usize,
    mut coeffs: Vec<Rational>,
) -> (Vec<Poly>, bool) {
    trim_trailing_zeros(&mut coeffs);
    let original_degree = coeffs.len().saturating_sub(1);
    let mut factors = Vec::new();
    let mut changed = false;

    while coeffs.len() > 1 {
        let Some(root) = find_rational_root(&coeffs) else {
            break;
        };
        factors.push(linear_factor(vars, var_index, root.clone()));
        coeffs = divide_by_linear_root(&coeffs, root);
        trim_trailing_zeros(&mut coeffs);
        changed = true;
    }

    if changed && coeffs.len() > 1 {
        factors.push(poly_from_univariate_coefficients(vars, var_index, &coeffs).monic());
    }

    if !changed && original_degree > 0 {
        factors.push(poly_from_univariate_coefficients(vars, var_index, &coeffs).monic());
    }

    (factors, changed)
}

fn find_rational_root(coeffs: &[Rational]) -> Option<Rational> {
    if coeffs.first().is_some_and(|coeff| coeff.is_zero()) {
        return Some(Rational::zero());
    }

    let int_coeffs = clear_denominators(coeffs)?;
    let constant = checked_abs_i128(*int_coeffs.first()?)?;
    let leading = checked_abs_i128(*int_coeffs.last()?)?;
    if constant == 0 || leading == 0 {
        return None;
    }

    let numerators = divisors(constant);
    let denominators = divisors(leading);
    let mut candidates = BTreeSet::new();
    for numerator in numerators {
        for denominator in &denominators {
            candidates.insert(Rational::new(numerator, *denominator));
            candidates.insert(Rational::new(-numerator, *denominator));
        }
    }

    candidates
        .into_iter()
        .find(|candidate| evaluate_univariate(coeffs, candidate.clone()).is_zero())
}

fn clear_denominators(coeffs: &[Rational]) -> Option<Vec<i128>> {
    let mut common_denominator = 1i128;
    for coeff in coeffs {
        common_denominator = lcm_i128(common_denominator, coeff.denominator_i128()?)?;
    }
    Some(
        coeffs
            .iter()
            .map(|coeff| {
                coeff
                    .numerator_i128()?
                    .checked_mul(common_denominator.checked_div(coeff.denominator_i128()?)?)
            })
            .collect::<Option<Vec<_>>>()?,
    )
}

fn divide_by_linear_root(coeffs: &[Rational], root: Rational) -> Vec<Rational> {
    let degree = coeffs.len() - 1;
    let mut quotient = vec![Rational::zero(); degree];
    quotient[degree - 1] = coeffs[degree].clone();
    for index in (1..degree).rev() {
        quotient[index - 1] = coeffs[index].clone() + root.clone() * quotient[index].clone();
    }
    quotient
}

fn evaluate_univariate(coeffs: &[Rational], value: Rational) -> Rational {
    coeffs.iter().rev().fold(Rational::zero(), |acc, coeff| {
        acc * value.clone() + coeff.clone()
    })
}

fn univariate_coefficients(poly: &Poly) -> Option<(usize, Vec<Rational>)> {
    let used_vars = used_variable_indices(poly);
    if used_vars.len() != 1 {
        return None;
    }

    let var_index = used_vars[0];
    let degree = poly
        .terms
        .keys()
        .map(|exp| exp[var_index])
        .max()
        .unwrap_or(0);
    let mut coeffs = vec![Rational::zero(); degree + 1];
    for (exp, coeff) in &poly.terms {
        coeffs[exp[var_index]] = coeffs[exp[var_index]].clone() + coeff.clone();
    }
    Some((var_index, coeffs))
}

fn poly_from_univariate_coefficients(
    vars: &[String],
    var_index: usize,
    coeffs: &[Rational],
) -> Poly {
    let terms = coeffs.iter().enumerate().filter_map(|(power, coeff)| {
        if coeff.is_zero() {
            None
        } else {
            let mut exp = vec![0; vars.len()];
            exp[var_index] = power;
            Some((exp, coeff.clone()))
        }
    });
    Poly::from_terms(vars, terms)
}

fn linear_factor(vars: &[String], var_index: usize, root: Rational) -> Poly {
    let var = Poly::var(vars, var_index);
    var.sub(&Poly::constant(vars, root)).monic()
}

fn ideal_power_generators(generators: &[Poly], power: usize) -> Vec<Poly> {
    if power == 0 || generators.is_empty() {
        return Vec::new();
    }

    let vars = generators[0].vars.clone();
    let mut products = Vec::new();
    ideal_power_products_recursive(
        generators,
        power,
        0,
        Poly::constant(&vars, Rational::one()),
        &mut products,
    );
    products
}

fn ideal_power_products_recursive(
    generators: &[Poly],
    remaining_power: usize,
    start_index: usize,
    current: Poly,
    products: &mut Vec<Poly>,
) {
    if remaining_power == 0 {
        if !current.is_zero() {
            products.push(current.monic());
        }
        return;
    }

    for index in start_index..generators.len() {
        let next = current.mul(&generators[index]);
        ideal_power_products_recursive(generators, remaining_power - 1, index, next, products);
    }
}

fn unique_factors(factors: Vec<Poly>) -> Vec<Poly> {
    let mut unique = BTreeMap::new();
    for factor in factors {
        if factor.is_zero() || factor.is_nonzero_constant() {
            continue;
        }
        let factor = factor.monic();
        unique.entry(factor.to_string()).or_insert(factor);
    }
    unique.into_values().collect()
}

fn same_up_to_unit(left: &Poly, right: &Poly) -> bool {
    left.monic() == right.monic()
}

fn common_monomial_exponent(poly: &Poly) -> Vec<usize> {
    let mut exponents = poly
        .terms
        .keys()
        .next()
        .cloned()
        .unwrap_or_else(|| vec![0; poly.vars.len()]);
    for exp in poly.terms.keys().skip(1) {
        for (common, power) in exponents.iter_mut().zip(exp.iter()) {
            *common = (*common).min(*power);
        }
    }
    exponents
}

fn used_variable_indices(poly: &Poly) -> Vec<usize> {
    (0..poly.vars.len())
        .filter(|index| poly.terms.keys().any(|exp| exp[*index] > 0))
        .collect()
}

fn monomial_divides(divisor: &[usize], multiple: &[usize]) -> bool {
    divisor
        .iter()
        .zip(multiple.iter())
        .all(|(left, right)| left <= right)
}

fn lcm_exponents(left: &[usize], right: &[usize]) -> Vec<usize> {
    left.iter()
        .zip(right.iter())
        .map(|(left, right)| (*left).max(*right))
        .collect()
}

fn subtract_exponents(left: &[usize], right: &[usize]) -> Vec<usize> {
    left.iter()
        .zip(right.iter())
        .map(|(left, right)| left - right)
        .collect()
}

fn support_mask(exp: &[usize]) -> usize {
    exp.iter().enumerate().fold(0usize, |mask, (index, power)| {
        if *power > 0 {
            mask | (1usize << index)
        } else {
            mask
        }
    })
}

fn constant_one(basis: &[Poly]) -> Poly {
    let vars = basis
        .first()
        .map_or_else(Vec::new, |poly| poly.vars.clone());
    Poly::constant(&vars, Rational::one())
}

fn trim_trailing_zeros(coeffs: &mut Vec<Rational>) {
    while coeffs.len() > 1 && coeffs.last().is_some_and(|coeff| coeff.is_zero()) {
        coeffs.pop();
    }
}

fn divisors(value: i128) -> Vec<i128> {
    if value == 0 {
        return vec![0];
    }

    let mut result = Vec::new();
    let mut candidate = 1i128;
    while candidate <= value / candidate {
        if value % candidate == 0 {
            result.push(candidate);
            let paired = value / candidate;
            if candidate != paired {
                result.push(paired);
            }
        }
        candidate += 1;
    }
    result.sort();
    result
}

fn lcm_i128(left: i128, right: i128) -> Option<i128> {
    if left == 0 || right == 0 {
        return Some(0);
    }
    let gcd = gcd_i128(checked_abs_i128(left)?, checked_abs_i128(right)?);
    checked_abs_i128(left.checked_div(gcd)?.checked_mul(right)?)
}

fn gcd_i128(mut left: i128, mut right: i128) -> i128 {
    while right != 0 {
        let rem = left % right;
        left = right;
        right = rem;
    }
    if left == 0 {
        1
    } else {
        checked_abs_i128(left).unwrap_or(left)
    }
}

fn checked_abs_i128(value: i128) -> Option<i128> {
    value.checked_abs()
}

fn rational_sqrt(value: Rational) -> Option<Rational> {
    if value.is_negative() {
        return None;
    }
    Some(Rational::new(
        integer_sqrt_exact(&value.num)?,
        integer_sqrt_exact(&value.den)?,
    ))
}

fn integer_sqrt_exact(value: &BigInt) -> Option<BigInt> {
    if value.is_negative() {
        return None;
    }

    let mut low = BigInt::zero();
    let mut high = value.clone();
    let two = BigInt::from(2u8);
    while low <= high {
        let mid = (&low + &high) / &two;
        let square = &mid * &mid;
        match square.cmp(&value) {
            Ordering::Equal => return Some(mid),
            Ordering::Less => low = mid + BigInt::one(),
            Ordering::Greater => {
                if mid.is_zero() {
                    break;
                }
                high = mid - BigInt::one();
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::{decompose_ideal, generic_multiplicity, ideal_dimension, reduced_groebner_basis};
    use crate::poly::Poly;
    use crate::rational::Rational;

    fn vars() -> Vec<String> {
        vec!["x".to_string(), "y".to_string(), "z".to_string()]
    }

    #[test]
    fn computes_reduced_groebner_basis() {
        let vars = vars();
        let x = Poly::var(&vars, 0);
        let y = Poly::var(&vars, 1);
        let basis = reduced_groebner_basis(&[x.mul(&y), y.pow(2).sub(&x)]);

        let rendered = basis.iter().map(ToString::to_string).collect::<Vec<_>>();
        assert!(rendered.iter().any(|poly| poly == "x^2"));
        assert!(rendered.iter().any(|poly| poly == "x*y"));
        assert!(rendered.iter().any(|poly| poly == "-x + y^2"));
    }

    #[test]
    fn keeps_one_copy_of_duplicate_essential_generators() {
        let vars = vars();
        let x = Poly::var(&vars, 0);
        let y = Poly::var(&vars, 1);
        let basis = reduced_groebner_basis(&[x.clone(), y.clone(), x]);
        let rendered = basis.iter().map(ToString::to_string).collect::<Vec<_>>();

        assert_eq!(basis.len(), 2);
        assert!(rendered.contains(&"x".to_string()));
        assert!(rendered.contains(&"y".to_string()));
    }

    #[test]
    fn decomposes_coordinate_axes() {
        let vars = vars();
        let x = Poly::var(&vars, 0);
        let y = Poly::var(&vars, 1);
        let z = Poly::var(&vars, 2);
        let components = decompose_ideal(&[x.mul(&y), x.mul(&z), y.mul(&z)]);
        let rendered = components
            .iter()
            .map(|basis| {
                basis
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .collect::<Vec<_>>();

        assert_eq!(components.len(), 3);
        assert!(rendered.iter().any(|basis| basis == "x,y"));
        assert!(rendered.iter().any(|basis| basis == "x,z"));
        assert!(rendered.iter().any(|basis| basis == "y,z"));
        assert!(
            components
                .iter()
                .all(|basis| ideal_dimension(basis, vars.len()) == Some(1))
        );
    }

    #[test]
    fn splits_univariate_rational_roots() {
        let vars = vars();
        let x = Poly::var(&vars, 0);
        let components = decompose_ideal(&[x.pow(2).sub(&Poly::constant(&vars, Rational::one()))]);
        let rendered = components
            .iter()
            .map(|basis| basis[0].to_string())
            .collect::<Vec<_>>();

        assert_eq!(components.len(), 2);
        assert!(rendered.contains(&"x - 1".to_string()));
        assert!(rendered.contains(&"x + 1".to_string()));
    }

    #[test]
    fn splits_binomial_difference_of_squares() {
        let vars = vars();
        let x = Poly::var(&vars, 0);
        let y = Poly::var(&vars, 1);
        let components = decompose_ideal(&[x.pow(2).sub(&y.pow(2))]);
        let rendered = components
            .iter()
            .map(|basis| basis[0].to_string())
            .collect::<Vec<_>>();

        assert_eq!(components.len(), 2);
        assert!(rendered.contains(&"x - y".to_string()));
        assert!(rendered.contains(&"x + y".to_string()));
    }

    #[test]
    fn computes_generic_multiplicity_from_ideal_powers() {
        let vars = vars();
        let x = Poly::var(&vars, 0);
        let y = Poly::var(&vars, 1);
        let z = Poly::var(&vars, 2);
        let hypersurface = x.pow(2).mul(&z).add(&y.pow(2));
        let component = reduced_groebner_basis(&[x, y]);

        assert_eq!(generic_multiplicity(&hypersurface, &component), 2);
    }
}
