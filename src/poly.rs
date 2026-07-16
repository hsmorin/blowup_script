use crate::rational::Rational;
use num_bigint::BigInt;
use std::collections::BTreeMap;
use std::fmt;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Poly {
    pub vars: Vec<String>,
    pub terms: BTreeMap<Vec<usize>, Rational>,
}

impl Poly {
    pub fn zero(vars: &[String]) -> Self {
        Self {
            vars: vars.to_vec(),
            terms: BTreeMap::new(),
        }
    }

    pub fn constant(vars: &[String], coeff: Rational) -> Self {
        if coeff.is_zero() {
            Self::zero(vars)
        } else {
            let mut terms = BTreeMap::new();
            terms.insert(vec![0; vars.len()], coeff);
            Self {
                vars: vars.to_vec(),
                terms,
            }
        }
    }

    pub fn var(vars: &[String], index: usize) -> Self {
        let mut exp = vec![0; vars.len()];
        exp[index] = 1;
        Self::monomial(vars, Rational::one(), exp)
    }

    pub fn monomial(vars: &[String], coeff: Rational, exponents: Vec<usize>) -> Self {
        assert_eq!(vars.len(), exponents.len());
        if coeff.is_zero() {
            return Self::zero(vars);
        }

        let mut terms = BTreeMap::new();
        terms.insert(exponents, coeff);
        Self {
            vars: vars.to_vec(),
            terms,
        }
    }

    pub fn from_terms(
        vars: &[String],
        terms: impl IntoIterator<Item = (Vec<usize>, Rational)>,
    ) -> Self {
        let mut result = Self::zero(vars);
        for (exp, coeff) in terms {
            assert_eq!(vars.len(), exp.len());
            result.add_term(exp, coeff);
        }
        result
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn is_nonzero_constant(&self) -> bool {
        self.terms.len() == 1
            && self
                .terms
                .keys()
                .next()
                .is_some_and(|exp| exp.iter().all(|power| *power == 0))
    }

    pub fn add(&self, rhs: &Self) -> Self {
        self.assert_same_vars(rhs);
        let mut result = self.clone();
        for (exp, coeff) in &rhs.terms {
            result.add_term(exp.clone(), coeff.clone());
        }
        result
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        self.assert_same_vars(rhs);
        let mut result = self.clone();
        for (exp, coeff) in &rhs.terms {
            result.add_term(exp.clone(), -coeff.clone());
        }
        result
    }

    pub fn neg(&self) -> Self {
        let terms = self
            .terms
            .iter()
            .map(|(exp, coeff)| (exp.clone(), -coeff.clone()))
            .collect();
        Self {
            vars: self.vars.clone(),
            terms,
        }
    }

    pub fn scale(&self, factor: Rational) -> Self {
        if factor.is_zero() {
            return Self::zero(&self.vars);
        }

        Self {
            vars: self.vars.clone(),
            terms: self
                .terms
                .iter()
                .map(|(exp, coeff)| (exp.clone(), coeff.clone() * factor.clone()))
                .collect(),
        }
    }

    pub fn monic(&self) -> Self {
        if let Some((_, coeff)) = self.leading_term() {
            self.scale(Rational::one() / coeff)
        } else {
            self.clone()
        }
    }

    pub fn mul_monomial(&self, coeff: Rational, exponents: &[usize]) -> Self {
        assert_eq!(self.vars.len(), exponents.len());
        if coeff.is_zero() || self.is_zero() {
            return Self::zero(&self.vars);
        }

        let terms = self.terms.iter().map(|(exp, term_coeff)| {
            let next_exp = exp
                .iter()
                .zip(exponents.iter())
                .map(|(left, right)| left + right)
                .collect::<Vec<_>>();
            (next_exp, term_coeff.clone() * coeff.clone())
        });
        Self::from_terms(&self.vars, terms)
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        self.assert_same_vars(rhs);
        let mut result = Self::zero(&self.vars);
        for (left_exp, left_coeff) in &self.terms {
            for (right_exp, right_coeff) in &rhs.terms {
                let exp = left_exp
                    .iter()
                    .zip(right_exp.iter())
                    .map(|(left, right)| left + right)
                    .collect();
                result.add_term(exp, left_coeff.clone() * right_coeff.clone());
            }
        }
        result
    }

    pub fn pow(&self, exponent: usize) -> Self {
        if exponent == 0 {
            return Self::constant(&self.vars, Rational::one());
        }

        let mut result = Self::constant(&self.vars, Rational::one());
        let mut base = self.clone();
        let mut exp = exponent;
        while exp > 0 {
            if exp % 2 == 1 {
                result = result.mul(&base);
            }
            exp /= 2;
            if exp > 0 {
                base = base.mul(&base);
            }
        }
        result
    }

    pub fn derivative(&self, var_index: usize) -> Self {
        let mut result = Self::zero(&self.vars);
        for (exp, coeff) in &self.terms {
            let power = exp[var_index];
            if power == 0 {
                continue;
            }
            let mut derived_exp = exp.clone();
            derived_exp[var_index] -= 1;
            result.add_term(
                derived_exp,
                coeff.clone() * Rational::from_i128(power as i128),
            );
        }
        result
    }

    pub fn partials(&self) -> Vec<Self> {
        (0..self.vars.len())
            .map(|index| self.derivative(index))
            .collect()
    }

    pub fn evaluate(&self, values: &[Rational]) -> Result<Rational, String> {
        if values.len() != self.vars.len() {
            return Err(format!(
                "expected {} values, received {}",
                self.vars.len(),
                values.len()
            ));
        }

        let mut total = Rational::zero();
        for (exp, coeff) in &self.terms {
            let mut term = coeff.clone();
            for (power, value) in exp.iter().zip(values.iter()) {
                term = term * value.pow(*power);
            }
            total = total + term;
        }
        Ok(total)
    }

    pub fn substitute(&self, replacements: &[Self]) -> Self {
        assert_eq!(self.vars.len(), replacements.len());
        let target_vars = if let Some(first) = replacements.first() {
            first.vars.clone()
        } else {
            Vec::new()
        };
        for replacement in replacements {
            assert_eq!(replacement.vars, target_vars);
        }

        let maximum_powers = (0..self.vars.len())
            .map(|variable_index| {
                self.terms
                    .keys()
                    .map(|exponent| exponent[variable_index])
                    .max()
                    .unwrap_or(0)
            })
            .collect::<Vec<_>>();
        let replacement_powers = replacements
            .iter()
            .zip(maximum_powers)
            .map(|(replacement, maximum_power)| {
                let mut powers = Vec::with_capacity(maximum_power + 1);
                powers.push(Self::constant(&target_vars, Rational::one()));
                for power in 1..=maximum_power {
                    powers.push(powers[power - 1].mul(replacement));
                }
                powers
            })
            .collect::<Vec<_>>();

        let mut result = Self::zero(&target_vars);
        for (exp, coeff) in &self.terms {
            let mut term = Self::constant(&target_vars, coeff.clone());
            for (var_index, power) in exp.iter().enumerate() {
                if *power > 0 {
                    term = term.mul(&replacement_powers[var_index][*power]);
                }
            }
            result = result.add(&term);
        }
        result
    }

    pub fn specialize(&self, assignments: &BTreeMap<usize, Rational>) -> Self {
        let mut result = Self::zero(&self.vars);
        for (exp, coeff) in &self.terms {
            let mut new_exp = exp.clone();
            let mut new_coeff = coeff.clone();
            for (index, value) in assignments {
                new_coeff = new_coeff * value.pow(exp[*index]);
                new_exp[*index] = 0;
            }
            result.add_term(new_exp, new_coeff);
        }
        result
    }

    pub fn translated_by(&self, assignments: &BTreeMap<usize, Rational>) -> Self {
        if assignments.values().all(Rational::is_zero) {
            return self.clone();
        }
        let replacements = (0..self.vars.len())
            .map(|index| {
                let variable = Self::var(&self.vars, index);
                if let Some(value) = assignments.get(&index) {
                    variable.add(&Self::constant(&self.vars, value.clone()))
                } else {
                    variable
                }
            })
            .collect::<Vec<_>>();
        self.substitute(&replacements)
    }

    pub fn translated_jet_by(
        &self,
        assignments: &BTreeMap<usize, Rational>,
        max_degree: usize,
    ) -> Self {
        if assignments.values().all(Rational::is_zero) {
            return Self::from_terms(
                &self.vars,
                self.terms.iter().filter_map(|(exponent, coefficient)| {
                    (exponent.iter().sum::<usize>() <= max_degree)
                        .then_some((exponent.clone(), coefficient.clone()))
                }),
            );
        }
        let integer_denominator = BigInt::from(1usize);
        if self
            .terms
            .values()
            .all(|coefficient| coefficient.den == integer_denominator)
            && assignments
                .values()
                .all(|value| value.den == integer_denominator)
        {
            return self.translated_integer_jet_by(assignments, max_degree);
        }

        let mut result = Self::zero(&self.vars);
        for (exponent, coefficient) in &self.terms {
            let mut partial =
                BTreeMap::from([(vec![0usize; self.vars.len()], coefficient.clone())]);
            for (index, power) in exponent.iter().copied().enumerate() {
                let Some(value) = assignments.get(&index) else {
                    if power > 0 {
                        partial = partial
                            .into_iter()
                            .filter_map(|(mut translated_exponent, translated_coefficient)| {
                                let degree = translated_exponent.iter().sum::<usize>() + power;
                                (degree <= max_degree).then(|| {
                                    translated_exponent[index] += power;
                                    (translated_exponent, translated_coefficient)
                                })
                            })
                            .collect();
                    }
                    continue;
                };

                let mut next = BTreeMap::<Vec<usize>, Rational>::new();
                for (translated_exponent, translated_coefficient) in partial {
                    let current_degree = translated_exponent.iter().sum::<usize>();
                    for translated_power in 0..=power.min(max_degree - current_degree) {
                        let factor = binomial_rational(power, translated_power)
                            * value.pow(power - translated_power);
                        if factor.is_zero() {
                            continue;
                        }
                        let mut next_exponent = translated_exponent.clone();
                        next_exponent[index] += translated_power;
                        let next_coefficient = translated_coefficient.clone() * factor;
                        let accumulated = next
                            .get(&next_exponent)
                            .cloned()
                            .unwrap_or_else(Rational::zero)
                            + next_coefficient;
                        if accumulated.is_zero() {
                            next.remove(&next_exponent);
                        } else {
                            next.insert(next_exponent, accumulated);
                        }
                    }
                }
                partial = next;
            }
            for (translated_exponent, translated_coefficient) in partial {
                result.add_term(translated_exponent, translated_coefficient);
            }
        }
        result
    }

    fn translated_integer_jet_by(
        &self,
        assignments: &BTreeMap<usize, Rational>,
        max_degree: usize,
    ) -> Self {
        let factors = (0..self.vars.len())
            .map(|index| {
                let maximum_power = self
                    .terms
                    .keys()
                    .map(|exponent| exponent[index])
                    .max()
                    .unwrap_or(0);
                let value = assignments
                    .get(&index)
                    .map(|value| value.num.clone())
                    .unwrap_or_else(|| BigInt::from(0usize));
                (0..=maximum_power)
                    .map(|power| {
                        (0..=power.min(max_degree))
                            .map(|translated_power| {
                                binomial_bigint(power, translated_power)
                                    * bigint_power(&value, power - translated_power)
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        if self.vars.len() == 3 {
            return self.translated_integer_trivariate_jet(&factors, max_degree);
        }
        let mut coefficients = BTreeMap::<Vec<usize>, BigInt>::new();
        for (exponent, coefficient) in &self.terms {
            accumulate_integer_jet_term(
                exponent,
                &factors,
                0,
                max_degree,
                &mut vec![0usize; self.vars.len()],
                coefficient.num.clone(),
                &mut coefficients,
            );
        }
        Self::from_terms(
            &self.vars,
            coefficients.into_iter().map(|(exponent, coefficient)| {
                (exponent, Rational::new(coefficient, BigInt::from(1usize)))
            }),
        )
    }

    fn translated_integer_trivariate_jet(
        &self,
        factors: &[Vec<Vec<BigInt>>],
        max_degree: usize,
    ) -> Self {
        let side = max_degree + 1;
        let mut coefficients = vec![BigInt::from(0usize); side * side * side];
        for (exponent, coefficient) in &self.terms {
            for first in 0..=exponent[0].min(max_degree) {
                let first_factor = &factors[0][exponent[0]][first];
                if first_factor == &BigInt::from(0usize) {
                    continue;
                }
                let first_coefficient = coefficient.num.clone() * first_factor;
                for second in 0..=exponent[1].min(max_degree - first) {
                    let second_factor = &factors[1][exponent[1]][second];
                    if second_factor == &BigInt::from(0usize) {
                        continue;
                    }
                    let second_coefficient = first_coefficient.clone() * second_factor;
                    for third in 0..=exponent[2].min(max_degree - first - second) {
                        let third_factor = &factors[2][exponent[2]][third];
                        if third_factor == &BigInt::from(0usize) {
                            continue;
                        }
                        let index = (first * side + second) * side + third;
                        coefficients[index] += second_coefficient.clone() * third_factor;
                    }
                }
            }
        }

        let mut terms = Vec::new();
        for first in 0..=max_degree {
            for second in 0..=max_degree - first {
                for third in 0..=max_degree - first - second {
                    let index = (first * side + second) * side + third;
                    if coefficients[index] != BigInt::from(0usize) {
                        terms.push((
                            vec![first, second, third],
                            Rational::new(coefficients[index].clone(), BigInt::from(1usize)),
                        ));
                    }
                }
            }
        }
        Self::from_terms(&self.vars, terms)
    }

    pub fn center_order(&self, center_indices: &[usize]) -> Option<usize> {
        self.terms
            .keys()
            .map(|exp| center_indices.iter().map(|index| exp[*index]).sum())
            .min()
    }

    pub fn homogeneous_part(&self, degree: usize) -> Self {
        let terms = self
            .terms
            .iter()
            .filter_map(|(exp, coeff)| {
                (exp.iter().sum::<usize>() == degree).then_some((exp.clone(), coeff.clone()))
            })
            .collect::<Vec<_>>();
        Self::from_terms(&self.vars, terms)
    }

    pub fn divide_by_var_power(&self, var_index: usize, power: usize) -> Result<Self, String> {
        if power == 0 {
            return Ok(self.clone());
        }

        let mut result = Self::zero(&self.vars);
        for (exp, coeff) in &self.terms {
            if exp[var_index] < power {
                return Err(format!(
                    "term {} is not divisible by {}^{}",
                    format_monomial(&self.vars, exp),
                    self.vars[var_index],
                    power
                ));
            }
            let mut new_exp = exp.clone();
            new_exp[var_index] -= power;
            result.add_term(new_exp, coeff.clone());
        }
        Ok(result)
    }

    pub fn leading_term(&self) -> Option<(Vec<usize>, Rational)> {
        self.terms
            .iter()
            .max_by(|(left_exp, _), (right_exp, _)| compare_monomials(left_exp, right_exp))
            .map(|(exp, coeff)| (exp.clone(), coeff.clone()))
    }

    pub fn total_degree(&self) -> usize {
        self.terms
            .keys()
            .map(|exp| exp.iter().sum::<usize>())
            .max()
            .unwrap_or(0)
    }

    pub fn homogeneous_degree(&self) -> Option<usize> {
        let mut degrees = self.terms.keys().map(|exp| exp.iter().sum::<usize>());
        let first = degrees.next()?;
        if degrees.all(|degree| degree == first) {
            Some(first)
        } else {
            None
        }
    }

    pub fn variable_index(&self, name: &str) -> Option<usize> {
        self.vars.iter().position(|var| var == name)
    }

    fn add_term(&mut self, exp: Vec<usize>, coeff: Rational) {
        if coeff.is_zero() {
            return;
        }

        let next = self.terms.get(&exp).cloned().unwrap_or_else(Rational::zero) + coeff;
        if next.is_zero() {
            self.terms.remove(&exp);
        } else {
            self.terms.insert(exp, next);
        }
    }

    fn assert_same_vars(&self, rhs: &Self) {
        assert_eq!(self.vars, rhs.vars, "polynomial variable mismatch");
    }
}

fn binomial_rational(n: usize, k: usize) -> Rational {
    Rational::new(binomial_bigint(n, k), BigInt::from(1usize))
}

fn binomial_bigint(n: usize, k: usize) -> BigInt {
    let k = k.min(n - k);
    let mut value = BigInt::from(1usize);
    for index in 1..=k {
        value = value * BigInt::from(n + 1 - index) / BigInt::from(index);
    }
    value
}

fn bigint_power(value: &BigInt, exponent: usize) -> BigInt {
    let mut result = BigInt::from(1usize);
    for _ in 0..exponent {
        result *= value;
    }
    result
}

fn accumulate_integer_jet_term(
    exponent: &[usize],
    factors: &[Vec<Vec<BigInt>>],
    variable_index: usize,
    remaining_degree: usize,
    translated_exponent: &mut [usize],
    coefficient: BigInt,
    coefficients: &mut BTreeMap<Vec<usize>, BigInt>,
) {
    if variable_index == exponent.len() {
        *coefficients
            .entry(translated_exponent.to_vec())
            .or_insert_with(|| BigInt::from(0usize)) += coefficient;
        return;
    }

    let power = exponent[variable_index];
    for translated_power in 0..=power.min(remaining_degree) {
        let factor = &factors[variable_index][power][translated_power];
        if factor == &BigInt::from(0usize) {
            continue;
        }
        translated_exponent[variable_index] = translated_power;
        accumulate_integer_jet_term(
            exponent,
            factors,
            variable_index + 1,
            remaining_degree - translated_power,
            translated_exponent,
            coefficient.clone() * factor,
            coefficients,
        );
    }
    translated_exponent[variable_index] = 0;
}

pub fn compare_monomials(left: &[usize], right: &[usize]) -> std::cmp::Ordering {
    let left_total = left.iter().sum::<usize>();
    let right_total = right.iter().sum::<usize>();
    left_total
        .cmp(&right_total)
        .then_with(|| left.iter().cmp(right.iter()))
}

impl fmt::Display for Poly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.terms.is_empty() {
            return write!(f, "0");
        }

        let mut rendered = String::new();
        for (position, (exp, coeff)) in self.terms.iter().rev().enumerate() {
            let negative = coeff.is_negative();
            let abs_coeff = coeff.abs();
            let monomial = format_monomial(&self.vars, exp);

            if position == 0 {
                if negative {
                    rendered.push('-');
                }
            } else if negative {
                rendered.push_str(" - ");
            } else {
                rendered.push_str(" + ");
            }

            if monomial == "1" {
                rendered.push_str(&abs_coeff.to_string());
            } else if abs_coeff.is_one() {
                rendered.push_str(&monomial);
            } else {
                rendered.push_str(&format!("{abs_coeff}*{monomial}"));
            }
        }

        write!(f, "{rendered}")
    }
}

fn format_monomial(vars: &[String], exp: &[usize]) -> String {
    let factors = vars
        .iter()
        .zip(exp.iter())
        .filter_map(|(var, power)| match *power {
            0 => None,
            1 => Some(var.clone()),
            _ => Some(format!("{var}^{power}")),
        })
        .collect::<Vec<_>>();

    if factors.is_empty() {
        "1".to_string()
    } else {
        factors.join("*")
    }
}

#[cfg(test)]
mod tests {
    use super::Poly;
    use crate::rational::Rational;
    use std::collections::BTreeMap;

    fn vars() -> Vec<String> {
        vec!["x".to_string(), "y".to_string()]
    }

    #[test]
    fn differentiates_polynomial() {
        let vars = vars();
        let x = Poly::var(&vars, 0);
        let y = Poly::var(&vars, 1);
        let poly = x.pow(2).mul(&y).add(&y.pow(3));

        assert_eq!(poly.derivative(0).to_string(), "2*x*y");
        assert_eq!(poly.derivative(1).to_string(), "x^2 + 3*y^2");
    }

    #[test]
    fn substitutes_polynomial() {
        let vars = vars();
        let x = Poly::var(&vars, 0);
        let y = Poly::var(&vars, 1);
        let poly = x.pow(2).add(&y);
        let replacements = vec![x.add(&y), y.clone()];

        assert_eq!(
            poly.substitute(&replacements).to_string(),
            "x^2 + 2*x*y + y^2 + y"
        );
    }

    #[test]
    fn evaluates_polynomial() {
        let vars = vars();
        let x = Poly::var(&vars, 0);
        let y = Poly::var(&vars, 1);
        let poly = x.pow(2).sub(&y);

        assert_eq!(
            poly.evaluate(&[Rational::from_i128(3), Rational::from_i128(9)])
                .unwrap(),
            Rational::zero()
        );
    }

    #[test]
    fn extracts_homogeneous_part() {
        let vars = vars();
        let x = Poly::var(&vars, 0);
        let y = Poly::var(&vars, 1);
        let poly = x.pow(2).add(&x.mul(&y)).add(&y.pow(3));

        assert_eq!(poly.homogeneous_part(2).to_string(), "x^2 + x*y");
    }

    #[test]
    fn integer_translated_jet_matches_full_translation() {
        let vars = vars();
        let x = Poly::var(&vars, 0);
        let y = Poly::var(&vars, 1);
        let poly = x.pow(5).mul(&y.pow(3)).add(&x.pow(2)).sub(&y);
        let assignments =
            BTreeMap::from([(0, Rational::from_i128(2)), (1, Rational::from_i128(-1))]);
        let full = poly.translated_by(&assignments);
        let expected = Poly::from_terms(
            &vars,
            full.terms
                .into_iter()
                .filter(|(exponent, _)| exponent.iter().sum::<usize>() <= 4),
        );

        assert_eq!(poly.translated_jet_by(&assignments, 4), expected);
    }
}
