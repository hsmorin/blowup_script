use std::cmp::Ordering;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct Rational {
    pub num: BigInt,
    pub den: BigInt,
}

impl Rational {
    pub fn new(num: impl Into<BigInt>, den: impl Into<BigInt>) -> Self {
        let mut num = num.into();
        let mut den = den.into();

        assert!(!den.is_zero(), "zero denominator");
        if num.is_zero() {
            return Self {
                num: BigInt::zero(),
                den: BigInt::one(),
            };
        }

        if den.is_negative() {
            num = -num;
            den = -den;
        }

        let gcd = num.abs().gcd(&den);
        Self {
            num: num / &gcd,
            den: den / gcd,
        }
    }

    pub fn zero() -> Self {
        Self {
            num: BigInt::zero(),
            den: BigInt::one(),
        }
    }

    pub fn one() -> Self {
        Self {
            num: BigInt::one(),
            den: BigInt::one(),
        }
    }

    pub fn from_i128(value: i128) -> Self {
        Self {
            num: BigInt::from(value),
            den: BigInt::one(),
        }
    }

    pub fn is_zero(&self) -> bool {
        self.num.is_zero()
    }

    pub fn is_one(&self) -> bool {
        self.num.is_one() && self.den.is_one()
    }

    pub fn is_negative(&self) -> bool {
        self.num.is_negative()
    }

    pub fn abs(&self) -> Self {
        Self::new(self.num.abs(), self.den.clone())
    }

    pub fn pow(&self, exp: usize) -> Self {
        let mut result = Self::one();
        for _ in 0..exp {
            result = result * self.clone();
        }
        result
    }

    pub fn to_nonnegative_usize(&self) -> Option<usize> {
        (self.den.is_one() && !self.num.is_negative())
            .then(|| self.num.to_usize())
            .flatten()
    }

    pub fn numerator_i128(&self) -> Option<i128> {
        self.num.to_i128()
    }

    pub fn denominator_i128(&self) -> Option<i128> {
        self.den.to_i128()
    }

    pub fn parse(input: &str) -> Result<Self, String> {
        let trimmed = input.trim();
        if trimmed.is_empty() {
            return Err("empty rational".to_string());
        }

        if let Some((left, right)) = trimmed.split_once('/') {
            let num = parse_bigint(left)?;
            let den = parse_bigint(right)?;
            if den.is_zero() {
                return Err("zero denominator".to_string());
            }
            Ok(Self::new(num, den))
        } else {
            Ok(Self::new(parse_bigint(trimmed)?, BigInt::one()))
        }
    }
}

impl Add for Rational {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(
            &self.num * &rhs.den + &rhs.num * &self.den,
            &self.den * &rhs.den,
        )
    }
}

impl Sub for Rational {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(
            &self.num * &rhs.den - &rhs.num * &self.den,
            &self.den * &rhs.den,
        )
    }
}

impl Mul for Rational {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(&self.num * &rhs.num, &self.den * &rhs.den)
    }
}

impl Div for Rational {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        assert!(!rhs.is_zero(), "division by zero rational");
        Self::new(&self.num * &rhs.den, &self.den * &rhs.num)
    }
}

impl Neg for Rational {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            num: -self.num,
            den: self.den,
        }
    }
}

impl Ord for Rational {
    fn cmp(&self, other: &Self) -> Ordering {
        (&self.num * &other.den).cmp(&(&other.num * &self.den))
    }
}

impl PartialOrd for Rational {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.den.is_one() {
            write!(f, "{}", self.num)
        } else {
            write!(f, "{}/{}", self.num, self.den)
        }
    }
}

fn parse_bigint(input: &str) -> Result<BigInt, String> {
    input
        .trim()
        .parse::<BigInt>()
        .map_err(|_| format!("invalid integer '{input}'"))
}

#[cfg(test)]
mod tests {
    use super::Rational;

    #[test]
    fn normalizes_sign_and_gcd() {
        assert_eq!(Rational::new(2, -4).to_string(), "-1/2");
        assert_eq!(Rational::new(0, -4).to_string(), "0");
    }

    #[test]
    fn parses_integer_and_fraction() {
        assert_eq!(Rational::parse("-3").unwrap(), Rational::new(-3, 1));
        assert_eq!(Rational::parse("6/9").unwrap(), Rational::new(2, 3));
    }

    #[test]
    fn adds_with_common_denominator_before_multiplying() {
        let value = Rational::new(i128::MAX / 2, 3);

        assert_eq!(value.clone() + value, Rational::new(i128::MAX - 1, 3));
    }

    #[test]
    fn adds_beyond_i128_range() {
        let value = Rational::from_i128(i128::MAX);

        assert_eq!(
            (value.clone() + value).to_string(),
            "340282366920938463463374607431768211454"
        );
    }

    #[test]
    fn multiplies_after_cross_cancellation() {
        let numerator = 1i128 << 100;
        let denominator = (1i128 << 99) - 1;

        assert_eq!(
            Rational::new(numerator, denominator) * Rational::new(denominator, numerator),
            Rational::one()
        );
    }
}
