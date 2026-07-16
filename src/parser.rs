use crate::poly::Poly;
use crate::rational::Rational;

#[derive(Clone, Debug, PartialEq, Eq)]
enum Token {
    Number(Rational),
    Ident(String),
    Plus,
    Minus,
    Star,
    Caret,
    LParen,
    RParen,
    End,
}

pub fn parse_polynomial(input: &str, vars: &[String]) -> Result<Poly, String> {
    let tokens = tokenize(input, vars)?;
    let mut parser = Parser {
        tokens,
        position: 0,
        vars: vars.to_vec(),
    };
    let poly = parser.parse_expr()?;
    if parser.peek() != &Token::End {
        return Err(format!("unexpected token {:?}", parser.peek()));
    }
    Ok(poly)
}

struct Parser {
    tokens: Vec<Token>,
    position: usize,
    vars: Vec<String>,
}

impl Parser {
    fn parse_expr(&mut self) -> Result<Poly, String> {
        self.parse_sum()
    }

    fn parse_sum(&mut self) -> Result<Poly, String> {
        let mut result = self.parse_product()?;
        loop {
            match self.peek() {
                Token::Plus => {
                    self.advance();
                    result = result.add(&self.parse_product()?);
                }
                Token::Minus => {
                    self.advance();
                    result = result.sub(&self.parse_product()?);
                }
                _ => break,
            }
        }
        Ok(result)
    }

    fn parse_product(&mut self) -> Result<Poly, String> {
        let mut result = self.parse_power()?;
        loop {
            match self.peek() {
                Token::Star => {
                    self.advance();
                    result = result.mul(&self.parse_power()?);
                }
                token if starts_implicit_factor(token) => {
                    result = result.mul(&self.parse_power()?);
                }
                _ => break,
            }
        }
        Ok(result)
    }

    fn parse_power(&mut self) -> Result<Poly, String> {
        let base = self.parse_unary()?;
        if self.peek() != &Token::Caret {
            return Ok(base);
        }

        self.advance();
        let exponent = match self.advance() {
            Token::Number(value) => value
                .to_nonnegative_usize()
                .ok_or_else(|| format!("expected nonnegative integer exponent, got {value:?}"))?,
            token => {
                return Err(format!(
                    "expected nonnegative integer exponent, got {token:?}"
                ));
            }
        };
        Ok(base.pow(exponent))
    }

    fn parse_unary(&mut self) -> Result<Poly, String> {
        match self.peek() {
            Token::Plus => {
                self.advance();
                self.parse_unary()
            }
            Token::Minus => {
                self.advance();
                Ok(self.parse_unary()?.neg())
            }
            _ => self.parse_factor(),
        }
    }

    fn parse_factor(&mut self) -> Result<Poly, String> {
        match self.advance() {
            Token::Number(value) => Ok(Poly::constant(&self.vars, value)),
            Token::Ident(name) => {
                let index = self
                    .vars
                    .iter()
                    .position(|var| var == &name)
                    .ok_or_else(|| format!("unknown variable '{name}'"))?;
                Ok(Poly::var(&self.vars, index))
            }
            Token::LParen => {
                let result = self.parse_expr()?;
                match self.advance() {
                    Token::RParen => Ok(result),
                    token => Err(format!("expected ')', got {token:?}")),
                }
            }
            token => Err(format!("expected number, variable, or '(', got {token:?}")),
        }
    }

    fn peek(&self) -> &Token {
        &self.tokens[self.position]
    }

    fn advance(&mut self) -> Token {
        let token = self.tokens[self.position].clone();
        self.position += 1;
        token
    }
}

fn starts_implicit_factor(token: &Token) -> bool {
    matches!(token, Token::Number(_) | Token::Ident(_) | Token::LParen)
}

fn tokenize(input: &str, vars: &[String]) -> Result<Vec<Token>, String> {
    let chars = input.chars().collect::<Vec<_>>();
    let mut tokens = Vec::new();
    let mut position = 0;
    let mut known_vars = vars.to_vec();
    known_vars.sort_by_key(|var| std::cmp::Reverse(var.len()));

    while position < chars.len() {
        let ch = chars[position];
        match ch {
            ' ' | '\t' | '\n' | '\r' => position += 1,
            '+' => {
                tokens.push(Token::Plus);
                position += 1;
            }
            '-' => {
                tokens.push(Token::Minus);
                position += 1;
            }
            '*' => {
                tokens.push(Token::Star);
                position += 1;
            }
            '^' => {
                tokens.push(Token::Caret);
                position += 1;
            }
            '(' => {
                tokens.push(Token::LParen);
                position += 1;
            }
            ')' => {
                tokens.push(Token::RParen);
                position += 1;
            }
            '0'..='9' => {
                let start = position;
                while position < chars.len() && chars[position].is_ascii_digit() {
                    position += 1;
                }
                if position < chars.len() && chars[position] == '/' {
                    position += 1;
                    let denom_start = position;
                    while position < chars.len() && chars[position].is_ascii_digit() {
                        position += 1;
                    }
                    if denom_start == position {
                        return Err("expected denominator after '/'".to_string());
                    }
                }
                let text = chars[start..position].iter().collect::<String>();
                tokens.push(Token::Number(Rational::parse(&text)?));
            }
            _ if is_ident_start(ch) => {
                let remaining = chars[position..].iter().collect::<String>();
                if let Some(var) = known_vars.iter().find(|var| remaining.starts_with(*var)) {
                    tokens.push(Token::Ident(var.clone()));
                    position += var.chars().count();
                } else {
                    let start = position;
                    position += 1;
                    while position < chars.len() && is_ident_continue(chars[position]) {
                        position += 1;
                    }
                    let ident = chars[start..position].iter().collect::<String>();
                    tokens.push(Token::Ident(ident));
                }
            }
            _ => return Err(format!("unexpected character '{ch}'")),
        }
    }

    tokens.push(Token::End);
    Ok(tokens)
}

fn is_ident_start(ch: char) -> bool {
    ch.is_ascii_alphabetic() || ch == '_'
}

fn is_ident_continue(ch: char) -> bool {
    ch.is_ascii_alphanumeric() || ch == '_'
}

#[cfg(test)]
mod tests {
    use super::parse_polynomial;

    #[test]
    fn parses_implicit_multiplication() {
        let vars = vec!["x0".to_string(), "x1".to_string()];
        let poly = parse_polynomial("2x0x1 + (x0+x1)^2", &vars).unwrap();
        assert_eq!(poly.to_string(), "x0^2 + 4*x0*x1 + x1^2");
    }

    #[test]
    fn rejects_unknown_variable() {
        let vars = vec!["x0".to_string()];
        assert!(parse_polynomial("x1", &vars).is_err());
    }
}
