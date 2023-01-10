use regex::Regex;
use std::error::Error;
use std::fmt::Debug;
use std::sync::LazyLock;
use std::str::FromStr;
use thiserror::Error;

#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
pub enum Op {
    Eq,
    Neq,
    Lt,
    Leq,
    Gt,
    Geq,
}

impl Op {
    pub fn apply<T: PartialOrd>(self, a: &T, b: &T) -> bool {
        match self {
            Op::Eq => a == b,
            Op::Neq => a != b,
            Op::Lt => a < b,
            Op::Leq => a <= b,
            Op::Gt => a > b,
            Op::Geq => a >= b,
        }
    }
}

impl FromStr for Op {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use Op::*;

        match s {
            "<" => Ok(Lt),
            "<=" => Ok(Leq),
            ">" => Ok(Gt),
            ">=" => Ok(Geq),
            "=" | "==" => Ok(Eq),
            "!=" => Ok(Neq),
            _ => Err(()),
        }
    }
}

#[derive(Error)]
pub enum FilterParseError<F: FromStrFilter> {
    #[error("Could not parse filter {}", .0)]
    ParseError(String),

    #[error("Error building filter from parts of {}: {}", .0, .1)]
    BuildError(String, #[source] F::Err),
}

impl<F> Debug for FilterParseError<F>
where
    F: FromStrFilter,
    F::Err: Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ParseError(arg0) => f.debug_tuple("ParseError").field(arg0).finish(),
            Self::BuildError(arg0, arg1) => {
                f.debug_tuple("BuildError").field(arg0).field(arg1).finish()
            }
        }
    }
}

pub trait FromStrFilter: Sized {
    type Err: Error + 'static;
    fn try_from_parts(key: &str, op: Op, value: &str) -> Result<Self, Self::Err>;

    fn parse_filter(filter_str: &str) -> Result<Self, FilterParseError<Self>> {
        use FilterParseError::*;

        static RE: LazyLock<Regex> =
            LazyLock::new(|| Regex::new(r"^(\w+)\s*([><]?=?)\s*(\S+)$").unwrap());

        let cap = RE
            .captures(filter_str)
            .ok_or_else(|| ParseError(filter_str.to_owned()))?;

        let key = &cap[1];

        let op = cap[2]
            .parse()
            .map_err(|()| ParseError(filter_str.to_owned()))?;

        let value_str = &cap[3];

        Self::try_from_parts(key, op, value_str).map_err(|e| BuildError(filter_str.to_owned(), e))
    }

    fn parse_filters<I, S>(filter_strs: I) -> Result<Vec<Self>, FilterParseError<Self>>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        filter_strs
            .into_iter()
            .map(|f| Self::parse_filter(f.as_ref()))
            .collect()
    }
}
