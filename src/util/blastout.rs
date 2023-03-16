use std::num::{ParseFloatError, ParseIntError};
use std::path::Path;
use std::str::{ParseBoolError, FromStr};
use std::sync::Arc;

use itertools::Itertools;
use polars::lazy::prelude::{Expr, LazyCsvReader};
use polars::prelude::{col, LazyFrame, LiteralValue, NullValues, Schema, DataType};

use thiserror::Error;

use crate::util::filter::{FromStrFilter, Op};

pub mod fields {
    use std::{collections::HashMap, sync::LazyLock};

    use polars::datatypes::DataType;

    pub const QSEQID: (&str, DataType) = ("qseqid", DataType::Utf8);
    pub const QGI: (&str, DataType) = ("qgi", DataType::Utf8);
    pub const QACC: (&str, DataType) = ("qacc", DataType::Utf8);
    pub const QACCVER: (&str, DataType) = ("qaccver", DataType::Utf8);
    pub const QLEN: (&str, DataType) = ("qlen", DataType::Int32);
    pub const SSEQID: (&str, DataType) = ("sseqid", DataType::Utf8);
    pub const SALLSEQID: (&str, DataType) = ("sallseqid", DataType::Utf8);
    pub const SGI: (&str, DataType) = ("sgi", DataType::Utf8);
    pub const SALLGI: (&str, DataType) = ("sallgi", DataType::Utf8);
    pub const SACC: (&str, DataType) = ("sacc", DataType::Utf8);
    pub const SACCVER: (&str, DataType) = ("saccver", DataType::Utf8);
    pub const SALLACC: (&str, DataType) = ("sallacc", DataType::Utf8);
    pub const SLEN: (&str, DataType) = ("slen", DataType::Int32);
    pub const QSTART: (&str, DataType) = ("qstart", DataType::Int32);
    pub const QEND: (&str, DataType) = ("qend", DataType::Int32);
    pub const SSTART: (&str, DataType) = ("sstart", DataType::Int32);
    pub const SEND: (&str, DataType) = ("send", DataType::Int32);
    pub const QSEQ: (&str, DataType) = ("qseq", DataType::Utf8);
    pub const SSEQ: (&str, DataType) = ("sseq", DataType::Utf8);
    pub const EVALUE: (&str, DataType) = ("evalue", DataType::Float32);
    pub const BITSCORE: (&str, DataType) = ("bitscore", DataType::Float32);
    pub const SCORE: (&str, DataType) = ("score", DataType::Float32);
    pub const LENGTH: (&str, DataType) = ("length", DataType::Int32);
    pub const PIDENT: (&str, DataType) = ("pident", DataType::Float32);
    pub const NIDENT: (&str, DataType) = ("nident", DataType::Int32);
    pub const MISMATCH: (&str, DataType) = ("mismatch", DataType::Int32);
    pub const POSITIVE: (&str, DataType) = ("positive", DataType::Int32);
    pub const GAPOPEN: (&str, DataType) = ("gapopen", DataType::Int32);
    pub const GAPS: (&str, DataType) = ("gaps", DataType::Int32);
    pub const PPOS: (&str, DataType) = ("ppos", DataType::Float32);
    pub const FRAMES: (&str, DataType) = ("frames", DataType::Utf8);
    pub const QFRAME: (&str, DataType) = ("qframe", DataType::Utf8);
    pub const SFRAME: (&str, DataType) = ("sframe", DataType::Utf8);
    pub const BTOP: (&str, DataType) = ("btop", DataType::Utf8);
    pub const STAXID: (&str, DataType) = ("staxid", DataType::Int64);
    pub const SSCINAME: (&str, DataType) = ("ssciname", DataType::Utf8);
    pub const SCOMNAME: (&str, DataType) = ("scomname", DataType::Utf8);
    pub const SBLASTNAME: (&str, DataType) = ("sblastname", DataType::Utf8);
    pub const SSKINGDOM: (&str, DataType) = ("sskingdom", DataType::Utf8);
    pub const STAXIDS: (&str, DataType) = ("staxids", DataType::Utf8);
    pub const SSCINAMES: (&str, DataType) = ("sscinames", DataType::Utf8);
    pub const SCOMNAMES: (&str, DataType) = ("scomnames", DataType::Utf8);
    pub const SBLASTNAMES: (&str, DataType) = ("sblastnames", DataType::Utf8);
    pub const SSKINGDOMS: (&str, DataType) = ("sskingdoms", DataType::Utf8);
    pub const STITLE: (&str, DataType) = ("stitle", DataType::Utf8);
    pub const SALLTITLES: (&str, DataType) = ("salltitles", DataType::Utf8);
    pub const SSTRAND: (&str, DataType) = ("sstrand", DataType::Utf8);
    pub const QCOVS: (&str, DataType) = ("qcovs", DataType::Float32);
    pub const QCOVHSP: (&str, DataType) = ("qcovhsp", DataType::Float32);
    pub const QCOVUS: (&str, DataType) = ("qcovus", DataType::Float32);

    pub static FIELD_TYPES: LazyLock<HashMap<&'static str, DataType>> = LazyLock::new(|| {
        HashMap::from([
            QSEQID,
            QGI,
            QACC,
            QACCVER,
            QLEN,
            SSEQID,
            SALLSEQID,
            SGI,
            SALLGI,
            SACC,
            SACCVER,
            SALLACC,
            SLEN,
            QSTART,
            QEND,
            SSTART,
            SEND,
            QSEQ,
            SSEQ,
            EVALUE,
            BITSCORE,
            SCORE,
            LENGTH,
            PIDENT,
            NIDENT,
            MISMATCH,
            POSITIVE,
            GAPOPEN,
            GAPS,
            PPOS,
            FRAMES,
            QFRAME,
            SFRAME,
            BTOP,
            STAXID,
            SSCINAME,
            SCOMNAME,
            SBLASTNAME,
            SSKINGDOM,
            STAXIDS,
            SSCINAMES,
            SCOMNAMES,
            SBLASTNAMES,
            SSKINGDOMS,
            STITLE,
            SALLTITLES,
            SSTRAND,
            QCOVS,
            QCOVHSP,
            QCOVUS
        ])
    });

    pub const DEFAULT_BLASTN_COLUMNS: &[&str] = &[
        QACCVER.0, SACCVER.0, PIDENT.0, LENGTH.0, MISMATCH.0, GAPOPEN.0, QSTART.0, QEND.0,
        SSTART.0, SEND.0,
    ];
}

pub enum BlastOutFmt {
    // TSV without header
    Six {
        delimiter: u8,
        present_columns: Vec<String>,
    },

    // like Six, but with comments
    Seven {
        delimiter: u8,
        present_columns: Vec<String>,
    },
}

impl FromStr for BlastOutFmt {
    type Err = anyhow::Error;

    fn from_str(fmt: &str) -> Result<Self, Self::Err> {
        let split_fmt = fmt.split_whitespace().collect_vec();

        fn parse_blast_outfmt_6_or_7(mut fmt_args: &[&str]) -> anyhow::Result<(u8, Vec<String>)> {
            let delim = if let Some(delim_str) = fmt_args.first().and_then(|s| s.strip_prefix("delim=")) {
                if let Ok(delim_chr) = delim_str.chars().exactly_one() {
                    fmt_args = &fmt_args[1..];
                    delim_chr as u8
                } else {
                    anyhow::bail!("Error parsing delimiter specifier")
                }
            } else {
                b'\t'
            };

            let present_columns = if fmt_args.is_empty() {
                fields::DEFAULT_BLASTN_COLUMNS
                    .iter()
                    .map(|&s| s.to_owned())
                    .collect_vec()
            } else {
                fmt_args.iter().map(|&s| s.to_owned()).collect_vec()
            };

            Ok((delim, present_columns))
        }

        match split_fmt.first().and_then(|s| s.parse::<u32>().ok()) {
            Some(6) => {
                let (delimiter, present_columns) = parse_blast_outfmt_6_or_7(&split_fmt[1..])?;
                Ok(BlastOutFmt::Six {
                    delimiter,
                    present_columns,
                })
            }
            Some(7) => {
                let (delimiter, present_columns) = parse_blast_outfmt_6_or_7(&split_fmt[1..])?;
                Ok(BlastOutFmt::Seven {
                    delimiter,
                    present_columns,
                })
            }
            Some(n) => anyhow::bail!("Unsupported outfmt type: {n} (only 6 and 7 are supported)"),
            None => anyhow::bail!("Unable to parse outfmt type"),
        }
    }
}

#[derive(Clone)]
pub struct BlastOutFilter {
    column: String,
    expr: Expr,
}

impl BlastOutFilter {
    pub fn get_column(&self) -> &str {
        &self.column
    }

    pub fn get_expr(&self) -> &Expr {
        &self.expr
    }
}

#[derive(Debug, Error)]
pub enum BlastoutFilterParseError {
    #[error("Unknown column {:?}", .0)]
    UnknownColumn(String),

    #[error("Error parsing filter value")]
    BoolParseError(#[from] ParseBoolError),

    #[error("Error parsing filter value")]
    FloatParseError(#[from] ParseFloatError),

    #[error("Error parsing filter value")]
    IntParseError(#[from] ParseIntError),

    #[error("Unsupported filter data type for column {:?}: {:?}", .0, .1)]
    UnsupportedDataType(String, &'static DataType),
}

impl FromStrFilter for BlastOutFilter {
    type Err = BlastoutFilterParseError;

    fn try_from_parts(key: &str, op: Op, value: &str) -> Result<Self, BlastoutFilterParseError> {
        use BlastoutFilterParseError::*;

        let lhs = col(key);

        let dtype = fields::FIELD_TYPES
            .get(key)
            .ok_or_else(|| UnknownColumn(key.to_owned()))?;

        let rhs = match dtype {
            DataType::Boolean => LiteralValue::Boolean(value.parse()?),
            DataType::Float32 => LiteralValue::Float32(value.parse()?),
            DataType::Float64 => LiteralValue::Float64(value.parse()?),
            DataType::Int32 => LiteralValue::Int32(value.parse()?),
            DataType::Int64 => LiteralValue::Int32(value.parse()?),
            DataType::UInt32 => LiteralValue::Int32(value.parse()?),
            DataType::UInt64 => LiteralValue::Int32(value.parse()?),
            DataType::Utf8 => LiteralValue::Utf8(value.to_owned()),
            _ => return Err(UnsupportedDataType(key.to_owned(), dtype)),
        };

        let rhs = Expr::Literal(rhs);

        let expr = match op {
            Op::Eq => lhs.eq(rhs),
            Op::Neq => lhs.neq(rhs),
            Op::Lt => lhs.lt(rhs),
            Op::Leq => lhs.lt_eq(rhs),
            Op::Gt => lhs.gt(rhs),
            Op::Geq => lhs.gt_eq(rhs),
        };

        Ok(BlastOutFilter {
            column: key.to_owned(),
            expr,
        })
    }
}

pub fn apply_filters(df: LazyFrame, filters: Vec<BlastOutFilter>) -> LazyFrame {
    if let Some(filter) = filters.into_iter().map(|filter| filter.expr).reduce(Expr::and) {
        df.filter(filter)
    } else {
        df
    }
}

pub fn load_blastout(
    path: impl AsRef<Path>,
    blast_outfmt: &BlastOutFmt,
    wanted_columns: Option<Vec<&str>>,
) -> anyhow::Result<LazyFrame> {
    let (delim, comment_char, present_columns) = match blast_outfmt {
        BlastOutFmt::Six {
            delimiter,
            present_columns,
        } => (*delimiter, None, present_columns),
        BlastOutFmt::Seven {
            delimiter,
            present_columns,
        } => (*delimiter, Some(b'#'), present_columns),
    };

    let schema = {
        let mut schema = Schema::new();

        let mut n_ignored = 0;

        for present in present_columns {
            if present.starts_with('_') {
                schema.with_column(format!("_{n_ignored}"), DataType::Utf8);
                n_ignored += 1;
            } else {
                let dtype = if let Some(dtype) = fields::FIELD_TYPES.get(present.as_str()) {
                    dtype
                } else {
                    anyhow::bail!("Unrecognized BLAST+ column name {present}")
                };

                schema.with_column(present.clone(), dtype.clone());
            }
        }

        schema
    };

    let result = LazyCsvReader::new(path)
        .has_header(false)
        .with_delimiter(delim)
        .with_comment_char(comment_char)
        .with_schema(Arc::new(schema))
        .with_null_values(Some(NullValues::AllColumnsSingle("N/A".to_owned())))
        .finish()?;

    let wanted_columns: Vec<&str> = wanted_columns.unwrap_or_else(|| {
        present_columns
            .iter()
            .filter(|p| !p.starts_with('_'))
            .map(|s| s.as_ref())
            .collect()
    });

    let result = result.select(wanted_columns.iter().map(|name| col(name)).collect_vec());

    Ok(result)
}
