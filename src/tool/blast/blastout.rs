use std::num::{ParseFloatError, ParseIntError};
use std::path::Path;
use std::str::{ParseBoolError, FromStr};
use std::sync::Arc;

use csv::StringRecord;
use itertools::Itertools;
use polars::lazy::prelude::{Expr, LazyCsvReader};
use polars::prelude::{col, DataType, LazyFileListReader, LazyFrame, LiteralValue, NullValues, Schema};

use thiserror::Error;

use crate::util::filter::{FromStrFilter, Op};

pub mod fields {
    use std::{collections::HashMap, sync::LazyLock};

    use polars::datatypes::DataType;

    pub const QSEQID: (&str, DataType) = ("qseqid", DataType::String);
    pub const QGI: (&str, DataType) = ("qgi", DataType::String);
    pub const QACC: (&str, DataType) = ("qacc", DataType::String);
    pub const QACCVER: (&str, DataType) = ("qaccver", DataType::String);
    pub const QLEN: (&str, DataType) = ("qlen", DataType::Int32);
    pub const SSEQID: (&str, DataType) = ("sseqid", DataType::String);
    pub const SALLSEQID: (&str, DataType) = ("sallseqid", DataType::String);
    pub const SGI: (&str, DataType) = ("sgi", DataType::String);
    pub const SALLGI: (&str, DataType) = ("sallgi", DataType::String);
    pub const SACC: (&str, DataType) = ("sacc", DataType::String);
    pub const SACCVER: (&str, DataType) = ("saccver", DataType::String);
    pub const SALLACC: (&str, DataType) = ("sallacc", DataType::String);
    pub const SLEN: (&str, DataType) = ("slen", DataType::Int32);
    pub const QSTART: (&str, DataType) = ("qstart", DataType::Int32);
    pub const QEND: (&str, DataType) = ("qend", DataType::Int32);
    pub const SSTART: (&str, DataType) = ("sstart", DataType::Int32);
    pub const SEND: (&str, DataType) = ("send", DataType::Int32);
    pub const QSEQ: (&str, DataType) = ("qseq", DataType::String);
    pub const SSEQ: (&str, DataType) = ("sseq", DataType::String);
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
    pub const FRAMES: (&str, DataType) = ("frames", DataType::String);
    pub const QFRAME: (&str, DataType) = ("qframe", DataType::String);
    pub const SFRAME: (&str, DataType) = ("sframe", DataType::String);
    pub const BTOP: (&str, DataType) = ("btop", DataType::String);
    pub const STAXID: (&str, DataType) = ("staxid", DataType::Int64);
    pub const SSCINAME: (&str, DataType) = ("ssciname", DataType::String);
    pub const SCOMNAME: (&str, DataType) = ("scomname", DataType::String);
    pub const SBLASTNAME: (&str, DataType) = ("sblastname", DataType::String);
    pub const SSKINGDOM: (&str, DataType) = ("sskingdom", DataType::String);
    pub const STAXIDS: (&str, DataType) = ("staxids", DataType::String);
    pub const SSCINAMES: (&str, DataType) = ("sscinames", DataType::String);
    pub const SCOMNAMES: (&str, DataType) = ("scomnames", DataType::String);
    pub const SBLASTNAMES: (&str, DataType) = ("sblastnames", DataType::String);
    pub const SSKINGDOMS: (&str, DataType) = ("sskingdoms", DataType::String);
    pub const STITLE: (&str, DataType) = ("stitle", DataType::String);
    pub const SALLTITLES: (&str, DataType) = ("salltitles", DataType::String);
    pub const SSTRAND: (&str, DataType) = ("sstrand", DataType::String);
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
        delimiter: char,
        present_columns: Vec<String>,
    },

    // like Six, but with comments
    Seven {
        delimiter: char,
        present_columns: Vec<String>,
    },
}

impl BlastOutFmt {
    pub fn comment_char(&self) -> Option<char> {
        match &self {
            Self::Six { .. } => None,
            Self::Seven { .. } => Some('#'),
        }
    }

    pub fn comment_prefix(&self) -> Option<String> {
        Some(self.comment_char()?.to_string())
    }

    pub fn delimiter(&self) -> char {
        match &self {
            Self::Six { delimiter, .. } => *delimiter,
            Self::Seven { delimiter, .. } => *delimiter,
        }
    }

    pub fn present_columns(&self) -> &[String] {
        match &self {
            Self::Six { present_columns, .. } => present_columns,
            Self::Seven { present_columns, .. } => present_columns,
        }
    }

    pub fn column_index(&self, column: &str) -> Option<usize> {
        self.present_columns().iter().position(|c| c == column)
    }

    pub fn csv_reader_builder(&self) -> csv::ReaderBuilder {
        let mut builder = csv::ReaderBuilder::new();

        builder
            .delimiter(self.delimiter() as u8)
            .has_headers(false)
            .comment(self.comment_char().map(|c| c as u8));

        builder
    }
}

impl FromStr for BlastOutFmt {
    type Err = anyhow::Error;

    fn from_str(fmt: &str) -> Result<Self, Self::Err> {
        let split_fmt = fmt.split_whitespace().collect_vec();

        fn parse_blast_outfmt_6_or_7(mut fmt_args: &[&str]) -> anyhow::Result<(char, Vec<String>)> {
            let delim = if let Some(delim_str) = fmt_args.first().and_then(|s| s.strip_prefix("delim=")) {
                if let Ok(delim_chr) = delim_str.chars().exactly_one() {
                    fmt_args = &fmt_args[1..];
                    delim_chr
                } else {
                    anyhow::bail!("Error parsing delimiter specifier")
                }
            } else {
                '\t'
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
pub enum Literal {
    String(String),
    Int64(i64),
    UInt64(u64),
    Float64(f64),
}

impl From<Literal> for LiteralValue {
    fn from(value: Literal) -> Self {
        match value {
            Literal::String(s) => LiteralValue::String(s),
            Literal::Int64(x) => LiteralValue::Int64(x),
            Literal::UInt64(x) => LiteralValue::UInt64(x),
            Literal::Float64(x) => LiteralValue::Float64(x),
        }
    }
}

#[derive(Clone)]
pub struct BlastOutFilter {
    column: String,
    op: Op,
    value: Literal,
}

type StringRecordPredicate = Box<dyn for<'a> Fn(&'a StringRecord) -> Option<bool>>;

impl BlastOutFilter {
    pub fn get_column(&self) -> &str {
        &self.column
    }

    pub fn into_polars_expr(self) -> Expr {
        let lhs = col(&self.column);

        let rhs = Expr::Literal(self.value.into());

        match self.op {
            Op::Eq => lhs.eq(rhs),
            Op::Neq => lhs.neq(rhs),
            Op::Lt => lhs.lt(rhs),
            Op::Leq => lhs.lt_eq(rhs),
            Op::Gt => lhs.gt(rhs),
            Op::Geq => lhs.gt_eq(rhs),
        }
    }

    pub fn into_string_record_predicate<S: AsRef<str>>(
        self,
        schema: &[S],
    ) -> Option<StringRecordPredicate> {
        let i = schema.iter().position(|col| col.as_ref() == self.column)?;

        let op = self.op;

        let f: StringRecordPredicate = match self.value {
            Literal::String(s) => Box::new(move |sr| Some(op.apply(&sr[i], &s))),
            Literal::Int64(x) => Box::new(move |sr| Some(op.apply(&sr[i].parse().ok()?, &x))),
            Literal::UInt64(x) => Box::new(move |sr| Some(op.apply(&sr[i].parse().ok()?, &x))),
            Literal::Float64(x) => Box::new(move |sr| Some(op.apply(&sr[i].parse().ok()?, &x))),
        };

        Some(f)
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

        let dtype = fields::FIELD_TYPES
            .get(key)
            .ok_or_else(|| UnknownColumn(key.to_owned()))?;

        let value = match dtype {
            DataType::Float32 => Literal::Float64(value.parse()?),
            DataType::Float64 => Literal::Float64(value.parse()?),
            DataType::Int32 => Literal::Int64(value.parse()?),
            DataType::Int64 => Literal::Int64(value.parse()?),
            DataType::UInt32 => Literal::Int64(value.parse()?),
            DataType::UInt64 => Literal::Int64(value.parse()?),
            DataType::String => Literal::String(value.to_owned()),
            _ => return Err(UnsupportedDataType(key.to_owned(), dtype)),
        };

        Ok(BlastOutFilter {
            column: key.to_owned(),
            op,
            value,
        })
    }
}

pub fn apply_filters(df: LazyFrame, filters: Vec<BlastOutFilter>) -> LazyFrame {
    if let Some(filter) = filters.into_iter().map(|filter| filter.into_polars_expr()).reduce(Expr::and) {
        df.filter(filter)
    } else {
        df
    }
}

pub fn filters_into_string_record_predicate(
    filters: Vec<BlastOutFilter>,
    outfmt: &BlastOutFmt,
) -> Option<impl for<'a> Fn(&'a StringRecord) -> Option<bool>> {
    fn and(preds: Vec<StringRecordPredicate>) -> impl for<'a> Fn(&'a StringRecord) -> Option<bool> {
        move |sr| {
            for pred in preds.iter() {
                if !pred(sr)? {
                    return Some(false);
                }
            }
            Some(true)
        }
    }

    let preds = filters
        .into_iter()
        .map(|f| f.into_string_record_predicate(outfmt.present_columns()))
        .collect::<Option<_>>()?;

    Some(and(preds))
}

pub fn load_blastout(
    path: impl AsRef<Path>,
    blast_outfmt: &BlastOutFmt,
    wanted_columns: Option<Vec<&str>>,
) -> anyhow::Result<LazyFrame> {
    let delim = blast_outfmt.delimiter();
    let comment_prefix = blast_outfmt.comment_prefix();
    let present_columns = blast_outfmt.present_columns();

    let schema = {
        let mut schema = Schema::new();

        let mut n_ignored = 0;

        for present in present_columns {
            if present.starts_with('_') {
                schema.with_column(format!("_{n_ignored}").into(), DataType::String);
                n_ignored += 1;
            } else {
                let dtype = if let Some(dtype) = fields::FIELD_TYPES.get(present.as_str()) {
                    dtype
                } else {
                    anyhow::bail!("Unrecognized BLAST+ column name {present}")
                };

                schema.with_column(present.into(), dtype.clone());
            }
        }

        schema
    };

    let result = LazyCsvReader::new(path)
        .with_has_header(false)
        .with_separator(delim as u8)
        .with_comment_prefix(comment_prefix.as_deref())
        .with_schema(Some(Arc::new(schema)))
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
