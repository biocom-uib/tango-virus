use std::num::{ParseFloatError, ParseIntError};
use std::str::ParseBoolError;
use std::sync::Arc;

use anyhow::Context;
use clap::{ArgEnum, Args};
use itertools::Itertools;

use polars::datatypes::DataType;
use polars::lazy::prelude::LazyCsvReader;
use polars::prelude::{
    col, CsvWriter, DistinctKeepStrategy, Expr, GetOutput, IntoSeries, LazyFrame, LiteralValue,
    NullValues, Schema,
};
use polars_io::SerWriter;
use thiserror::Error;

use crate::filter::{FromStrFilter, Op};
use crate::util::writing_new_file_or_stdout;

pub mod fields {
    use std::collections::HashMap;

    use lazy_static::lazy_static;
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

    lazy_static! {
        pub static ref FIELD_TYPES: HashMap<&'static str, DataType> = HashMap::from([
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
        ]);
    }

    pub const DEFAULT_BLASTN_COLUMNS: &[&str] = &[
        QACCVER.0, SACCVER.0, PIDENT.0, LENGTH.0, MISMATCH.0, GAPOPEN.0, QSTART.0, QEND.0,
        SSTART.0, SEND.0,
    ];
}

#[derive(ArgEnum, Debug, Copy, Clone)]
pub enum WeightColAgg {
    Max,
    Mean,
    Median,
    Min,
    NumUnique,
    Product,
    Sum,
}

impl WeightColAgg {
    fn into_expr_fn(self) -> fn(Expr) -> Expr {
        use WeightColAgg::*;

        match self {
            Max => Expr::max,
            Mean => Expr::mean,
            Median => Expr::median,
            Min => Expr::min,
            NumUnique => Expr::n_unique,
            Product => Expr::product,
            Sum => Expr::sum,
        }
    }
}

#[derive(Args)]
pub struct PreprocessBlastOutArgs {
    /// Path to the BLAST+ output file. STDIN is not supported.
    blastout_path: String,

    /// Preprocessed output file. Use '-' to write to stdout.
    #[clap(short, long, default_value = "-")]
    output: String,

    /// BLAST+ output format. Only 6 and 7 are supported at the moment. Columns can be specified
    /// just like in BLAST+, like '7 qaccver saccver pident bitscore'. Irrelevant columns can be
    /// ignored using _ instead, e.g. '7 qaccver saccver _ bitscore'.
    #[clap(long)]
    blast_outfmt: String,

    /// Column name (as specified in outfmt) to use as query identifier.
    #[clap(long)]
    query_id_col: String,

    /// Column name (as specified in outfmt) to use as subject identifier.
    #[clap(long)]
    subject_id_col: String,

    /// Column name (as specified in outfmt) to use as weight.
    #[clap(long)]
    weight_col: Option<String>,

    /// Aggregation function to apply to the selected weight column (if needed). [requires
    /// --weight-col]
    #[clap(long, arg_enum, requires = "weight-col")]
    weight_col_agg: Option<WeightColAgg>,

    /// Filters to apply before processing. Example: --filter 'evalue<=1e-3'
    #[clap(long)]
    filter: Vec<String>,
}

enum BlastOutFmt {
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

fn parse_blast_outfmt(fmt: &str) -> anyhow::Result<BlastOutFmt> {
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

struct BlastOutFilter(Expr);

#[derive(Debug, Error)]
enum BlastoutFilterParseError {
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

        Ok(BlastOutFilter(expr))
    }
}

fn apply_filters(df: LazyFrame, filters: impl IntoIterator<Item = BlastOutFilter>) -> LazyFrame {
    if let Some(filter) = filters.into_iter().map(|filter| filter.0).reduce(Expr::and) {
        df.filter(filter)
    } else {
        df
    }
}

fn load_blastout(
    path: String,
    blast_outfmt: &BlastOutFmt,
    wanted_columns: Option<Vec<String>>,
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
        .with_null_values(Some(NullValues::AllColumns("N/A".to_owned())))
        .finish()?;

    let wanted_columns = wanted_columns.unwrap_or_else(|| {
        present_columns
            .iter()
            .filter(|p| !p.starts_with('_'))
            .cloned()
            .collect()
    });

    let result = result.select(wanted_columns.iter().map(|name| col(name)).collect_vec());

    Ok(result)
}

fn group_blast_hits(hits: LazyFrame, query_id_col: &str, subject_id_col: &str) -> LazyFrame {
    hits.filter(
        col(query_id_col)
            .is_not_null()
            .and(col(subject_id_col).is_not_null()),
    )
    .select([col(query_id_col), col(subject_id_col).cast(DataType::Utf8)])
    .distinct(None, DistinctKeepStrategy::First)
    .groupby([col(query_id_col)])
    .agg([col(subject_id_col).list()])
}

fn group_blast_hits_with_weights(
    hits: LazyFrame,
    query_id_col: &str,
    subject_id_col: &str,
    weight_col: &str, // may not be present, it's ok
    weight_col_agg: impl FnOnce(Expr) -> Expr,
) -> LazyFrame {
    let zipped_subject_col = format!("{subject_id_col}/{weight_col}");

    hits.filter(
        col(query_id_col)
            .is_not_null()
            .and(col(subject_id_col).is_not_null())
            .and(col(weight_col).is_not_null()),
    )
    .select([
        col(query_id_col),
        col(subject_id_col).cast(DataType::Utf8),
        col(weight_col),
    ])
    .groupby([col(query_id_col), col(subject_id_col)])
    .agg([weight_col_agg(col(weight_col))])
    .select([
        col(query_id_col),
        col(subject_id_col)
            .cast(DataType::Utf8)
            .map_many(
                |params| {
                    Ok(itertools::izip!(params[0].utf8()?, params[1].f32()?)
                        .map(|(id, w)| format!("{}/{}", id.unwrap_or(""), w.unwrap_or(f32::NAN)))
                        .collect())
                },
                &[col(weight_col)],
                GetOutput::from_type(DataType::Utf8),
            )
            .alias(&zipped_subject_col),
    ])
    .groupby([col(query_id_col)])
    .agg([col(&zipped_subject_col).list()])
}

pub fn preprocess_blastout(args: PreprocessBlastOutArgs) -> anyhow::Result<()> {
    let df = {
        let format = parse_blast_outfmt(&args.blast_outfmt)
            .context("Error parsing blast outfmt specifier")?;

        let df = load_blastout(args.blastout_path, &format, None)
            .context("Error loading blast output")?;

        eprintln!("Loaded BLAST+ output with schema {:?}", df.schema());

        let filters: Vec<_> = args
            .filter
            .iter()
            .map(|f| BlastOutFilter::parse_filter(f))
            .try_collect()?;

        apply_filters(df, filters)
    };

    let mut grouped = if let Some(weight_col) = &args.weight_col {
        let weight_col_agg = args
            .weight_col_agg
            .unwrap_or(WeightColAgg::Max)
            .into_expr_fn();

        group_blast_hits_with_weights(
            df,
            &args.query_id_col,
            &args.subject_id_col,
            weight_col,
            weight_col_agg,
        )
        .collect()?
    } else {
        group_blast_hits(df, &args.query_id_col, &args.subject_id_col).collect()?
    };

    for col in grouped.get_columns_mut() {
        if let DataType::List(inner) = col.dtype() {
            let inner_is_utf8 = matches!(**inner, DataType::Utf8);

            if !inner_is_utf8 {
                *col = col.cast(&DataType::List(Box::new(DataType::Utf8)))?;
            }

            *col = col.list()?.lst_join(";")?.into_series();
        }
    }

    writing_new_file_or_stdout!(&args.output, writer => {
        CsvWriter::new(writer)
            .with_delimiter(b'\t')
            .finish(&mut grouped)?;
    });

    Ok(())
}
