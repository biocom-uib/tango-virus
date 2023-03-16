use anyhow::Context;
use clap::{Args, ValueEnum};

use polars::datatypes::DataType;
use polars::lazy::prelude::Expr;
use polars::prelude::{
    col, CsvWriter, GetOutput, IntoSeries, LazyFrame, ListNameSpaceImpl, SerWriter, UniqueKeepStrategy,
};

use crate::util::blastout;
use crate::util::filter::FromStrFilter;
use crate::util::writing_new_file_or_stdout;


#[derive(ValueEnum, Debug, Copy, Clone)]
pub enum WeightColAgg {
    Max,
    Mean,
    Median,
    Min,
    Count,
    Product,
    Sum,
}

impl Default for WeightColAgg {
    fn default() -> Self {
        WeightColAgg::Count
    }
}

impl WeightColAgg {
    fn into_expr_fn(self) -> fn(Expr) -> Expr {
        use WeightColAgg::*;

        match self {
            Max => Expr::max,
            Mean => Expr::mean,
            Median => Expr::median,
            Min => Expr::min,
            Count => Expr::count,
            Product => Expr::product,
            Sum => Expr::sum,
        }
    }
}

/// Load the output of NCBI's BLAST+ (blastn) to produce a suitable input file for the assign
/// subcommand.
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

    /// Aggregation function to apply to the selected weight column (if needed).
    #[clap(long, value_enum, default_value_t)]
    weight_col_agg: WeightColAgg,

    /// Filters to apply before processing. Example: --filter 'evalue<=1e-3'
    #[clap(long)]
    filter: Vec<String>,
}

fn group_blast_hits(hits: LazyFrame, query_id_col: &str, subject_id_col: &str) -> LazyFrame {
    hits.filter(
        col(query_id_col)
            .is_not_null()
            .and(col(subject_id_col).is_not_null()),
    )
    .select([col(query_id_col), col(subject_id_col).cast(DataType::Utf8)])
    .unique(None, UniqueKeepStrategy::First)
    .groupby([col(query_id_col)])
    .agg([col(subject_id_col).list()])
}

fn group_blast_hits_with_weights(
    hits: LazyFrame,
    query_id_col: &str,
    subject_id_col: &str,
    weight_col: &str,
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
        let format = args
            .blast_outfmt
            .parse()
            .context("Error parsing blast outfmt specifier")?;

        let df = blastout::load_blastout(&args.blastout_path, &format, None)
            .context("Error loading blast output")?;

        eprintln!("Loaded BLAST+ output with schema {:?}", df.schema());

        let filters =
            blastout::BlastOutFilter::parse_filters(&args.filter).context("Error parsing --filter")?;

        blastout::apply_filters(df, filters)
    };

    let mut grouped = if let Some(weight_col) = &args.weight_col {
        let weight_col_agg = args.weight_col_agg.into_expr_fn();

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
        let writer = writer.context("Error creating grouped hits file")?;

        CsvWriter::new(writer)
            .with_delimiter(b'\t')
            .finish(&mut grouped)?;
    });

    Ok(())
}
