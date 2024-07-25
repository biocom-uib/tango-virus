use anyhow::Context;
use clap::{Args, ValueEnum};

use polars::datatypes::DataType;
use polars::lazy::prelude::Expr;
use polars::prelude::{
    col, CsvWriter, GetOutput, LazyFrame, SerWriter, UniqueKeepStrategy
};
use polars::series::Series;

use crate::tool::blast::blastout;
use crate::util::{self, filter::FromStrFilter, writing_new_file_or_stdout};


#[derive(ValueEnum, Debug, Default, Copy, Clone)]
pub enum WeightColAgg {
    Max,
    Mean,
    Median,
    Min,
    #[default]
    Count,
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

fn aggregate_subjects(subject_id: Expr) -> Expr {
    subject_id.str().join(";", true)
}

fn group_blast_hits(hits: LazyFrame, query_id_col: &str, subject_id_col: &str) -> LazyFrame {
    hits.filter(
        col(query_id_col)
            .is_not_null()
            .and(col(subject_id_col).is_not_null()),
    )
    .select([
        col(query_id_col),
        col(subject_id_col).cast(DataType::String)
    ])
    .unique(None, UniqueKeepStrategy::First)
    .group_by([col(query_id_col)])
    .agg([aggregate_subjects(col(subject_id_col))])
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
        col(subject_id_col).cast(DataType::String),
        col(weight_col),
    ])
    .group_by([col(query_id_col), col(subject_id_col)])
    .agg([weight_col_agg(col(weight_col))])
    .select([
        col(query_id_col),
        col(subject_id_col)
            .map_many(
                |params: &mut [Series]| {
                    let ids = params[0].str()?;

                    let ws = params[1].cast(&DataType::String)?;
                    let ws = ws.str()?;

                    let result = itertools::izip!(ids, ws)
                        .map(|(id, w)| format!("{}/{}", id.unwrap_or(""), w.unwrap_or("")))
                        .collect::<Series>();

                    Ok(Some(result))
                },
                &[col(weight_col)],
                GetOutput::from_type(DataType::String),
            )
            .alias(&zipped_subject_col),
    ])
    .group_by([col(query_id_col)])
    .agg([aggregate_subjects(col(&zipped_subject_col))])
}

pub fn preprocess_blastout(args: PreprocessBlastOutArgs) -> anyhow::Result<()> {
    let df = {
        let format = args
            .blast_outfmt
            .parse()
            .context("Error parsing blast outfmt specifier")?;

        let mut df = blastout::load_blastout(&args.blastout_path, &format, None)
            .context("Error loading blast output")?;

        eprintln!("Loaded BLAST+ output with schema {:?}", df.schema());

        let filters = blastout::BlastOutFilter::parse_filters(&args.filter)
            .context("Error parsing --filter")?;

        blastout::apply_filters(df, filters)
    };

    let grouped = if let Some(weight_col) = &args.weight_col {
        let weight_col_agg = args.weight_col_agg.into_expr_fn();

        group_blast_hits_with_weights(
            df,
            &args.query_id_col,
            &args.subject_id_col,
            weight_col,
            weight_col_agg,
        )
    } else {
        group_blast_hits(df, &args.query_id_col, &args.subject_id_col)
    };

    let mut grouped = grouped.collect()?;

    //for col_name in grouped.get_column_names_owned() {
    //    if let DataType::List(inner_dtype) = grouped.column(&col_name)?.dtype() {
    //        let inner_dtype_is_string = inner_dtype.is_string();

    //        grouped.try_apply(&col_name, |series| {
    //            let mut series = Cow::Borrowed(series);

    //            if !inner_dtype_is_string {
    //                series = Cow::Owned(series.cast(&DataType::List(Box::new(DataType::String)))?);
    //            }

    //            Ok(series.list()?.join_literal(";", true)?.into_series())
    //        })?;
    //    }
    //}

    writing_new_file_or_stdout!(&args.output, writer => {
        let writer = writer.context("Error creating grouped hits file")?;

        util::ignore_broken_pipe(
            CsvWriter::new(writer)
                .with_separator(b'\t')
                .finish(&mut grouped)
        )?;
    });

    Ok(())
}
