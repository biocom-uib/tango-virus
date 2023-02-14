#![feature(
    associated_type_defaults,
    hash_raw_entry,
    hash_set_entry,
    iterator_try_reduce,
    once_cell,
    never_type
)]

use clap::{Parser, Subcommand};

pub mod taxonomy;
pub(crate) mod preprocessed_taxonomy;

pub(crate) mod filter;
pub(crate) mod util;

mod crispr;
mod get_lineage;
mod preprocess_blastout;
mod preprocess_taxonomy;
mod refine_vpf_class;
mod tango_assign;

#[derive(Subcommand)]
enum Commands {
    PreprocessTaxonomy(preprocess_taxonomy::PreprocessTaxonomyArgs),
    GetLineage(get_lineage::GetLineageArgs),
    PreprocessBlastout(preprocess_blastout::PreprocessBlastOutArgs),
    TangoAssign(tango_assign::TangoAssignArgs),
    RefineVpfClass(refine_vpf_class::RefineVpfClassArgs),
}

/// METEOR: Metagenome and Metavirome Joint Analysis
#[derive(Parser)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    match args.command {
        Commands::PreprocessTaxonomy(args) => {
            preprocess_taxonomy::preprocess_taxonomy(args)?;
        }
        Commands::GetLineage(args) => {
            get_lineage::get_lineage(args)?;
        }
        Commands::PreprocessBlastout(args) => {
            preprocess_blastout::preprocess_blastout(args)?;
        }
        Commands::TangoAssign(args) => {
            tango_assign::tango_assign(args)?;
        }
        Commands::RefineVpfClass(args) => {
            refine_vpf_class::refine_vpf_class(args)?;
        }
    }

    Ok(())
}
