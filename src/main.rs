#![feature(
    associated_type_defaults,
    hash_raw_entry,
    hash_set_entry,
    iterator_try_reduce,
    never_type,
    result_flattening,
)]

use clap::{Parser, Subcommand};

pub mod taxonomy;
pub(crate) mod preprocessed_taxonomy;

pub(crate) mod util;

mod crispr_match;
mod fetch;
#[cfg(feature = "ppi")]
mod ppin;
mod get_lineage;
mod preprocess_blastout;
mod preprocess_taxonomy;
mod refine_vpf_class;
mod tango_assign;
mod tool;

#[derive(Subcommand)]
enum Commands {
    Fetch(fetch::FetchArgs),
    PreprocessTaxonomy(preprocess_taxonomy::PreprocessTaxonomyArgs),
    PreprocessBlastout(preprocess_blastout::PreprocessBlastOutArgs),
    TangoAssign(tango_assign::TangoAssignArgs),
    RefineVpfClass(refine_vpf_class::RefineVpfClassArgs),
    CrisprMatch(crispr_match::CrisprMatchArgs),
    GetLineage(get_lineage::GetLineageArgs),
    #[cfg(feature = "ppi")]
    Ppi(ppin::PpiArgs),
}

/// METEOR: Metagenome and Metavirome Joint Analysis
#[derive(Parser)]
#[command(author, version, about)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    match args.command {
        Commands::Fetch(args) => {
            fetch::fetch(args)?;
        }
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
        #[cfg(feature = "ppi")]
        Commands::Ppi(args) => {
            ppin::ppi(args)?;
        }
        Commands::CrisprMatch(args) => {
            crispr_match::crispr_match(args)?;
        }
    }

    Ok(())
}
