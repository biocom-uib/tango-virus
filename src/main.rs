#![feature(
    absolute_path,
    associated_type_defaults,
    hash_raw_entry,
    hash_set_entry,
    iterator_try_reduce,
    never_type,
    once_cell,
    result_flattening,
)]

use clap::{Parser, Subcommand};

pub mod taxonomy;
pub(crate) mod preprocessed_taxonomy;

pub(crate) mod util;

mod crispr_match;
mod ppin;
mod get_lineage;
mod preprocess_blastout;
mod preprocess_taxonomy;
mod refine_vpf_class;
mod tango_assign;

#[derive(Subcommand)]
enum Commands {
    PreprocessTaxonomy(preprocess_taxonomy::PreprocessTaxonomyArgs),
    PreprocessBlastout(preprocess_blastout::PreprocessBlastOutArgs),
    TangoAssign(tango_assign::TangoAssignArgs),
    RefineVpfClass(refine_vpf_class::RefineVpfClassArgs),
    CrisprMatch(crispr_match::CrisprMatchArgs),
    GetLineage(get_lineage::GetLineageArgs),
    PrepareUniprotBlastdb(ppin::prepare_uniprot_blastdb::PrepareUniProtBlastDBArgs),
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
        Commands::PrepareUniprotBlastdb(args) => {
            ppin::prepare_uniprot_blastdb::prepare_uniprot_blastdb(args)?;
        }
        Commands::CrisprMatch(args) => {
            crispr_match::crispr_match(args)?;
        }
    }

    Ok(())
}
