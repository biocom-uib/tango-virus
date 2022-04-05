#![feature(
    associated_type_bounds,
    associated_type_defaults,
    generic_associated_types,
    io_read_to_string,
    iterator_try_reduce,
    label_break_value,
    never_type,
)]

use clap::{Parser, Subcommand};

pub mod taxonomy;

mod preprocess_blastout;
mod preprocess_taxonomy;
mod assign;
mod refine_vpf_class;

#[derive(Subcommand)]
enum Commands {
    PreprocessTaxonomy(preprocess_taxonomy::PreprocessTaxonomyArgs),
    PreprocessBlastout(preprocess_blastout::PreprocessBlastOutArgs),
    Assign(assign::AssignArgs),
    RefineVpfClass(refine_vpf_class::RefineVpfClassArgs),
}

/// Rust port of TANGO: Taxonomic Assignment in Metagenomics
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
        Commands::PreprocessBlastout(args) => {
            preprocess_blastout::preprocess_blastout(args)?;
        }
        Commands::Assign(args) => {
            assign::assign(args)?;
        }
        Commands::RefineVpfClass(args) => {
            refine_vpf_class::refine_vpf_class(args)?;
        }
    }

    Ok(())
}
