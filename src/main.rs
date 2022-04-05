#![feature(
    associated_type_bounds,
    associated_type_defaults,
    generic_associated_types,
    io_read_to_string,
    iterator_try_reduce,
    label_break_value,
    never_type,
    stdio_locked,
    type_alias_impl_trait
)]

use clap::{Parser, Subcommand};

pub mod taxonomy;
mod preprocess_blastout;
mod preprocess_taxonomy;
mod assign;

#[derive(Subcommand)]
enum Commands {
    PreprocessTaxonomy(preprocess_taxonomy::PreprocessTaxonomyArgs),
    PreprocessBlastout(preprocess_blastout::PreprocessBlastOutArgs),
    Assign(assign::AssignArgs),
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
    }

    Ok(())
}
