#![feature(
    generic_associated_types,
    never_type,
    io_read_to_string,
    type_alias_impl_trait,
    label_break_value,
)]

use anyhow::{Error, bail};
use clap::{Parser, Args, Subcommand};

pub mod taxonomy;
mod preprocess;

#[derive(Args)]
struct AssignArgs {}

#[derive(Subcommand)]
enum Commands {
    Preprocess(preprocess::PreprocessArgs),
    Assign(AssignArgs),
}

/// Rust port of TANGO: Taxonomic Assignment in Metagenomics
#[derive(Parser)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

fn main() -> Result<(), Error> {
    let args = Cli::parse();

    match args.command {
        Commands::Preprocess(args) => {
            preprocess::preprocess(args)?;
        },
        Commands::Assign(AssignArgs {}) => {
            bail!("assign is not yet implemented");
        },
    }

    Ok(())
}
