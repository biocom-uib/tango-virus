use clap::{Args, Subcommand};

pub mod ncbi_taxonomy;
pub mod string_viruses;
pub mod uniprot;

#[derive(Subcommand)]
enum Commands {
    NcbiTaxonomy(ncbi_taxonomy::NcbiTaxonomyFetchArgs),
    StringViruses(string_viruses::StringVirusesFetchArgs),
    Uniprot(uniprot::UniProtFetchArgs),
}

/// Download external files to be used by Meteor.
#[derive(Args)]
pub struct FetchArgs {
    #[clap(subcommand)]
    command: Commands,
}

pub fn fetch(args: FetchArgs) -> anyhow::Result<()> {
    match args.command {
        Commands::NcbiTaxonomy(args) => {
            ncbi_taxonomy::fetch(&args)?;
        }
        Commands::StringViruses(args) => {
            string_viruses::fetch(&args)?;
        }
        Commands::Uniprot(args) => {
            uniprot::fetch(&args)?;
        }
    }

    Ok(())
}

