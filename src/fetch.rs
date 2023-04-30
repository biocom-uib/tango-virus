use clap::{Args, Subcommand};

pub mod ncbi_taxonomy;
#[cfg(feature = "ppi")]
pub mod string_viruses;
#[cfg(feature = "ppi")]
pub mod uniprot;

#[derive(Subcommand)]
enum Commands {
    NcbiTaxonomy(ncbi_taxonomy::NcbiTaxonomyFetchArgs),
    #[cfg(feature = "ppi")]
    StringViruses(string_viruses::StringVirusesFetchArgs),
    #[cfg(feature = "ppi")]
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
        #[cfg(feature = "ppi")]
        Commands::StringViruses(args) => {
            string_viruses::fetch(&args)?;
        }
        #[cfg(feature = "ppi")]
        Commands::Uniprot(args) => {
            uniprot::fetch(&args)?;
        }
    }

    Ok(())
}

