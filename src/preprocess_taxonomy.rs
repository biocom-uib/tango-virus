use std::{path::{PathBuf, Path}, fs::File};

use anyhow::anyhow;
use clap::{ArgEnum, Args};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::taxonomy::{
    formats::{
        ncbi::{self, NcbiTaxonomy, AllNames},
        newick::NewickTaxonomy,
    },
    TaxonomyMut,
};

#[derive(ArgEnum, Debug, Copy, Clone)]
pub enum TaxonomyFormat {
    NCBI,
    Newick,
}

#[derive(Deserialize, Serialize)]
pub enum SomeTaxonomy {
    NcbiTaxonomyWithSingleClassNames(NcbiTaxonomy<ncbi::SingleClassNames>),
    NcbiTaxonomyWithManyNames(NcbiTaxonomy<ncbi::AllNames>),
    NewickTaxonomy(NewickTaxonomy),
}

#[derive(Deserialize, Serialize)]
pub struct PreprocessedTaxonomy {
    pub tree: SomeTaxonomy,
}

#[derive(ArgEnum, Debug, Copy, Clone)]
pub enum PreprocessedTaxonomyFormat {
    CBOR,
    JSON,
}

impl PreprocessedTaxonomy {
    pub fn deserialize_with_format<P: AsRef<Path>>(
        path: P,
        format: PreprocessedTaxonomyFormat,
    ) -> anyhow::Result<Self> {
        match format {
            PreprocessedTaxonomyFormat::CBOR =>
                Ok(serde_cbor::from_reader(File::open(path)?)?),
            PreprocessedTaxonomyFormat::JSON =>
                Ok(serde_json::from_reader(File::open(path)?)?),
        }
    }

    pub fn serialize_with_format<P: AsRef<Path>>(
        &self,
        path: P,
        format: PreprocessedTaxonomyFormat,
    ) -> anyhow::Result<()> {
        match format {
            PreprocessedTaxonomyFormat::CBOR =>
                Ok(serde_cbor::to_writer(File::create(path)?, self)?),
            PreprocessedTaxonomyFormat::JSON =>
                Ok(serde_json::to_writer(File::create(path)?, self)?),
        }
    }
}

const DEFAULT_CONTRACTION_RANKS: &'static str = "superkingdom,phylum,class,order,family,genus,species";

#[derive(Args)]
pub struct PreprocessTaxonomyArgs {
    /// Enables contraction. If provided, RANKS should be a comma-separated list of the ranks to keep after contraction.
    #[clap(short, long, value_name = "RANKS")]
    contract: Option<Option<String>>,

    /// Input taxonomy format
    #[clap(long, arg_enum, default_value_t = TaxonomyFormat::NCBI)]
    input_format: TaxonomyFormat,

    /// Path to the input taxonomy (in the case of NCBI, specify the extracted directory of taxdump)
    input_taxonomy: String,

    /// Comma-separated list of name classes to load (if INPUT_FORMAT is ncbi)
    #[clap(long, default_value = "scientific name,synonym")]
    ncbi_name_classes: String,

    /// Comma-separated list of rank names to associate to each level of the Newick taxonomy (if
    /// INPUT_FORMAT is newick)
    #[clap(long)]
    newick_ranks: Option<String>,

    /// Output taxonomy format
    #[clap(long, arg_enum, default_value_t = PreprocessedTaxonomyFormat::CBOR)]
    output_format: PreprocessedTaxonomyFormat,

    /// Output taxonomy path
    output_taxonomy: String,
}

impl PreprocessTaxonomyArgs {
    pub fn get_contraction_ranks(&self) -> Option<Vec<&str>> {
        self.contract.as_ref().map(|ranks| {
            ranks
                .as_deref()
                .unwrap_or(DEFAULT_CONTRACTION_RANKS)
                .split(',')
                .map(str::trim)
                .collect()
        })
    }
}

fn preprocess_ncbi(args: PreprocessTaxonomyArgs) -> anyhow::Result<PreprocessedTaxonomy> {
    let mut path = PathBuf::from(&args.input_taxonomy);

    let name_classes = args.ncbi_name_classes.split(',').collect_vec();

    let mut taxo = {
        path.push("nodes.dmp");
        let taxo = NcbiTaxonomy::load_nodes(&path)?;
        path.pop();
        eprintln!("Loaded nodes.dmp");

        path.push("merged.dmp");
        let taxo = taxo.with_merged_taxids(&path)?;
        path.pop();
        eprintln!("Loaded merged.dmp");

        path.push("names.dmp");
        let names = AllNames::load_filtered_names_dmp(&name_classes, &path)?;
        path.pop();
        eprintln!("Loaded names.dmp with classes {name_classes:?}");

        taxo.with_names(names)
    };

    if let Some(ranks) = &args.get_contraction_ranks() {
        taxo.contract(&ranks);
        todo!("Review contraction, seems broken: see 208964 before/after");
    }

    Ok(PreprocessedTaxonomy {
        tree: if name_classes.len() == 1 {
            match std::mem::take(&mut taxo.names).only_of_class(name_classes[0]) {
                Ok(names) => SomeTaxonomy::NcbiTaxonomyWithSingleClassNames(taxo.with_names(names)),
                Err(names) => SomeTaxonomy::NcbiTaxonomyWithManyNames(taxo.with_names(names)),
            }
        } else {
            SomeTaxonomy::NcbiTaxonomyWithManyNames(taxo)
        }
    })
}

fn preprocess_newick(args: PreprocessTaxonomyArgs) -> anyhow::Result<PreprocessedTaxonomy> {
    let ranks = args
        .newick_ranks
        .as_ref()
        .ok_or_else(|| anyhow!("Flag --newick-ranks is required for --input-format newick"))?
        .split(',')
        .map(|s| s.trim().to_owned())
        .collect_vec();

    let taxo = NewickTaxonomy::load_newick(&args.input_taxonomy, ranks)?;

    if let Some(_ranks) = &args.get_contraction_ranks() {
        anyhow!("Contraction is not supported for Newick taxonomies");
    }

    Ok(PreprocessedTaxonomy {
        tree: SomeTaxonomy::NewickTaxonomy(taxo)
    })
}

pub fn preprocess_taxonomy(args: PreprocessTaxonomyArgs) -> anyhow::Result<()> {
    let output_taxonomy = args.output_taxonomy.clone();
    let output_format = args.output_format;

    let taxonomy = match args.input_format {
        TaxonomyFormat::NCBI => preprocess_ncbi(args),
        TaxonomyFormat::Newick => preprocess_newick(args),
    }?;

    taxonomy.serialize_with_format(&output_taxonomy, output_format)?;

    Ok(())
}
