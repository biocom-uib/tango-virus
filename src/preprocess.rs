use std::{path::PathBuf, collections::HashMap, fs::File};

use anyhow::anyhow;
use clap::{ArgEnum, Args};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::taxonomy::{formats::{ncbi::{self, NcbiTaxonomy}, newick::{NewickTaxonomy, load_newick_simple_tree, self}}, NodeId};

#[derive(ArgEnum, Debug, Copy, Clone)]
pub enum TaxonomyFormat {
    NCBI,
    Newick,
}

#[derive(Deserialize, Serialize)]
pub enum SomeTaxonomy {
    NcbiTaxonomyWithNames(NcbiTaxonomy<ncbi::NoNames>),
    NcbiTaxonomyWithSingleClassNames(NcbiTaxonomy<ncbi::SingleClassNames>),
    NewickTaxonomy(NewickTaxonomy),
}

#[derive(Deserialize, Serialize)]
pub struct PreprocessedTaxonomy {
    tree: SomeTaxonomy,
}

#[derive(ArgEnum, Debug, Copy, Clone)]
pub enum PreprocessedTaxonomyFormat {
    CBOR,
    JSON,
}

const DEFAULT_CONTRACTION_RANKS: &'static str = "superkingdom,phylum,class,order,family,genus,species";

#[derive(Args)]
pub struct PreprocessArgs {
    /// Enables contraction. If provided, RANKS should be a comma-separated list of the ranks to keep after contraction.
    #[clap(short, long, value_name = "RANKS", default_value = DEFAULT_CONTRACTION_RANKS)]
    contract: Option<Option<String>>,

    /// Input taxonomy format
    #[clap(long, arg_enum, default_value_t = TaxonomyFormat::NCBI)]
    input_format: TaxonomyFormat,

    /// Path to the input taxonomy (in the case of NCBI, specify the extracted directory of taxdump)
    input_taxonomy: String,

    /// Rank names to associate to each level of the Newick taxonomy (if INPUT_FORMAT is newick)
    #[clap(long)]
    newick_ranks: Option<String>,

    /// Output taxonomy format
    #[clap(long, arg_enum, default_value_t = PreprocessedTaxonomyFormat::CBOR)]
    output_format: PreprocessedTaxonomyFormat,

    /// Output taxonomy path
    output_taxonomy: String,
}

impl PreprocessArgs {
    pub fn get_contraction_ranks(&self) -> Option<Vec<&str>> {
        self.contract.map(|ranks| {
            ranks.as_deref().unwrap_or(DEFAULT_CONTRACTION_RANKS).split(',').map(str::trim).collect()
        })
    }
}

fn preprocess_ncbi(args: PreprocessArgs) -> anyhow::Result<()> {
    let mut path = PathBuf::from(args.input_taxonomy);
    path.push("nodes.dmp");

    let mut taxo = NcbiTaxonomy::load_nodes(&path)?;

    if let Some(ranks) = &args.get_contraction_ranks() {
        taxo.contract(&ranks);
    }

    Ok(())
}

fn preprocess_newick(args: PreprocessArgs) -> anyhow::Result<()> {
    let ranks = args
        .newick_ranks
        .ok_or_else(|| anyhow!("Flag --newick-ranks is required for --input-format newick"))?
        .split(',')
        .map(|s| s.trim().to_owned())
        .collect_vec();

    let mut taxo = NewickTaxonomy::load_newick(args.input_taxonomy.as_ref(), ranks)?;

    if let Some(ranks) = &args.get_contraction_ranks() {
        taxo.contract(ranks);
    }

    Ok(())
}

pub fn preprocess(args: PreprocessArgs) -> anyhow::Result<()> {
    let taxonomy = match args.input_format {
        TaxonomyFormat::NCBI => preprocess_ncbi(args),
        TaxonomyFormat::Newick => preprocess_newick(args),
    };

    Ok(())
}
