use std::{path::{PathBuf, Path}, fs::File};

use anyhow::{anyhow, Context};
use clap::{ArgEnum, Args};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::taxonomy::{
    formats::{
        ncbi::{self, AllNames, NamesAssoc, NcbiTaxonomy},
        newick::NewickTaxonomy,
    },
    Taxonomy, TaxonomyMut,
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
    pub ordered_ranks: Option<Vec<String>>,
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

fn rank_syms_to_strings<'a, Tax: Taxonomy>(tax: &Tax, syms: &Vec<Tax::RankSym>) -> anyhow::Result<Vec<String>> {
    syms.iter()
        .map(|&sym| tax.rank_sym_str(sym).map(ToOwned::to_owned))
        .collect::<Option<Vec<String>>>()
        .ok_or_else(|| anyhow!("Unable to translate ordered_ranks to Vec<String>. This is a bug"))

}

fn contract_and_order_ranks<Names: 'static + NamesAssoc + Send>(
    args: &PreprocessTaxonomyArgs,
    tax: &mut NcbiTaxonomy<Names>,
) -> anyhow::Result<Option<Vec<String>>> {
    let ordered_ranks = if let Some(ranks) = &args.get_contraction_ranks() {
        let contraction_ranks_syms: Vec<_> = ranks.iter()
            .map(|rank| tax.lookup_rank_sym(rank).ok_or_else(|| anyhow!("Unrecognized rank {rank}")))
            .try_collect()
            .context("Error during contraction")?;

        tax.contract(contraction_ranks_syms.iter().copied().collect());

        let rank_ordering = tax.rank_ordering();

        if let Ok(rank_ordering) = &rank_ordering {
            assert_eq!(
                rank_syms_to_strings(tax, rank_ordering)?,
                rank_syms_to_strings(tax, &contraction_ranks_syms)?,
            );
        }

        assert!(tax.topology_health_check());
        //todo!("Review contraction, seems broken: see 208964 before/after");

        rank_ordering
    } else {
        tax.rank_ordering()
    };

    let ordered_ranks = match ordered_ranks {
        Ok(ordered_ranks) => {
            Some(
                rank_syms_to_strings(tax, &ordered_ranks)
                    .context("Error mapping rank symbols back to strings. This is bug")?
            )
        },

        Err(topsort_loop) => {
            eprintln!(
                "Warning: the taxonomy contains rank loops: {:?}",
                rank_syms_to_strings(tax, &topsort_loop.0),
            );
            None
        }
    };

    Ok(ordered_ranks)
}

fn preprocess_ncbi(args: PreprocessTaxonomyArgs) -> anyhow::Result<PreprocessedTaxonomy> {
    let mut path = PathBuf::from(&args.input_taxonomy);

    let name_classes = args.ncbi_name_classes.split(',').collect_vec();

    let mut tax = {
        path.push("nodes.dmp");
        let tax = NcbiTaxonomy::load_nodes(&path)?;
        path.pop();
        eprintln!("Loaded nodes.dmp");

        path.push("merged.dmp");
        let tax = tax.with_merged_taxids(&path)?;
        path.pop();
        eprintln!("Loaded merged.dmp");

        path.push("names.dmp");
        let names = AllNames::load_filtered_names_dmp(&name_classes, &path)?;
        path.pop();
        eprintln!("Loaded names.dmp with classes {name_classes:?}");

        tax.with_names(names)
    };

    let ordered_ranks = contract_and_order_ranks(&args, &mut tax)?;

    Ok(PreprocessedTaxonomy {
        tree: if name_classes.len() == 1 {
            match std::mem::take(&mut tax.names).only_of_class(name_classes[0]) {
                Ok(names) => SomeTaxonomy::NcbiTaxonomyWithSingleClassNames(tax.with_names(names)),
                Err(names) => SomeTaxonomy::NcbiTaxonomyWithManyNames(tax.with_names(names)),
            }
        } else {
            SomeTaxonomy::NcbiTaxonomyWithManyNames(tax)
        },
        ordered_ranks,
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

    let tax = NewickTaxonomy::load_newick(&args.input_taxonomy, ranks.clone())?;

    if let Some(_ranks) = &args.get_contraction_ranks() {
        anyhow!("Contraction is not supported for Newick taxonomies");
    }

    Ok(PreprocessedTaxonomy {
        tree: SomeTaxonomy::NewickTaxonomy(tax),
        ordered_ranks: Some(ranks),
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
