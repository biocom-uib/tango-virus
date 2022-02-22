use std::{path::PathBuf, collections::HashMap};

use clap::{ArgEnum, Args};
use itertools::Itertools;
use string_interner::{backend::{BufferBackend, Backend}, StringInterner};

use crate::taxonomy::{formats::ncbi::NcbiTaxonomy, NodeId};

#[derive(ArgEnum, Debug, Copy, Clone)]
pub enum TaxonomyFormat {
    NCBI,
    Newick,
}

#[derive(ArgEnum, Debug, Copy, Clone)]
pub enum PreprocessedTaxonomyFormat {
    CBOR,
    JSON,
}

#[derive(Args)]
pub struct PreprocessArgs {
    /// Enables contraction.
    #[clap(short, long)]
    contract: bool,

    /// Comma-separated list of the ranks to keep after contraction.
    #[clap(long, value_name = "RANKS", default_value = "superkingdom,phylum,class,order,family,genus,species")]
    contraction_ranks: String,

    /// Input taxonomy format
    #[clap(long, arg_enum, default_value_t = TaxonomyFormat::NCBI)]
    input_format: TaxonomyFormat,

    /// Path to the input taxonomy (in the case of NCBI, specify the extracted directory of taxdump)
    input_taxonomy: String,

    /// Output taxonomy format
    #[clap(long, arg_enum, default_value_t = PreprocessedTaxonomyFormat::CBOR)]
    output_format: PreprocessedTaxonomyFormat,

    /// Output taxonomy path
    output_taxonomy: String,
}

trait AssocNodeContainer {
    type Item;

    fn get_node(&self, node: NodeId) -> Option<&Self::Item>;
    fn get_node_mut(&mut self, node: NodeId) -> Option<&mut Self::Item>;
}

trait AssocNodeStorage {
    type Container<T>: AssocNodeContainer<Item = T>;
}

struct CompactStorage;

impl<T> AssocNodeContainer for Vec<T> {
    type Item = T;

    fn get_node(&self, node: NodeId) -> Option<&T> {
        self.get(node.0)
    }

    fn get_node_mut(&mut self, node: NodeId) -> Option<&mut T> {
        self.get_mut(node.0)
    }
}

impl AssocNodeStorage for CompactStorage {
    type Container<T> = Vec<T>;
}

struct SparseStorage;

impl<T> AssocNodeContainer for HashMap<NodeId, T> {
    type Item = T;

    fn get_node(&self, node: NodeId) -> Option<&T> {
        self.get(&node)
    }

    fn get_node_mut(&mut self, node: NodeId) -> Option<&mut T> {
        self.get_mut(&node)
    }
}

impl AssocNodeStorage for SparseStorage {
    type Container<T> = HashMap<NodeId, T>;
}

pub type Symbol = <BufferBackend as Backend>::Symbol;

struct PreprocessedTaxonomy<Storage: AssocNodeStorage> {
    parent_ids: Storage::Container<NodeId>,
    children_lookup: Storage::Container<Vec<NodeId>>,

    ranks: Storage::Container<Symbol>,
    ranks_interner: StringInterner<BufferBackend>,

    preorder: Storage::Container<usize>,
}

enum SomePreprocessedTaxonomy {
    Compact(PreprocessedTaxonomy<CompactStorage>),
    Sparse(PreprocessedTaxonomy<SparseStorage>),
}

fn preprocess_ncbi(args: PreprocessArgs) -> anyhow::Result<()> {
    let mut path = PathBuf::from(args.input_taxonomy);
    path.push("nodes.dmp");

    let mut taxo = NcbiTaxonomy::load_nodes(&path)?;

    if args.contract {
        let ranks = args.contraction_ranks.split(',').collect_vec();
        taxo.contract(&ranks);
    }

    let preorder: HashMap<NodeId, usize> = taxo
        .preorder_descendants(taxo.root)
        .enumerate()
        .map(|(i, n)| (n, i))
        .collect();

    Ok(())
}

fn preprocess_newick(args: PreprocessArgs) -> anyhow::Result<()> {
    Ok(())
}

pub fn preprocess(args: PreprocessArgs) -> anyhow::Result<()> {
    let taxonomy = match args.input_format {
        TaxonomyFormat::NCBI => preprocess_ncbi(args),
        TaxonomyFormat::Newick => preprocess_newick(args),
    };

    Ok(())
}
