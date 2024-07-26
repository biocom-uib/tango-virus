use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path, iter::Flatten,
};

use crate::taxonomy::{
    self, generic::GenericTaxonomy, ContractedNodes, LabelledTaxonomy, NodeId, Taxonomy,
    TaxonomyMut, TopologyReplacer,
};

use newick_rs::SimpleTree;
use serde::{Deserialize, Serialize};
use string_interner::backend::{Backend, StringBackend};

use self::dmp::parse_fields;
pub use self::{
    dmp::DmpError,
    names::{AllNames, NamesAssoc, NoNames, SingleClassNames},
};

pub mod dmp;
pub mod names;

pub type TaxId = NodeId;

const FORMAT_VERSION: u32 = GenericTaxonomy::FORMAT_VERSION;

#[derive(Deserialize, Serialize)]
pub struct NcbiTaxonomy<Names> {
    tree: GenericTaxonomy,

    // old to new
    merged_taxids: Option<HashMap<TaxId, TaxId>>,
    contracted_taxids: Option<ContractedNodes>,

    pub(crate) names: Names,
}

impl<Names> NcbiTaxonomy<Names> {
    pub fn node_count(&self) -> usize {
        1 /* root */ + self.tree.parent_ids.len()
    }

    pub fn with_merged_taxids(self, merged_dmp: &Path) -> Result<Self, DmpError> {
        let mut merged = HashMap::new();

        for line in BufReader::new(File::open(merged_dmp)?).lines() {
            let line = line?;

            parse_fields! {
                &line => {
                    let old_taxid  = next as "old_tax_id";
                    let new_taxid  = next as "new_tax_id";
                };

                let old_taxid = NodeId(old_taxid.parse()?);
                let new_taxid = NodeId(new_taxid.parse()?);

                merged.insert(old_taxid, new_taxid);
            }
        }

        Ok(Self {
            merged_taxids: Some(merged),
            ..self
        })
    }

    pub fn topology_health_check(&self) -> bool {
        self.tree.topology_health_check()
    }

    pub fn with_names<NewNames>(self, new_names: NewNames) -> NcbiTaxonomy<NewNames> {
        NcbiTaxonomy {
            tree: self.tree,
            merged_taxids: self.merged_taxids,
            contracted_taxids: self.contracted_taxids,
            names: new_names,
        }
    }

    pub fn map_names<F, NewNames>(self, f: F) -> NcbiTaxonomy<NewNames>
    where
        F: FnOnce(Names) -> NewNames,
    {
        NcbiTaxonomy {
            tree: self.tree,
            merged_taxids: self.merged_taxids,
            contracted_taxids: self.contracted_taxids,
            names: f(self.names),
        }
    }

    pub fn without_names(self) -> NcbiTaxonomy<NoNames> {
        self.with_names(NoNames {})
    }
}

impl NcbiTaxonomy<NoNames> {
    pub const FORMAT_VERSION: u32 = FORMAT_VERSION + NoNames::FORMAT_VERSION;

    pub fn load_nodes<P: AsRef<Path>>(nodes_dmp: P) -> Result<Self, DmpError> {
        let mut builder = GenericTaxonomy::builder();

        for line in BufReader::new(File::open(nodes_dmp)?).lines() {
            let line = line?;

            parse_fields! {
                &line => {
                    let taxid  = next as "tax id";
                    let ptaxid = next as "parent tax_id";
                    let rank   = next as "rank";
                };

                let taxid = NodeId(taxid.parse()?);
                let ptaxid = NodeId(ptaxid.parse()?);

                if rank != "no rank" {
                    builder.set_rank(taxid, rank);
                }

                if taxid == ptaxid {
                    builder.set_root(taxid)?;
                } else {
                    builder.insert_edge(ptaxid, taxid)?;
                }
            }
        }

        let tree = builder.build()?;

        Ok(NcbiTaxonomy {
            tree,
            merged_taxids: None,
            contracted_taxids: None,
            names: NoNames::new(),
        })
    }
}

impl NcbiTaxonomy<SingleClassNames> {
    pub const FORMAT_VERSION: u32 = FORMAT_VERSION + SingleClassNames::FORMAT_VERSION;
}

impl NcbiTaxonomy<AllNames> {
    pub const FORMAT_VERSION: u32 = FORMAT_VERSION + AllNames::FORMAT_VERSION;
}

impl<Names: 'static> Taxonomy for NcbiTaxonomy<Names> {
    fn get_root(&self) -> TaxId {
        self.tree.get_root()
    }

    // try to fix bad nodes into actual nodes either through merged NCBI taxids or known contracted
    // nodes
    fn fixup_node(&self, taxid: u64) -> Option<NodeId> {
        let mut node = taxid;

        let merged = self.merged_taxids.as_ref();
        let get_merged = |n| merged.and_then(|m| m.get(&NodeId(n)));

        let contracted = self.contracted_taxids.as_ref();
        let get_contracted = |n| contracted.and_then(|c| c.0.get(&NodeId(n)));

        loop {
            if let Some(merged_node) = get_merged(node) {
                node = merged_node.0;
            } else if let Some(contracted_node) = get_contracted(node) {
                node = contracted_node.0;
            } else if let Some(tree_node) = self.tree.fixup_node(node) {
                return Some(tree_node);
            } else {
                return None;
            }
        }
    }

    fn is_leaf(&self, node: TaxId) -> bool {
        self.tree.is_leaf(node)
    }

    fn has_uniform_depths(&self) -> Option<usize> {
        self.tree.has_uniform_depths()
    }

    type Children<'a> = <GenericTaxonomy as Taxonomy>::Children<'a>;

    fn iter_children(&self, node: TaxId) -> Self::Children<'_> {
        self.tree.iter_children(node)
    }

    fn find_parent(&self, node: TaxId) -> Option<TaxId> {
        self.tree.find_parent(node)
    }

    type RankSym = <StringBackend as Backend>::Symbol;

    fn rank_sym_str(&self, rank_sym: Self::RankSym) -> Option<&str> {
        self.tree.rank_sym_str(rank_sym)
    }

    fn lookup_rank_sym(&self, rank: &str) -> Option<Self::RankSym> {
        self.tree.lookup_rank_sym(rank)
    }

    fn find_rank(&self, node: TaxId) -> Option<Self::RankSym> {
        self.tree.find_rank(node)
    }

    type NodeRanks<'a> = <GenericTaxonomy as Taxonomy>::NodeRanks<'a>;

    fn node_ranks(&self) -> Self::NodeRanks<'_> {
        self.tree.node_ranks()
    }
}

impl<Names: NamesAssoc + Send + 'static> LabelledTaxonomy for NcbiTaxonomy<Names> {
    type Labels<'a> = Flatten<<Option<Names::NamesLookupIter<'a>> as IntoIterator>::IntoIter>;

    fn labels_of(&self, node: NodeId) -> Self::Labels<'_> {
        self.names
            .lookup_names(node)
            .map(|lookup| Names::iter_lookup_names(lookup))
            .into_iter()
            .flatten()
    }

    type NodesWithLabel<'a> = NodesWithLabel<'a, Names>;

    fn nodes_with_label<'a>(&'a self, label: &'a str) -> Self::NodesWithLabel<'a> {
        NodesWithLabel {
            tax: self,
            iter: self
                .names
                .lookup_taxids(label)
                .map(Names::iter_lookup_taxids),
        }
    }
}

pub struct NodesWithLabel<'a, Names: NamesAssoc> {
    tax: &'a NcbiTaxonomy<Names>,
    iter: Option<Names::TaxIdsLookupIter<'a>>,
}

impl<'a, Names: NamesAssoc + 'static> Iterator for NodesWithLabel<'a, Names> {
    type Item = TaxId;

    fn next(&mut self) -> Option<Self::Item> {
        let next = self.iter.as_mut()?.next()?;

        self.tax.fixup_node(next.0)
    }
}

impl<Names: NamesAssoc + Send + 'static> TaxonomyMut for NcbiTaxonomy<Names> {
    type UnderlyingTopology = <GenericTaxonomy as TaxonomyMut>::UnderlyingTopology;

    fn replace_topology_with<Replacer>(&mut self, replacer: Replacer) -> Replacer::Result
    where
        Replacer: TopologyReplacer<Self::UnderlyingTopology>,
    {
        self.tree.replace_topology_with(replacer)
    }

    fn contract(&mut self, new_ranks: std::collections::HashSet<Self::RankSym>) -> ContractedNodes
    where
        Self::UnderlyingTopology:
            for<'b> Taxonomy<Children<'b> = Self::Children<'b>, RankSym = Self::RankSym>,
    {
        let contracted = self.replace_topology_with(taxonomy::Contractor {
            ranks_syms: new_ranks,
        });

        if let Some(old_contracted) = &mut self.contracted_taxids {
            old_contracted.0.extend(contracted.0.iter());
        } else {
            self.contracted_taxids = Some(contracted.clone());
        }

        contracted
    }
}

fn make_simple_tree(
    node_id: TaxId,
    children_lookup: &mut HashMap<TaxId, Vec<TaxId>>,
    names: &mut HashMap<TaxId, String>,
) -> SimpleTree {
    let children = children_lookup
        .get_mut(&node_id)
        .map(std::mem::take)
        .unwrap_or_default()
        .into_iter()
        .map(|child_id| make_simple_tree(child_id, children_lookup, names))
        .collect();

    SimpleTree {
        name: names.remove(&node_id).unwrap_or_default(),
        children,
        length: None,
    }
}

impl<Names> From<NcbiTaxonomy<Names>> for GenericTaxonomy {
    fn from(taxonomy: NcbiTaxonomy<Names>) -> Self {
        taxonomy.tree
    }
}

impl From<NcbiTaxonomy<SingleClassNames>> for SimpleTree {
    fn from(taxonomy: NcbiTaxonomy<SingleClassNames>) -> Self {
        let root = taxonomy.tree.root;
        let mut children_lookup = taxonomy.tree.children_lookup;
        let mut names = taxonomy.names.taxid_to_name;

        make_simple_tree(root, &mut children_lookup, &mut names)
    }
}
