use core::fmt;
use std::collections::{HashMap, hash_map::Entry};

use serde::{Deserialize, Serialize, Deserializer, Serializer, ser::SerializeMap};
use string_interner::{backend::{StringBackend, Backend}, StringInterner, Symbol};
use thiserror::Error;

use super::{NodeId, Taxonomy, TaxonomyMut};


pub type RankSymbol = <StringBackend as Backend>::Symbol;

#[derive(Deserialize, Serialize)]
pub struct GenericTaxonomy {
    pub root: NodeId,

    pub parent_ids: HashMap<NodeId, NodeId>,
    pub children_lookup: HashMap<NodeId, Vec<NodeId>>,

    #[serde(deserialize_with = "deserialize_ranks", serialize_with = "serialize_ranks")]
    pub ranks: HashMap<NodeId, RankSymbol>,
    pub ranks_interner: StringInterner,
}

impl GenericTaxonomy {
    pub fn builder() -> GenericTaxonomyBuilder {
        GenericTaxonomyBuilder::default()
    }

    pub fn node_count(&self) -> usize {
        1 /* root */ + self.parent_ids.len()
    }

    pub fn has_node(&self, node: NodeId) -> bool {
        node == self.root || self.parent_ids.contains_key(&node)
    }
}

impl Taxonomy for GenericTaxonomy {
    fn get_root(&self) -> NodeId {
        self.root
    }

    fn is_leaf(&self, node: NodeId) -> NodeId {
        if let Some(ch) = self.children_lookup.get(&node) {
            ch.is_empty()
        } else {
            false
        }
    }

    type Children<'a> = std::iter::Copied<std::slice::Iter<'a, NodeId>>;

    fn iter_children<'a>(&'a self, node: NodeId) -> Self::Children<'a> {
        self.children_lookup
            .get(&node)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
            .iter()
            .copied()
    }

    fn find_parent(&self, node: NodeId) -> Option<NodeId> {
        if node == self.root {
            None
        } else {
            self.parent_ids.get(&node).copied()
        }
    }

    type RankSym = <StringBackend as Backend>::Symbol;

    fn lookup_rank_sym(&self, rank: &str) -> Option<Self::RankSym> {
        self.ranks_interner.get(rank)
    }

    fn find_rank(&self, node: NodeId) -> Option<Self::RankSym> {
        self.ranks.get(&node).copied()
    }
}

impl TaxonomyMut for GenericTaxonomy {
    type AddEdgeFn<'b> = &'b mut (dyn 'b + FnMut(NodeId, NodeId));

    fn replace_topology_with<'a, Body>(&'a mut self, body: Body)
    where
        Body: for<'b> FnOnce(&'a mut Self, Self::AddEdgeFn<'b>) {

        let mut parent_ids = HashMap::new();
        let mut children_lookup = HashMap::new();

        body(self, &mut |parent, child| {
            parent_ids.insert(child, parent);

            children_lookup
                .entry(parent)
                .or_insert_with(|| Vec::new())
                .push(child);
        });

        self.parent_ids = parent_ids;
        self.children_lookup = children_lookup;

        self.ranks.retain(|node, _| self.has_node(*node));
    }
}

pub fn serialize_ranks<S: Serializer>(ranks: &HashMap<NodeId, RankSymbol>, serializer: S) -> Result<S::Ok, S::Error> {
    let mut ranks_s = serializer.serialize_map(Some(ranks.len()))?;

    for (node, rank_sym) in ranks.iter() {
        ranks_s.serialize_entry(node, &rank_sym.to_usize())?;
    }

    ranks_s.end()
}

pub fn deserialize_ranks<'de, D: Deserializer<'de>>(
    deserializer: D,
) -> Result<HashMap<NodeId, RankSymbol>, D::Error> {
    use serde::de::{MapAccess, Visitor};

    struct RanksVisitor {}

    impl<'de> Visitor<'de> for RanksVisitor {
        type Value = HashMap<NodeId, RankSymbol>;

        fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
            formatter.write_str("ranks dictionary")
        }

        fn visit_map<M: MapAccess<'de>>(self, mut access: M) -> Result<Self::Value, M::Error> {
            let mut res = HashMap::with_capacity(access.size_hint().unwrap_or(0));

            while let Some((node, rank_sym_usize)) = access.next_entry()? {
                let rank_sym = RankSymbol::try_from_usize(rank_sym_usize).ok_or_else(|| {
                    <M::Error as serde::de::Error>::invalid_value(
                        serde::de::Unexpected::Unsigned(rank_sym_usize as u64),
                        &"a valid StringInterner symbol (as usize)",
                    )
                })?;

                res.insert(node, rank_sym);
            }

            Ok(res)
        }
    }

    deserializer.deserialize_map(RanksVisitor {})
}

#[derive(Default)]
pub struct GenericTaxonomyBuilder {
    root: Option<NodeId>,

    parent_ids: HashMap<NodeId, NodeId>,
    children_lookup: HashMap<NodeId, Vec<NodeId>>,

    ranks: HashMap<NodeId, RankSymbol>,
    ranks_interner: StringInterner,
}

#[derive(Debug, Error)]
pub enum TaxonomyBuildError {
    #[error("Root not found")]
    MissingRoot,

    #[error("Attempted to set multiple roots: {} and {}", .0, .1)]
    MultipleRoots(NodeId, NodeId),

    #[error("multiple parents found for {}: {} and {}", .0, .1, .2)]
    MultipleParents(NodeId, NodeId, NodeId),
}

impl GenericTaxonomyBuilder {
    pub fn set_root(&mut self, root: NodeId) -> Result<(), TaxonomyBuildError> {
        if let Some(r) = self.root {
            return Err(TaxonomyBuildError::MultipleRoots(r, root));
        } else {
            self.root = Some(root);
        }

        Ok(())
    }

    pub fn set_rank(&mut self, node: NodeId, rank: &str) {
        self.ranks.insert(node, self.ranks_interner.get_or_intern(rank));
    }

    pub fn insert_edge(&mut self, parent: NodeId, child: NodeId) -> Result<(), TaxonomyBuildError> {
        match self.parent_ids.entry(child) {
            Entry::Vacant(vac) => {
                vac.insert(parent);
            },
            Entry::Occupied(occ) => {
                return Err(TaxonomyBuildError::MultipleParents(child, parent, *occ.get()));
            },
        }

        self.children_lookup.entry(parent).or_insert_with(|| Vec::new()).push(child);

        Ok(())
    }

    pub fn build(self) -> Result<GenericTaxonomy, TaxonomyBuildError> {
        Ok(GenericTaxonomy {
            root: self.root.ok_or(TaxonomyBuildError::MissingRoot)?,
            parent_ids: self.parent_ids,
            children_lookup: self.children_lookup,
            ranks: self.ranks,
            ranks_interner: self.ranks_interner,
        })
    }
}
