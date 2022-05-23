use std::{fmt::{Debug, Display}, hash::Hash, collections::{HashSet, HashMap}, str::FromStr};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use self::util::Loop;

#[derive(Debug, Copy, Clone, Eq, PartialEq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct NodeId(pub usize);

impl Display for NodeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl FromStr for NodeId {
    type Err = <usize as FromStr>::Err;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(NodeId(s.parse()?))
    }
}

pub trait Taxonomy: Sized {
    fn get_root(&self) -> NodeId;

    fn is_leaf(&self, node: NodeId) -> bool {
        self.iter_children(node).next().is_none()
    }

    fn find_parent(&self, node: NodeId) -> Option<NodeId>;

    fn has_uniform_depths(&self) -> Option<usize>;

    type Children<'a>: Iterator<Item = NodeId> + 'a
    where
        Self: 'a;

    fn iter_children(&self, node: NodeId) -> Self::Children<'_>;

    type RankSym: Eq + Copy + Debug + Hash;

    fn rank_sym_str(&self, rank_sym: Self::RankSym) -> Option<&str>;

    fn lookup_rank_sym(&self, rank: &str) -> Option<Self::RankSym>;

    fn find_rank(&self, node: NodeId) -> Option<Self::RankSym>;

    fn find_rank_str(&self, node: NodeId) -> Option<&str> {
        self.find_rank(node).and_then(|sym| self.rank_sym_str(sym))
    }

    type NodeRanks<'a>: Iterator<Item = (NodeId, Self::RankSym)> + 'a
    where
        Self: 'a;

    fn node_ranks(&self) -> Self::NodeRanks<'_>;

    fn rank_ordering(&self) -> Result<Vec<Self::RankSym>, Loop<Self::RankSym>> {
        let rank_graph = self
            .preorder_edges(self.get_root())
            .filter_map(|(parent, child)| {
                let parent_rank = self.find_rank(parent)?;
                let child_rank = self.find_rank(child)?;

                if parent_rank != child_rank {
                    Some((parent_rank, child_rank))
                } else {
                    None
                }
            })
            .into_grouping_map()
            .collect();

        // turbofish because rustc goes crazy
        util::topsort::<Self::RankSym, HashSet<Self::RankSym>>(&rank_graph)
    }

    fn ancestors(&self, node: NodeId) -> AncestorsIter<'_, Self> {
        AncestorsIter {
            taxo: self,
            current: node,
        }
    }

    fn lca(&self, node1: NodeId, node2: NodeId) -> Option<NodeId> {
        if node1 == node2 {
            Some(node1)

        } else if self.has_uniform_depths().is_some() && self.is_leaf(node1) && self.is_leaf(node2) {
            self.ancestors(node1)
                .zip(self.ancestors(node2))
                .find(|(p1, p2)| p1 == p2)
                .map(|(p, _)| p)

        } else {
            let parents1 = self.ancestors(node1).collect_vec().into_iter().rev().chain([node1]);
            let parents2 = self.ancestors(node2).collect_vec().into_iter().rev().chain([node2]);

            itertools::izip!(parents1, parents2)
                .take_while(|(p1, p2)| p1 == p2)
                .last()
                .map(|(p, _)| p)
        }
    }

    fn preorder_edges(&self, node: NodeId) -> PreOrderEdgesIter<'_, Self> {
        PreOrderEdgesIter::new(self, node)
    }

    fn postorder_descendants(&self, node: NodeId) -> PostOrderIter<'_, Self> {
        PostOrderIter::new(self, node)
    }
}

pub trait LabelledTaxonomy : Taxonomy {
    type Labels<'a>: Iterator<Item = &'a str>
        where Self: 'a;

    fn labels_of(&self, node: NodeId) -> Self::Labels<'_>;

    fn some_label_of(&self, node: NodeId) -> Option<&str> {
        self.labels_of(node).next()
    }

    type NodesWithLabel<'a>: Iterator<Item = NodeId> + 'a
        where Self: 'a;

    fn nodes_with_label<'a>(&'a self, label: &'a str) -> Self::NodesWithLabel<'a>;

    fn some_node_with_label(&self, label: &str) -> Option<NodeId> {
        self.nodes_with_label(label).next()
    }
}

pub struct MissingRanks<R>(pub Vec<R>);

pub trait TopologyReplacer<Tax> {
    type Result = ();

    fn replace_topology_with<AddEdge>(self, tax: &Tax, add_edge: AddEdge) -> Self::Result
    where
        AddEdge: FnMut(NodeId, NodeId);
}

pub struct Contractor<Tax: Taxonomy> {
    ranks_syms: HashSet<Tax::RankSym>,
}

impl<Tax: Taxonomy> Contractor<Tax> {
    pub fn new(ranks_syms: HashSet<Tax::RankSym>) -> Self {
        Self { ranks_syms }
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct ContractedNodes(pub HashMap<NodeId, NodeId>);

impl<Tax: Taxonomy> TopologyReplacer<Tax> for Contractor<Tax> {
    type Result = ContractedNodes;

    fn replace_topology_with<AddEdge>(self, tax: &Tax, mut add_edge: AddEdge) -> Self::Result
    where
        AddEdge: FnMut(NodeId, NodeId),
    {
        struct StackFrame<'tax, Tax: Taxonomy + 'tax> {
            contraction_parent: NodeId,
            children_iter: Tax::Children<'tax>,
        }

        impl<'tax, Tax: Taxonomy> StackFrame<'tax, Tax> {
            fn root_new(tax: &'tax Tax) -> Self {
                let root = tax.get_root();
                Self::new(root, tax.iter_children(root))
            }

            fn new(contraction_parent: NodeId, children_iter: Tax::Children<'tax>) -> Self {
                StackFrame { contraction_parent, children_iter }
            }
        }

        let mut dropped = HashMap::new();
        let mut stack = vec![StackFrame::root_new(tax)];

        while let Some(frame) = stack.last_mut() {
            if let Some(child) = frame.children_iter.next() {
                let child_valid_rank = tax
                    .find_rank(child)
                    .map_or(false, |rank_sym| self.ranks_syms.contains(&rank_sym));

                let contraction_parent = frame.contraction_parent;

                let grandchildren = tax.iter_children(child);

                if child_valid_rank {
                    add_edge(frame.contraction_parent, child);
                    stack.push(StackFrame::new(child, grandchildren));
                } else {
                    dropped.insert(child, contraction_parent);
                    stack.push(StackFrame::new(contraction_parent, grandchildren));
                }
            } else {
                stack.pop();
            }
        }

        ContractedNodes(dropped)
    }
}

pub trait TaxonomyMut: Taxonomy {
    type UnderlyingTopology: Taxonomy = Self;

    /// Need not free space
    fn replace_topology_with<Replacer>(&mut self, replacer: Replacer) -> Replacer::Result
    where
        Replacer: TopologyReplacer<Self::UnderlyingTopology>;

    fn contract(&mut self, new_ranks: HashSet<Self::RankSym>) -> ContractedNodes
    where
        Self::UnderlyingTopology:
            for<'b> Taxonomy<Children<'b> = Self::Children<'b>, RankSym = Self::RankSym>,
    {
        self.replace_topology_with(Contractor { ranks_syms: new_ranks })
    }
}

#[derive(Clone)]
pub struct AncestorsIter<'tax, Tax> {
    taxo: &'tax Tax,
    current: NodeId,
}

impl<'tax, Tax: Taxonomy> Iterator for AncestorsIter<'tax, Tax> {
    type Item = NodeId;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(parent) = self.taxo.find_parent(self.current) {
            self.current = parent;
            Some(parent)
        } else {
            // root
            None
        }
    }
}

#[derive(Clone)]
pub struct PreOrderEdgesIter<'tax, Tax: Taxonomy> {
    taxo: &'tax Tax,
    stack: Vec<(NodeId, Tax::Children<'tax>)>,
}

impl<'tax, Tax: Taxonomy> PreOrderEdgesIter<'tax, Tax> {
    pub fn new(taxo: &'tax Tax, node: NodeId) -> Self {
        Self {
            taxo,
            stack: vec![(node, taxo.iter_children(node))],
        }
    }
}

impl<'tax, Tax: Taxonomy> Iterator for PreOrderEdgesIter<'tax, Tax> {
    type Item = (NodeId, NodeId);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((current, children)) = self.stack.last_mut() {
            if let Some(child) = children.next() {
                let current = *current;
                self.stack.push((child, self.taxo.iter_children(child)));
                Some((current, child))
            } else {
                self.stack.pop();
                self.next()
            }
        } else {
            None
        }
    }
}

#[derive(Clone)]
pub struct PostOrderIter<'tax, Tax: Taxonomy> {
    taxo: &'tax Tax,
    stack: Vec<(NodeId, Tax::Children<'tax>)>,
}

impl<'tax, Tax: Taxonomy> PostOrderIter<'tax, Tax> {
    pub fn new(taxo: &'tax Tax, node_id: NodeId) -> Self {
        Self {
            taxo,
            stack: vec![(node_id, taxo.iter_children(node_id))],
        }
    }
}

impl<'tax, Tax: Taxonomy> Iterator for PostOrderIter<'tax, Tax> {
    type Item = NodeId;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((current, children)) = self.stack.last_mut() {
            if let Some(child) = children.next() {
                let grand_children = self.taxo.iter_children(child);
                self.stack.push((child, grand_children));
                self.next()
            } else {
                let current = *current;
                self.stack.pop();
                Some(current)
            }
        } else {
            None
        }
    }
}


pub mod generic;
pub mod formats;
pub(crate) mod util;
