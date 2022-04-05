use std::{fmt::Display, hash::Hash, collections::HashSet, iter, str::FromStr};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

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

    fn iter_children<'a>(&'a self, node: NodeId) -> Self::Children<'a>;

    type RankSym: Eq + Copy + Hash;

    fn rank_sym_str(&self, rank_sym: Self::RankSym) -> Option<&str>;

    fn lookup_rank_sym(&self, rank: &str) -> Option<Self::RankSym>;

    fn find_rank(&self, node: NodeId) -> Option<Self::RankSym>;

    fn get_rank(&self, node: NodeId) -> Self::RankSym {
        self.find_rank(node)
            .unwrap_or_else(|| panic!("rank of {node} not found"))
    }

    type NodeRanks<'a>: Iterator<Item = (NodeId, Self::RankSym)> + 'a
    where
        Self: 'a;

    fn node_ranks<'a>(&'a self) -> Self::NodeRanks<'a>;

    fn ancestors<'a>(&'a self, node: NodeId) -> AncestorsIter<'a, Self> {
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
                .skip_while(|(p1, p2)| p1 != p2)
                .next()
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

    fn preorder_edges<'a>(&'a self, node: NodeId) -> PreOrderEdgesIter<'a, Self> {
        PreOrderEdgesIter::new(self, node)
    }

    fn postorder_descendants<'a>(&'a self, node: NodeId) -> PostOrderIter<'a, Self> {
        PostOrderIter::new(self, node)
    }
}

pub trait LabelledTaxonomy : Taxonomy {
    type Labels<'a>: Iterator<Item = &'a str>
        where Self: 'a;

    fn labels_of<'a>(&'a self, node: NodeId) -> Self::Labels<'a>;

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
    fn replace_topology_with<AddEdge>(self, tax: &Tax, add_edge: AddEdge)
    where
        AddEdge: FnMut(NodeId, NodeId);
}

struct Contractor<Tax: Taxonomy> {
    ranks_syms: HashSet<Tax::RankSym>,
}

impl<Tax: Taxonomy> TopologyReplacer<Tax> for Contractor<Tax> {
    fn replace_topology_with<AddEdge>(self, tax: &Tax, mut add_edge: AddEdge)
    where
        AddEdge: FnMut(NodeId, NodeId),
    {
        struct StackFrame<'tax, Tax: Taxonomy + 'tax> {
            contraction_parent: NodeId,
            children_iter: iter::Peekable<Tax::Children<'tax>>,
        }

        impl<'tax, Tax: Taxonomy> StackFrame<'tax, Tax> {
            fn root_new(tax: &'tax Tax) -> Self {
                let root = tax.get_root();
                Self::new(root, tax.iter_children(root).peekable())
            }

            fn new(contraction_parent: NodeId, children_iter: iter::Peekable<Tax::Children<'tax>>) -> Self {
                StackFrame { contraction_parent, children_iter }
            }
        }

        let mut stack = vec![StackFrame::root_new(tax)];

        while let Some(frame) = stack.last_mut() {
            if let Some(child) = frame.children_iter.next() {
                let mut grandchildren = tax.iter_children(child).peekable();

                // leaf
                if grandchildren.peek().is_none() {
                    add_edge(frame.contraction_parent, child);

                // internal
                } else {
                    match tax.find_rank(child) {
                        Some(rank_sym) if self.ranks_syms.contains(&rank_sym) => {
                            add_edge(frame.contraction_parent, child);
                            stack.push(StackFrame::new(child, grandchildren));
                        }
                        _ => {
                            // frame.node is not a valid parent
                            let contraction_parent = frame.contraction_parent;
                            stack.push(StackFrame::new(contraction_parent, grandchildren));
                        }
                    }
                }
            } else {
                stack.pop();
            }
        }
    }
}

pub trait TaxonomyMut: Taxonomy {
    type UnderlyingTopology: Taxonomy = Self;

    /// Need not free space
    fn replace_topology_with<Replacer>(&mut self, replacer: Replacer)
    where
        Replacer: TopologyReplacer<Self::UnderlyingTopology>;

    fn contract<'a, 'r, R>(&'a mut self, new_ranks: &'r [R]) -> MissingRanks<&'r R>
    where
        R: AsRef<str> + Clone,
        Self::UnderlyingTopology: for<'b> Taxonomy<Children<'b> = Self::Children<'b>, RankSym = Self::RankSym>,
    {
        let mut missing_ranks = MissingRanks(Vec::new());
        let mut ranks_syms = HashSet::new();

        for rank in new_ranks {
            if let Some(rank_sym) = self.lookup_rank_sym(rank.as_ref()) {
                ranks_syms.insert(rank_sym);
            } else {
                missing_ranks.0.push(rank);
            }
        }

        self.replace_topology_with(Contractor { ranks_syms });

        missing_ranks
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
