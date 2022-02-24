use std::{fmt::Display, hash::Hash, collections::HashSet};

use serde::{Deserialize, Serialize};

#[derive(Debug, Copy, Clone, Eq, PartialEq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct NodeId(pub usize);

impl Display for NodeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

pub trait Taxonomy: Sized {
    fn get_root(&self) -> NodeId;

    fn is_leaf(&self, node: NodeId) -> bool {
        self.iter_children(node).next().is_none()
    }

    fn find_parent(&self, node: NodeId) -> Option<NodeId>;

    type Children<'a>: Iterator<Item = NodeId> + 'a
    where
        Self: 'a;

    fn iter_children<'a>(&'a self, node: NodeId) -> Self::Children<'a>;

    type RankSym: Eq + Copy + Hash;

    fn lookup_rank_sym(&self, rank: &str) -> Option<Self::RankSym>;

    fn find_rank(&self, node: NodeId) -> Option<Self::RankSym>;

    fn get_rank(&self, node: NodeId) -> Self::RankSym {
        self.find_rank(node)
            .unwrap_or_else(|| panic!("rank of {node} not found"))
    }

    fn ancestors<'a>(&'a self, node: NodeId) -> AncestorsIter<'a, Self> {
        AncestorsIter {
            taxo: self,
            current: node,
        }
    }

    fn preorder_descendants<'a>(&'a self, node: NodeId) -> PreOrderIter<'a, Self> {
        PreOrderIter::new(self, node)
    }

    fn postorder_descendants<'a>(&'a self, node: NodeId) -> PostOrderIter<'a, Self> {
        PostOrderIter::new(self, node)
    }
}

pub struct MissingRanks<R>(pub Vec<R>);

pub trait TaxonomyMut : Taxonomy {
    type AddEdgeFn<'b>: FnMut(NodeId, NodeId) + 'b;

    /// Need not free space
    fn replace_topology_with<Body>(&mut self, body: Body)
    where
        Body: for<'b> FnOnce(Self::AddEdgeFn<'b>);

    fn contract<'a, 'r, R: AsRef<str> + Clone>(&'a mut self, new_ranks: &'r [R]) -> MissingRanks<&'r R>
    where
        Self::Children<'a>: Clone,
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

        self.replace_topology_with(|add_edge| {
            struct StackFrame<'tax, Tax: Taxonomy + 'tax> {
                contraction_parent: NodeId,
                children_iter: Tax::Children<'tax>,
            }

            impl<'tax, Tax: Taxonomy> StackFrame<'tax, Tax> {
                fn new(contraction_parent: NodeId, children_iter: Tax::Children<'tax>) -> Self {
                    StackFrame { contraction_parent, children_iter }
                }
            }

            let root = self.get_root();
            let mut stack = vec![StackFrame::<Self>::new(root, self.iter_children(root))];

            while let Some(frame) = stack.last_mut() {
                if let Some(child) = frame.children_iter.next() {
                    let mut grandchildren = self.iter_children(child);

                    // leaf
                    if grandchildren.clone().next().is_none() {
                        add_edge(frame.contraction_parent, child);

                    // internal
                    } else {
                        match self.find_rank(child) {
                            Some(rank_sym) if ranks_syms.contains(&rank_sym) => {
                                add_edge(frame.contraction_parent, child);
                                stack.push(StackFrame::new(child, grandchildren));
                            }
                            _ => {
                                // frame.node is not a valid parent
                                stack.push(StackFrame::new(frame.contraction_parent, grandchildren));
                            }
                        }
                    }
                } else {
                    stack.pop();
                }
            }
        });

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
pub struct PreOrderIter<'tax, Tax: Taxonomy> {
    taxo: &'tax Tax,
    stack: Vec<Tax::Children<'tax>>,
}

impl<'tax, Tax: Taxonomy> PreOrderIter<'tax, Tax> {
    pub fn new(taxo: &'tax Tax, node: NodeId) -> Self {
        Self {
            taxo,
            stack: vec![taxo.iter_children(node)],
        }
    }
}

impl<'tax, Tax: Taxonomy> Iterator for PreOrderIter<'tax, Tax> {
    type Item = NodeId;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(node) = self.stack.last_mut() {
            if let Some(child) = node.next() {
                self.stack.push(self.taxo.iter_children(child));
                Some(child)
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
