use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{self, Read, Write},
    iter,
    path::Path,
    slice,
    string::FromUtf8Error,
};

use itertools::Itertools;
use newick_rs::{newick, SimpleTree};
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::taxonomy::{LabelledTaxonomy, NodeId, Taxonomy};

#[derive(Debug, Serialize, Deserialize)]
pub struct NewickTaxonomy {
    pub root: NodeId,

    pub parent_ids: Vec<NodeId>,
    children_lookup: Vec<Vec<NodeId>>,

    pub labels: Vec<String>,
    label_lookup: HashMap<String, Vec<NodeId>>,

    depths: Vec<usize>,
    rank_storage: Vec<String>,
}

#[derive(Debug, Error)]
pub enum NewickLoadError {
    #[error("IO error loading newick")]
    IoError(#[from] io::Error),

    #[error("Encoding error of newick data")]
    EncodingError(#[from] FromUtf8Error),

    #[error("Newick parse error: {}", .0)]
    ParseError(#[source] Box<dyn Error + Send + Sync>),
}

impl NewickTaxonomy {
    pub fn from_simple_tree(value: SimpleTree, ranks: Vec<String>) -> Self {
        let rank_storage = ranks;

        let mut label_lookup = HashMap::new();
        let mut labels = Vec::new();

        let mut depths = Vec::new();

        let mut new_node = |label: String, depth: usize| {
            let node_id = NodeId(labels.len());

            labels.push(label.clone());

            label_lookup
                .entry(label)
                .or_insert_with(Vec::new)
                .push(node_id);

            depths.push(depth);

            assert!(depth < rank_storage.len());

            node_id
        };

        struct StackFrame {
            node_id: NodeId,
            children: Vec<NodeId>,
            orig_children_iter: std::vec::IntoIter<SimpleTree>,
        }

        impl StackFrame {
            fn new(node_id: NodeId, orig_children: Vec<SimpleTree>) -> Self {
                StackFrame {
                    node_id,
                    children: Vec::new(),
                    orig_children_iter: orig_children.into_iter(),
                }
            }
        }

        let root = new_node(value.name, 0);
        let mut stack = vec![StackFrame::new(root, value.children)];

        let mut parent_ids = vec![root];
        let mut children_lookup = vec![vec![]];

        while let (depth, Some(frame)) = (stack.len(), stack.last_mut()) {
            if let Some(child) = frame.orig_children_iter.next() {
                let child_id = new_node(child.name, depth);

                assert!(parent_ids.len() == child_id.0);
                parent_ids.push(frame.node_id);
                frame.children.push(child_id);

                assert!(children_lookup.len() == child_id.0);
                children_lookup.push(vec![]);

                stack.push(StackFrame::new(child_id, child.children));
            } else {
                let popped = stack.pop().unwrap();
                children_lookup[popped.node_id.0] = popped.children;
            }
        }

        NewickTaxonomy {
            root,
            parent_ids,
            children_lookup,
            labels,
            label_lookup,
            depths,
            rank_storage,
        }
    }

    pub fn load_newick(
        path: impl AsRef<Path>,
        ranks: Vec<String>,
    ) -> Result<Self, NewickLoadError> {
        let simple_tree = read_newick_simple_tree(File::open(path)?)?;

        Ok(Self::from_simple_tree(simple_tree, ranks))
    }

    pub fn all_nodes(&self) -> impl Iterator<Item = NodeId> {
        (0..self.parent_ids.len()).map(NodeId)
    }
}

impl Taxonomy for NewickTaxonomy {
    fn get_root(&self) -> NodeId {
        self.root
    }

    fn find_parent(&self, node: NodeId) -> Option<NodeId> {
        if node == self.root {
            None
        } else {
            self.parent_ids.get(node.0).copied()
        }
    }

    fn has_uniform_depths(&self) -> Option<usize> {
        let mut depths = self
            .all_nodes()
            .filter(|&node| self.is_leaf(node))
            .map(|node| self.depths[node.0])
            .peekable();

        let &depth = depths.peek().expect("Tree has no leaves");

        if depths.all_equal() {
            Some(depth)
        } else {
            None
        }
    }

    type Children<'a> = std::iter::Copied<std::slice::Iter<'a, NodeId>>;

    fn iter_children(&self, node: NodeId) -> Self::Children<'_> {
        self.children_lookup[node.0].iter().copied()
    }

    type RankSym = usize;

    fn rank_sym_str(&self, rank_sym: Self::RankSym) -> Option<&str> {
        self.rank_storage.get(rank_sym).map(|s| &**s)
    }

    fn lookup_rank_sym(&self, rank: &str) -> Option<Self::RankSym> {
        self.rank_storage.iter().position(|r| r == rank)
    }

    fn find_rank(&self, node: NodeId) -> Option<Self::RankSym> {
        self.depths.get(node.0).copied()
    }

    type NodeRanks<'a> = EnumerateAsNodeId<iter::Copied<slice::Iter<'a, Self::RankSym>>>;

    fn node_ranks(&self) -> Self::NodeRanks<'_> {
        EnumerateAsNodeId {
            inner: self.depths.iter().copied().enumerate(),
        }
    }
}

impl LabelledTaxonomy for NewickTaxonomy {
    type Labels<'a> = std::option::IntoIter<&'a str>;

    fn labels_of(&self, node: NodeId) -> Self::Labels<'_> {
        self.labels.get(node.0).map(String::as_str).into_iter()
    }

    type NodesWithLabel<'a> =
        std::iter::Copied<std::iter::Flatten<std::option::IntoIter<&'a Vec<NodeId>>>>;

    fn nodes_with_label<'a>(&'a self, label: &'a str) -> Self::NodesWithLabel<'a> {
        self.label_lookup.get(label).into_iter().flatten().copied()
    }
}

pub struct EnumerateAsNodeId<I> {
    inner: iter::Enumerate<I>,
}

impl<I: Iterator> Iterator for EnumerateAsNodeId<I> {
    type Item = (NodeId, <I as Iterator>::Item);

    fn next(&mut self) -> Option<Self::Item> {
        let (i, x) = self.inner.next()?;

        Some((NodeId(i), x))
    }
}

pub fn read_newick_simple_tree<R: Read>(mut reader: R) -> Result<SimpleTree, NewickLoadError> {
    let contents = io::read_to_string(&mut reader)?;

    let tree: SimpleTree = newick::from_newick(&contents)
        .map_err(|err| NewickLoadError::ParseError(Box::new(err.to_owned())))?;

    Ok(tree)
}

pub fn write_newick_simple_tree<W: Write>(
    tree: &SimpleTree,
    mut writer: W,
) -> Result<(), io::Error> {
    writer.write_all(newick::to_newick(tree).as_bytes())
}

fn make_simple_tree(
    node_id: NodeId,
    children_lookup: &mut Vec<Vec<NodeId>>,
    labels: &mut Vec<String>,
) -> SimpleTree {
    let children = std::mem::take(&mut children_lookup[node_id.0])
        .into_iter()
        .map(|child_id| make_simple_tree(child_id, children_lookup, labels))
        .collect();

    SimpleTree {
        name: std::mem::take(&mut labels[node_id.0]),
        children,
        length: None,
    }
}

impl From<NewickTaxonomy> for SimpleTree {
    fn from(taxonomy: NewickTaxonomy) -> Self {
        let root = taxonomy.root;
        let mut children_lookup = taxonomy.children_lookup;
        let mut labels = taxonomy.labels;

        make_simple_tree(root, &mut children_lookup, &mut labels)
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use newick_rs::SimpleTree;

    use crate::taxonomy::{formats::newick::NewickTaxonomy, LabelledTaxonomy, Taxonomy};

    macro_rules! newick_literal {
        ([ $($trees:tt),+ ] $name:expr) => {
            SimpleTree {
                name: $name.to_string(),
                children: vec![ $( newick_literal!($trees) ),+ ],
                length: None
            }
        };

        (($sub:ident)) => {
            $sub
        };

        ($name:expr) => {
            SimpleTree {
                name: $name.to_string(),
                children: vec![],
                length: None
            }
        };
    }

    #[test]
    fn ancestors() {
        fn chain(n: u32) -> SimpleTree {
            if n == 0 {
                newick_literal!(0u32)
            } else {
                let sub = chain(n - 1);
                newick_literal!( [ (sub) ] n )
            }
        }

        let simple_t = chain(5);

        let ranks = "a,b,c,d,e,f".split(',').map(|s| s.to_owned()).collect_vec();

        let t = NewickTaxonomy::from_simple_tree(simple_t, ranks);

        assert_eq!(
            t.ancestors(t.nodes_with_label("0").next().unwrap())
                .map(|node| t.labels_of(node).next().unwrap())
                .collect_vec(),
            vec!["1", "2", "3", "4", "5"]
        );
    }

    #[test]
    fn postorder_descendants() {
        let sub = newick_literal! { [1i32, 2i32] 3i32 };
        let simple_t = newick_literal! { [ 0i32, (sub) ] 4i32 };

        let ranks = "a,b,c".split(',').map(|s| s.to_owned()).collect_vec();

        let t = NewickTaxonomy::from_simple_tree(simple_t, ranks);

        println!("{t:?}");

        assert_eq!(
            t.postorder_descendants(t.get_root())
                .map(|node| t.some_label_of(node).unwrap())
                .collect_vec(),
            vec!["0", "1", "2", "3", "4"]
        );
    }

    #[test]
    fn preorder_edges() {
        let sub = newick_literal! { [1i32, 2i32] 3i32 };
        let simple_t = newick_literal! { [ 0i32, (sub) ] 4i32 };

        let ranks = "a,b,c".split(',').map(|s| s.to_owned()).collect_vec();

        let t = NewickTaxonomy::from_simple_tree(simple_t, ranks);

        assert_eq!(
            t.preorder_edges(t.get_root())
                .map(|(_parent, node)| t.some_label_of(node).unwrap())
                .collect_vec(),
            vec!["0", "3", "1", "2"]
        );
    }

    #[test]
    fn lca() {
        let sub = newick_literal! { [1i32, 2i32] 3i32 };
        let simple_t = newick_literal! { [ 0i32, (sub) ] 4i32 };

        // [0, [1, 2] 3] 4

        let ranks = "a,b,c".split(',').map(|s| s.to_owned()).collect_vec();

        let t = NewickTaxonomy::from_simple_tree(simple_t, ranks);

        assert_eq!(
            t.lca(
                t.some_node_with_label("1").unwrap(),
                t.some_node_with_label("2").unwrap(),
            ),
            Some(t.some_node_with_label("3").unwrap()),
        );

        assert_eq!(
            t.lca(
                t.some_node_with_label("0").unwrap(),
                t.some_node_with_label("1").unwrap(),
            ),
            Some(t.some_node_with_label("4").unwrap()),
        );

        assert_eq!(
            t.lca(
                t.some_node_with_label("1").unwrap(),
                t.some_node_with_label("3").unwrap(),
            ),
            Some(t.some_node_with_label("3").unwrap()),
        );

        assert_eq!(
            t.lca(
                t.some_node_with_label("1").unwrap(),
                t.some_node_with_label("1").unwrap(),
            ),
            Some(t.some_node_with_label("1").unwrap()),
        );
    }
}
