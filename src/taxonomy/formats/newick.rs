use std::{
    collections::HashMap,
    error::Error,
    io::{self, Read, Write},
    string::FromUtf8Error,
};

use newick_rs::{newick, SimpleTree};
use serde::{Serialize, Deserialize};
use thiserror::Error;

use crate::taxonomy::NodeId;

#[derive(Serialize, Deserialize)]
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
    ParseError(#[source] Box<dyn Error>),
}

impl NewickTaxonomy {
    pub fn from_simple_tree(value: SimpleTree, rank_storage: Vec<String>) -> Self {
        let mut label_lookup = HashMap::new();
        let mut labels = Vec::new();

        let mut depths = Vec::new();

        let mut new_node = |label: String, depth: usize| {
            let node_id = NodeId(labels.len());

            labels.push(label.clone());
            label_lookup
                .entry(label)
                .or_insert_with(|| Vec::new())
                .push(node_id);

            depths.push(depth);

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

        let mut parent_ids = Vec::new();
        let mut children_lookup = Vec::new();

        while let (depth, Some(frame)) = (stack.len(), stack.last_mut()) {
            if let Some(child) = frame.orig_children_iter.next() {
                let child_id = new_node(child.name, depth);
                parent_ids.push(frame.node_id);
                frame.children.push(child_id);
                stack.push(StackFrame::new(child_id, child.children));
            } else {
                children_lookup.push(stack.pop().unwrap().children);
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
}

pub fn load_newick_simple_tree<R: Read>(mut reader: R) -> Result<SimpleTree, NewickLoadError> {
    let contents = io::read_to_string(&mut reader)?;

    let tree: SimpleTree = newick::from_newick(&contents)
        .map_err(|err| NewickLoadError::ParseError(Box::new(err.to_owned())))?;

    Ok(tree.into())
}

pub fn save_newick_simple_tree<W: Write>(
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
