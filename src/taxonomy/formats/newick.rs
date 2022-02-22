use std::{
    collections::HashMap,
    error::Error,
    io::{self, Read, Write},
    string::FromUtf8Error,
};

use newick_rs::{newick, SimpleTree};
use thiserror::Error;

use crate::taxonomy::NodeId;

pub struct NewickTaxonomy {
    pub root: NodeId,

    pub parent_ids: Vec<NodeId>,
    children_lookup: Vec<Vec<NodeId>>,

    pub labels: Vec<String>,
    label_lookup: HashMap<String, Vec<NodeId>>,
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

impl From<SimpleTree> for NewickTaxonomy {
    fn from(value: SimpleTree) -> Self {
        let mut label_lookup = HashMap::new();
        let mut labels = Vec::new();

        let mut new_node = |label: String| {
            let node_id = NodeId(labels.len());

            labels.push(label.clone());
            label_lookup
                .entry(label)
                .or_insert_with(|| Vec::new())
                .push(node_id);

            node_id
        };

        let root = new_node(value.name);

        let mut stack = vec![(root, Vec::new(), value.children.into_iter())];

        let mut parent_ids = Vec::new();
        let mut children_lookup = Vec::new();

        while let Some((node_id, children_ids, children)) = stack.last_mut() {
            if let Some(child) = children.next() {
                let child_id = new_node(child.name);
                parent_ids.push(*node_id);
                children_ids.push(child_id);
                stack.push((child_id, Vec::new(), child.children.into_iter()));
            } else {
                children_lookup.push(stack.pop().unwrap().1);
            }
        }

        NewickTaxonomy {
            root,
            parent_ids,
            children_lookup,
            labels,
            label_lookup,
        }
    }
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
