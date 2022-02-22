#![feature(iterator_try_reduce)]
use std::{
    borrow::Borrow,
    collections::{HashMap, HashSet},
    fmt::{Debug, Display, Write},
    fs::File,
    hash::Hash,
    io::BufRead,
    io::BufReader,
    str::FromStr, cmp::Ordering,
};

use anyhow::{anyhow, bail, ensure, Result};
use clap::Parser;
use itertools::Itertools;
use ordered_float::NotNan;
use taxonomy::{GeneralTaxonomy, TaxRank, Taxonomy};
use thiserror::Error;

#[derive(Debug, Error)]
#[error("{}", .0)]
pub struct TaxonomyError(pub taxonomy::TaxonomyError);

impl From<taxonomy::TaxonomyError> for TaxonomyError {
    fn from(err: taxonomy::TaxonomyError) -> TaxonomyError {
        TaxonomyError(err)
    }
}

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub enum NodeType {
    Internal,
    Leaf,
}

pub struct TaxonomicTree<T, D> {
    pub root: usize,

    pub labels: Vec<T>,
    pub ranks: Vec<TaxRank>,

    pub parent_ids: Vec<usize>,
    pub parent_dists: Vec<D>,

    children_lookup: Vec<Vec<usize>>,
    label_lookup: HashMap<T, Vec<usize>>,
}

#[derive(Debug, Clone)]
pub struct TaxonomicTreeNode<'t, T, D> {
    pub parent: usize,
    pub label: &'t T,
    pub parent_distance: D,
    pub children: &'t [usize],
}

impl TaxonomicTree<String, f32> {
    pub fn from_general_taxonomy(
        tax: GeneralTaxonomy,
    ) -> Result<TaxonomicTree<String, f32>, TaxonomyError> {
        let mut tree = TaxonomicTree {
            children_lookup: (0..tax.tax_ids.len())
                .map(|id| tax.children(id))
                .try_collect()?,

            label_lookup: HashMap::with_capacity(tax.tax_ids.len()),

            root: tax.root(),
            labels: tax.tax_ids,

            ranks: tax.ranks,
            parent_ids: tax.parent_ids,
            parent_dists: tax.parent_dists,
        };

        tree.label_lookup = tree
            .labels
            .iter()
            .enumerate()
            .map(|(id, label)| (label.clone(), id))
            .into_group_map();

        Ok(tree)
    }
}

impl<T, D> TaxonomicTree<T, D> {
    pub fn lift(&mut self, new_root_label: T, lift_distance: D)
    where
        T: Eq + Hash + Clone,
    {
        let old_root = self.root;

        self.root = self.labels.len();

        self.labels.push(new_root_label.clone());
        self.ranks.push(TaxRank::Unspecified);

        self.parent_ids.push(self.root);
        self.parent_ids[old_root] = self.root;

        self.parent_dists.push(lift_distance);
        self.parent_dists.swap(old_root, self.root);

        self.children_lookup.push(vec![old_root]);
        self.label_lookup.insert(new_root_label, vec![self.root]);
    }

    pub fn is_root(&self, node_id: usize) -> bool {
        self.root == node_id
    }

    pub fn is_leaf(&self, node_id: usize) -> bool {
        self.children(node_id).is_empty()
    }

    pub fn node_type(&self, node_id: usize) -> NodeType {
        if self.is_leaf(node_id) {
            NodeType::Leaf
        } else {
            NodeType::Internal
        }
    }

    pub fn lookup_node_id<Q>(&self, label: &Q) -> Option<usize>
    where
        T: Eq + Hash + Borrow<Q>,
        Q: Eq + Hash + ?Sized,
    {
        self.label_lookup
            .get(label)
            .into_iter()
            .flatten()
            .exactly_one()
            .ok()
            .copied()
    }

    pub fn label(&self, node_id: usize) -> &T {
        &self.labels[node_id]
    }

    pub fn rank(&self, node_id: usize) -> TaxRank {
        self.ranks[node_id]
    }

    pub fn parent(&self, node_id: usize) -> usize {
        self.parent_ids[node_id]
    }

    pub fn parents(&self, node_id: usize) -> ParentsIter<'_, T, D> {
        ParentsIter::new(self, node_id)
    }

    pub fn lineage(&self, node_id: usize) -> impl Iterator<Item = (usize, TaxRank)> + '_ {
        [node_id]
            .into_iter()
            .chain(self.parents(node_id))
            .take_while(|&node_id| node_id != self.root)
            .map(|node_id| (node_id, self.rank(node_id)))
    }

    pub fn fill_lineage(&mut self, top_down_ranks: &[TaxRank])
    where
        T: Debug,
    {
        assert!(Some(top_down_ranks.len()) == self.check_equal_depths());

        let mut stack = Vec::new();

        stack.push((&self.children_lookup[self.root]).into_iter());

        let mut depth = 1;

        while let Some(iter) = stack.last_mut() {
            if let Some(&child_id) = iter.next() {
                self.ranks[child_id] = top_down_ranks[depth - 1];

                stack.push((&self.children_lookup[child_id]).into_iter());
                depth += 1;
            } else {
                stack.pop();
                depth -= 1;
            }
        }
    }

    pub fn depth(&self, node_id: usize) -> usize {
        self.parents(node_id).count()
    }

    pub fn lca(&self, node_id1: usize, node_id2: usize) -> Option<usize> {
        if self.is_leaf(node_id1) && self.is_leaf(node_id2) {
            self.leaf_lca(node_id1, node_id2)
        } else {
            let parents1 = self.parents(node_id1).collect_vec().into_iter().rev();
            let parents2 = self.parents(node_id2).collect_vec().into_iter().rev();

            itertools::izip!(parents1, parents2)
                .take_while(|(p1, p2)| p1 == p2)
                .last()
                .map(|(p, _)| p)
        }
    }

    pub fn leaf_lca(&self, node_id1: usize, node_id2: usize) -> Option<usize> {
        assert!(self.check_equal_depths().is_some());

        self.parents(node_id1)
            .zip(self.parents(node_id2))
            .skip_while(|(p1, p2)| p1 != p2)
            .next()
            .map(|(p, _)| p)
    }

    pub fn node(&self, node_id: usize) -> TaxonomicTreeNode<T, D>
    where
        D: Clone,
    {
        TaxonomicTreeNode {
            label: &self.labels[node_id],
            parent: self.parent_ids[node_id],
            parent_distance: self.parent_dists[node_id].clone(),
            children: self.children(node_id),
        }
    }

    pub fn children(&self, node_id: usize) -> &[usize] {
        &self.children_lookup[node_id]
    }

    pub fn subtree_height(&self, node_id: usize) -> usize {
        self.children(node_id)
            .iter()
            .map(|&child| self.subtree_height(child) + 1)
            .max()
            .unwrap_or(0)
    }

    pub fn preorder_descendants(&self, node_id: usize) -> PreOrderIter<T, D> {
        PreOrderIter::new(self, node_id)
    }

    pub fn postorder_descendants(&self, node_id: usize) -> PostOrderIter<T, D> {
        PostOrderIter::new(self, node_id)
    }

    pub fn descendant_leaves(&self, node_id: usize) -> impl Iterator<Item = usize> + '_ {
        self.preorder_descendants(node_id)
            .filter(|(_, ty)| *ty == NodeType::Leaf)
            .map(|(x, _)| x)
    }

    pub fn check_equal_depths(&self) -> Option<usize> {
        let mut known_depth = None;

        for depth in self
            .descendant_leaves(self.root)
            .map(|leaf_id| self.depth(leaf_id))
        {
            if let Some(d) = known_depth {
                if depth != d {
                    return None;
                }
            } else {
                known_depth = Some(depth);
            }
        }

        known_depth
    }
}

#[derive(Clone)]
pub struct ParentsIter<'t, T, D> {
    tree: &'t TaxonomicTree<T, D>,
    current: usize,
}

impl<'t, T, D> ParentsIter<'t, T, D> {
    pub fn new(tree: &'t TaxonomicTree<T, D>, node_id: usize) -> Self {
        Self {
            tree: tree,
            current: node_id,
        }
    }
}

impl<'t, T, D> Iterator for ParentsIter<'t, T, D> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.tree.is_root(self.current) {
            None
        } else {
            self.current = self.tree.parent(self.current);
            Some(self.current)
        }
    }
}

pub struct PreOrderIter<'t, T, D> {
    tree: &'t TaxonomicTree<T, D>,
    stack: Vec<std::slice::Iter<'t, usize>>,
}

impl<'t, T, D> PreOrderIter<'t, T, D> {
    pub fn new(tree: &'t TaxonomicTree<T, D>, node_id: usize) -> Self {
        Self {
            tree: tree,
            stack: vec![tree.children(node_id).iter()],
        }
    }

    pub fn without_node_types(self) -> impl Iterator<Item = usize> + 't {
        self.map(|(x, _)| x)
    }
}

impl<'t, T, D> Clone for PreOrderIter<'t, T, D> {
    fn clone(&self) -> Self {
        Self {
            tree: self.tree.clone(),
            stack: self.stack.clone(),
        }
    }
}

impl<'t, T, D> Iterator for PreOrderIter<'t, T, D> {
    type Item = (usize, NodeType);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(node) = self.stack.last_mut() {
            if let Some(&child) = node.next() {
                let grand_children = self.tree.children(child);

                if grand_children.is_empty() {
                    Some((child, NodeType::Leaf))
                } else {
                    self.stack.push(grand_children.iter());
                    Some((child, NodeType::Internal))
                }
            } else {
                self.stack.pop();
                self.next()
            }
        } else {
            None
        }
    }
}

pub struct PostOrderIter<'t, T, D> {
    tree: &'t TaxonomicTree<T, D>,
    stack: Vec<(usize, std::slice::Iter<'t, usize>)>,
}

impl<'t, T, D> PostOrderIter<'t, T, D> {
    pub fn new(tree: &'t TaxonomicTree<T, D>, node_id: usize) -> Self {
        Self {
            tree: tree,
            stack: vec![(node_id, tree.children(node_id).iter())],
        }
    }

    pub fn without_node_types(self) -> impl Iterator<Item = usize> + 't {
        self.map(|(x, _)| x)
    }
}

impl<'t, T, D> Clone for PostOrderIter<'t, T, D> {
    fn clone(&self) -> Self {
        Self {
            tree: self.tree.clone(),
            stack: self.stack.clone(),
        }
    }
}

impl<'t, T, D> Iterator for PostOrderIter<'t, T, D> {
    type Item = (usize, NodeType);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((current, children)) = self.stack.last_mut() {
            if let Some(&child) = children.next() {
                let grand_children = self.tree.children(child);
                self.stack.push((child, grand_children.iter()));
                self.next()
            } else {
                let current = *current;
                self.stack.pop();
                Some((current, self.tree.node_type(current)))
            }
        } else {
            None
        }
    }
}

#[derive(Default)]
pub struct TaxonAnnotations {
    pub descendants: usize,
    pub matches: usize,
    pub nonmatches: usize,
    pub tps: usize,
    pub fps: usize,
    pub tns: usize,
    pub fns: usize,
}

impl TaxonAnnotations {
    pub fn new() -> Self {
        Default::default()
    }
}

pub fn reads_to_nodes<S, SR, T, D, RI>(
    taxonomy: &TaxonomicTree<T, D>,
    reads: RI,
) -> Result<HashSet<usize>>
where
    S: Eq + Hash + Debug + ?Sized,
    SR: AsRef<S>,
    RI: IntoIterator<Item = SR>,
    T: Eq + Hash + Debug + Borrow<S>,
{
    reads
        .into_iter()
        .map(|read| {
            let id = taxonomy
                .lookup_node_id(read.as_ref())
                .ok_or_else(|| anyhow!("Read {:?} not found in taxonomy", read.as_ref()))?;

            ensure!(
                taxonomy.is_leaf(id),
                "Read doesn't map to a leaf: {:?}",
                read.as_ref()
            );

            Ok(id)
        })
        .try_collect()
}

pub fn annotate_match_count<T: Debug, D>(
    taxonomy: &TaxonomicTree<T, D>,
    reads_node_ids: &HashSet<usize>,
) -> Result<(usize, HashMap<usize, TaxonAnnotations>)> {
    let reads_lca_id: usize = reads_node_ids
        .iter()
        .copied()
        .try_reduce(|r1id, r2id| {
            taxonomy.lca(r1id, r2id).ok_or_else(|| {
                anyhow!(
                    "Error computing the LCA of {:?} and {:?}",
                    taxonomy.label(r1id),
                    taxonomy.label(r2id)
                )
            })
        })?
        .ok_or_else(|| anyhow!("No reads available to compute LCA"))?;

    let mut ann_map = HashMap::new();

    for (desc_id, node_type) in taxonomy.postorder_descendants(reads_lca_id) {
        let (descendants, matches) = match node_type {
            NodeType::Leaf => {
                let mut anns = TaxonAnnotations::new();
                if reads_node_ids.contains(&desc_id) {
                    anns.matches = 1;
                } else {
                    anns.matches = 0;
                }
                let matches = anns.matches;
                anns.nonmatches = 1 - matches;

                ann_map.insert(desc_id, anns);

                (0, matches)
            }
            NodeType::Internal => {
                let anns = ann_map
                    .get_mut(&desc_id)
                    .expect("Node annotations should have been created by its children");

                anns.nonmatches = anns.descendants - anns.matches;

                (anns.descendants, anns.matches)
            }
        };

        let parent_anns = ann_map
            .entry(taxonomy.parent(desc_id))
            .or_insert_with(TaxonAnnotations::new);

        parent_anns.descendants += descendants + 1;
        parent_anns.matches += matches;
    }

    Ok((reads_lca_id, ann_map))
}

pub fn annotate_precision_recall(
    ann_map: &mut HashMap<usize, TaxonAnnotations>,
    reads_lca_id: usize,
) {
    let (lca_m, lca_n) = {
        let &TaxonAnnotations {
            matches,
            nonmatches,
            ..
        } = ann_map
            .get(&reads_lca_id)
            .expect("LCA should already be annotated");

        (matches, nonmatches)
    };

    for ann in ann_map.values_mut() {
        let m = ann.matches;
        let n = ann.nonmatches;

        ann.tps = m;
        ann.fps = n;
        ann.tns = lca_n - n;
        ann.tps = lca_m - m;
    }
}

pub fn precision_recall_penalty(ann: &TaxonAnnotations, q: f64) -> f64 {
    if ann.tps != 0 {
        q * ann.fns as f64 / ann.tps as f64 + (1. - q) * ann.fps as f64 / ann.tps as f64
    } else {
        q * ann.fns as f64 + (1. - q) * ann.fps as f64
    }
}

pub fn assign_reads(ann_map: &HashMap<usize, TaxonAnnotations>, q: f64) -> (Vec<usize>, f64) {
    let mut min_nodes = Vec::new();
    let mut min_penalty = NotNan::new(f64::INFINITY).unwrap();

    for (node_id, ann) in ann_map {
        let penalty = NotNan::new(precision_recall_penalty(ann, q)).unwrap();

        match penalty.cmp(&min_penalty) {
            Ordering::Less => {
                min_nodes = vec![*node_id];
                min_penalty = penalty;
            },
            Ordering::Equal => {
                min_nodes.push(*node_id);
            },
            Ordering::Greater => {},
        }
    }

    (min_nodes, min_penalty.into_inner())
}

#[derive(Copy, Clone, Debug)]
pub enum TreeFormat {
    Newick,
}

#[derive(Debug, Error)]
pub enum TreeFormatParseError {
    #[error("Unrecognized tree format: {}", .0)]
    UnknownFormat(String),
}

impl FromStr for TreeFormat {
    type Err = TreeFormatParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Newick" => Ok(TreeFormat::Newick),
            _ => Err(TreeFormatParseError::UnknownFormat(s.to_owned())),
        }
    }
}

impl Display for TreeFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            TreeFormat::Newick => write!(f, "Newick"),
        }
    }
}

/// Rust port of TANGO: Taxonomic Assignment in Metagenomics
#[derive(Debug, Parser)]
#[clap(version)]
struct Args {
    /// Path to the taxonomic tree to be used
    #[clap(short, long)]
    tree: std::path::PathBuf,

    /// Format of the taxonomic tree
    #[clap(long, default_value_t = TreeFormat::Newick)]
    tree_format: TreeFormat,

    /// Text file containing the parsed output of a mapping program, in the format: [read_id]
    /// [species_id_1] ... [species_id_n]
    #[clap(short, long)]
    reads: std::path::PathBuf,

    /// Parameter that allows balancing the taxonomic assignment between precision (q = 0) and
    /// recall (q = 1), with q = 0.5 (default) corresponding to the F-measure (harmonic mean of
    /// precision and recall)
    #[clap(short, default_value_t = 0.5)]
    q: f64,
}

fn main() -> Result<()> {
    let args = Args::parse();

    ensure!(args.q >= 0. && args.q <= 1., "q must be between 0 and 1");

    let tax_tree = {
        let mut tree_reader = BufReader::new(File::open(&args.tree)?);

        let taxonomy = match args.tree_format {
            TreeFormat::Newick => {
                taxonomy::formats::newick::load_newick(&mut tree_reader).map_err(TaxonomyError)?
            }
        };

        let mut tax_tree = TaxonomicTree::from_general_taxonomy(taxonomy)?;

        let root_label = &mut tax_tree.labels[tax_tree.root];
        let mut root_label_borrow = root_label.trim();
        if root_label_borrow.ends_with(';') {
            root_label_borrow = root_label_borrow.strip_suffix(';').unwrap();
        }
        *root_label = root_label_borrow.to_owned();

        use TaxRank::*;
        tax_tree.lift("root".to_owned(), 0.);
        tax_tree.fill_lineage(&[Kingdom, Phylum, Class, Order, Family, Genus, Species]);
        tax_tree
    };

    for (line_num, line) in BufReader::new(File::open(&args.reads)?).lines().enumerate() {
        let line = line?;
        let mut reads_iter = line.trim().split_whitespace();
        let line_id = reads_iter
            .next()
            .ok_or_else(|| anyhow!("Invalid format at line {}: {:?}", line_num, line))?;

        let reads_node_ids = reads_to_nodes::<str, _, _, _, _>(&tax_tree, reads_iter)?;

        if reads_node_ids.is_empty() {
            continue;
        }

        let (reads_lca_id, mut anns) = annotate_match_count(&tax_tree, &reads_node_ids)?;
        annotate_precision_recall(&mut anns, reads_lca_id);

        let (assigned_nodes, penalty) = assign_reads(&anns, args.q);

        let mut output_line = line_id.to_owned();
        write!(output_line, " {}", penalty)?;

        for node_id in assigned_nodes {
            write!(
                output_line,
                " {} ({:?})",
                tax_tree.label(node_id),
                tax_tree.rank(node_id),
            )?;
            //}
            //None => {
                //bail!("BUG: Empty annotation set in assign_reads");
            //}
        }

        println!("{}", output_line);
    }

    Ok(())
}
