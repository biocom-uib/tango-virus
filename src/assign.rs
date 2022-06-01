use std::collections::{HashMap, HashSet};
use std::io;

use clap::Args;
use crossbeam::channel::{Receiver, Sender};
use csv::{StringRecord, WriterBuilder};
use extended_rational::Rational;
use itertools::Itertools;
use rayon::iter::{ParallelBridge, ParallelIterator};
use serde::{Deserialize, Serialize};

use crate::preprocess_blastout;
use crate::preprocess_taxonomy::{
    with_some_ncbi_or_newick_taxonomy, with_some_taxonomy, PreprocessedTaxonomy,
    PreprocessedTaxonomyFormat, SomeTaxonomy,
};
use crate::taxonomy::formats::ncbi::{NamesAssoc, NcbiTaxonomy};
use crate::taxonomy::{LabelledTaxonomy, NodeId, Taxonomy};
use crate::util::writing_new_file_or_stdout;

#[derive(Clone, Debug, Default)]
struct TaxonAnnotations {
    descendant_leaves: usize,
    matches: usize,
    nonmatches: usize,
    tps: usize,
    fps: usize,
    //tns: usize,
    fns: usize,
}

impl TaxonAnnotations {
    pub fn new() -> Self {
        Default::default()
    }
}

type TaxonomyAnnotations = HashMap<NodeId, TaxonAnnotations>;

fn annotate_match_count<Tax: LabelledTaxonomy>(
    taxonomy: &Tax,
    reads_node_ids: &HashSet<NodeId>,
) -> anyhow::Result<(NodeId, TaxonomyAnnotations)> {
    let reads_lca_id: NodeId = reads_node_ids
        .iter()
        .copied()
        .try_reduce(|r1id, r2id| {
            taxonomy.lca(r1id, r2id).ok_or_else(|| {
                let ancestors1 = taxonomy.ancestors(r1id).map(|n| n.0).collect_vec();
                let ancestors2 = taxonomy.ancestors(r2id).map(|n| n.0).collect_vec();
                anyhow::anyhow!(
                    "Error computing the LCA of {:?} ({r1id}) and {:?} ({r2id}) with ancestors\n{ancestors1:?}\nand\n{ancestors2:?}",
                    taxonomy.some_label_of(r1id).unwrap_or(""),
                    taxonomy.some_label_of(r2id).unwrap_or(""),
                )
            })
        })?
        .ok_or_else(|| anyhow::anyhow!("No reads available to compute LCA"))?;

    let mut ann_map = HashMap::new();

    for node in taxonomy.postorder_descendants(reads_lca_id) {
        let (descendant_leaves, matches) = if taxonomy.is_leaf(node) {
            let mut anns = TaxonAnnotations::new();

            anns.descendant_leaves = 1;

            if reads_node_ids.contains(&node) {
                anns.matches = 1;
            } else {
                anns.matches = 0;
            }

            let matches = anns.matches;
            anns.nonmatches = 1 - matches;

            assert!(ann_map.insert(node, anns).is_none());

            (1, matches)
        } else {
            let anns = ann_map
                .get_mut(&node)
                .expect("Node annotations should have been created by its children");

            anns.nonmatches = anns.descendant_leaves - anns.matches;

            (anns.descendant_leaves, anns.matches)
        };

        if node != taxonomy.get_root() && node != reads_lca_id {
            let parent_anns = ann_map
                .entry(taxonomy.find_parent(node).unwrap())
                .or_insert_with(TaxonAnnotations::new);

            parent_anns.descendant_leaves += descendant_leaves;
            parent_anns.matches += matches;
        }
    }

    Ok((reads_lca_id, ann_map))
}

fn annotate_precision_recall(ann_map: &mut TaxonomyAnnotations, reads_lca_id: NodeId) {
    let lca_ann = ann_map
        .get(&reads_lca_id)
        .expect("LCA should already be annotated")
        .clone();

    for ann in ann_map.values_mut() {
        ann.tps = ann.matches;
        ann.fps = ann.nonmatches;
        //ann.tns = lca_ann.nonmatches - ann.nonmatches;
        ann.fns = lca_ann.matches - ann.matches;
    }
}

fn precision_recall_penalty(ann: &TaxonAnnotations, q: Rational) -> Rational {
    let r = |num, denom| Rational::new(num as i64, denom as i64);

    if ann.tps != 0 {
        q * r(ann.fns, ann.tps) + (r(1, 1) - q) * r(ann.fps, ann.tps)
    } else {
        q * r(ann.fns, 1) + (r(1, 1) - q) * r(ann.fps, 1)
    }
}

fn assign_reads(ann_map: &TaxonomyAnnotations, q: Rational) -> (Vec<NodeId>, Rational) {
    let mut min_nodes = Vec::new();
    let mut min_penalty = Rational::new(-1, 1);

    for (node_id, ann) in ann_map {
        let penalty = precision_recall_penalty(ann, q);

        assert!(!penalty.is_infinity());

        if min_penalty.is_negative() || penalty < min_penalty {
            min_penalty = penalty;
            min_nodes = vec![*node_id];
        } else if penalty == min_penalty {
            min_nodes.push(*node_id);
        }
    }

    (min_nodes, min_penalty)
}

#[derive(Args)]
pub struct AssignArgs {
    #[clap(long, arg_enum, default_value_t)]
    taxonomy_format: PreprocessedTaxonomyFormat,

    /// Path to the preprocessed taxonomy (presumably from preprocess-taxonomy)
    preprocessed_taxonomy: String,

    /// Path to the preprocessed reads file (presumably from preprocess-blastout). It should be a
    /// file containing the parsed output of a mapping program, in the (tab-separated) format
    ///
    /// [read_id]    [species_id_1];...;[species_id_n]
    preprocessed_reads: String,

    /// Output path (use '-' for STDOUT)
    #[clap(short, long, default_value = "-")]
    output: String,

    /// Parameter that allows balancing the taxonomic assignment between precision (q = 0) and
    /// recall (q = 1), with q = 0.5 (default) corresponding to the F-measure (harmonic mean of
    /// precision and recall)
    #[clap(short, default_value_t = 0.5)]
    q: f64,
}

struct ReadsCsvHeader {
    //query_id_col: String,
    subjects_id_col: String,
    weights_col: Option<String>,
}

fn peek_reads_csv_header<R>(csv_reader: &mut csv::Reader<R>) -> anyhow::Result<ReadsCsvHeader>
where
    R: io::Read,
{
    let mut record = StringRecord::new();
    csv_reader.read_record(&mut record)?;
    assert!(record.len() == 2);
    let _query_id_col = record[0].to_owned();
    let subjects_col = &record[1];

    if let Some((subjects_id_col, weights_id)) = subjects_col.split_once('/') {
        Ok(ReadsCsvHeader {
            //query_id_col,
            subjects_id_col: subjects_id_col.to_owned(),
            weights_col: Some(weights_id.to_owned()),
        })
    } else {
        Ok(ReadsCsvHeader {
            //query_id_col,
            subjects_id_col: subjects_col.to_owned(),
            weights_col: None,
        })
    }
}

struct QueryAssignments {
    query_id: String,
    assigned_nodes: Vec<NodeId>,
    penalty: f32,
}

fn produce_assignments<R, Tax, F>(
    mut csv_reader: csv::Reader<R>,
    header: &ReadsCsvHeader,
    tax: &Tax,
    read_to_taxids: F,
    q: Rational,
    sender: Sender<QueryAssignments>,
) -> anyhow::Result<()>
where
    R: io::Read + Send,
    Tax: LabelledTaxonomy + Sync,
    F: Fn(&str) -> Vec<NodeId> + Sync,
{
    csv_reader.records().enumerate().par_bridge().try_for_each(|(record_id, record)| {
        let record = record?;
        assert!(record.len() == 2);

        let query_id = &record[0];
        let subject_ids = &record[1];

        if subject_ids.is_empty() {
            return Ok(());
        }

        let subject_ids: Vec<&str> = if header.weights_col.is_some() {
            let linenum = record_id + 1;

            record[1].split(';').map(|item| Some(item.split_once('/')?.0)).collect::<Option<_>>()
                .ok_or_else(|| anyhow::anyhow!("Error parsing reads at line {linenum}"))?
        } else {
            record[1].split(';').collect()
        };

        let subject_node_ids: HashSet<_> = subject_ids.iter().filter_map(|&name| {
            let taxids = read_to_taxids(name);

            match taxids.len() {
                0 => {
                    eprintln!("\nWarning: No taxid matching {name:?} found in the taxonomy, skipping");
                    None
                },
                1 => {
                    Some(taxids[0])
                },
                _ => {
                    eprintln!("\nWarning: Multiple taxids matched {name:?} in the taxonomy, returning the latest one");
                    taxids.iter().max().copied()
                },
            }
        })
        .collect();

        if subject_node_ids.is_empty() {
            eprintln!("\nWarning: Skipping {query_id} as none out of {} subject ids were recognized", subject_ids.len());

        } else {
            let (lca_id, mut ann_map) = annotate_match_count(tax, &subject_node_ids)?;

            annotate_precision_recall(&mut ann_map, lca_id);
            let (assigned_nodes, penalty) = assign_reads(&ann_map, q);

            sender.send(QueryAssignments {
                query_id: query_id.to_owned(),
                assigned_nodes,
                penalty: f32::from(penalty),
            })?;
        }

        anyhow::Ok(())
    })
}

fn ncbi_taxonomy_lookup_taxid_leaf<Names: 'static + NamesAssoc + Send>(
    tax: &NcbiTaxonomy<Names>,
) -> impl Fn(&str) -> Vec<NodeId> + '_ {
    |name: &str| {
        let mut v = name
            .parse()
            .into_iter()
            .map(NodeId)
            .filter_map(|node| tax.fixup_node(node))
            .filter(|&node| tax.is_leaf(node))
            .collect_vec();

        v.sort();
        v.dedup();
        v
    }
}

fn labelled_taxonomy_lookup_taxid_leaf<Tax>(tax: &Tax) -> impl Fn(&str) -> Vec<NodeId> + '_
where
    Tax: LabelledTaxonomy,
{
    |name: &str| {
        let mut v = tax
            .nodes_with_label(name)
            .filter(|&node| tax.is_leaf(node))
            .collect_vec();

        v.sort();
        v.dedup();
        v
    }
}

fn load_reads_and_produce_assignments(
    reads_path: &str,
    taxonomy: &PreprocessedTaxonomy,
    q: Rational,
    sender: Sender<QueryAssignments>,
) -> anyhow::Result<()> {
    let mut csv_reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(reads_path)?;

    let header = peek_reads_csv_header(&mut csv_reader)?;

    macro_rules! produce_assignments_with {
        ($tax:expr, $lookup_taxid:expr) => {
            produce_assignments(csv_reader, &header, $tax, $lookup_taxid, q, sender)
        };
    }

    with_some_ncbi_or_newick_taxonomy!(&taxonomy.tree,
        ncbi: tax => {
            if header.subjects_id_col == preprocess_blastout::fields::STAXID.0 {
                produce_assignments_with!(tax, ncbi_taxonomy_lookup_taxid_leaf(tax))?;
            } else if header.subjects_id_col == preprocess_blastout::fields::SSCINAME.0 {
                produce_assignments_with!(tax, labelled_taxonomy_lookup_taxid_leaf(tax))?;
            } else {
                anyhow::bail!("Unable to map read subject ID's ({}) to NCBI Taxonomy ID's", &header.subjects_id_col)
            }
        },
        newick: tax => {
            eprintln!("Warning: Unable to ensure that Newick taxonomy labels match the search subject ID key ({})", &header.subjects_id_col);

            produce_assignments_with!(tax, labelled_taxonomy_lookup_taxid_leaf(tax))?;
        }
    );

    Ok(())
}

#[derive(Serialize, Deserialize)]
pub struct AssignmentRecord {
    pub query_id: String,
    pub assigned_taxid: usize,
    pub assigned_name: String,
    pub assigned_rank: String,
    pub penalty: f32,
}

fn write_assignments(
    output: &str,
    taxonomy: &PreprocessedTaxonomy,
    assignments: Receiver<QueryAssignments>,
) -> anyhow::Result<()> {
    with_some_taxonomy!(&taxonomy.tree, tax => {
        writing_new_file_or_stdout!(output, writer => {
            let mut csv_writer = WriterBuilder::new()
                .delimiter(b'\t')
                .has_headers(true)
                .from_writer(writer);

            eprintln!("");

            let mut count = 0;

            for record in assignments {
                for node in record.assigned_nodes {
                    let assigned_name = tax.labels_of(node).next().unwrap_or_else(|| {
                        eprintln!("\nWarning: No label found for taxid {node}");
                        ""
                    });

                    let rank = tax
                        .find_rank(node)
                        .and_then(|rank_sym| tax.rank_sym_str(rank_sym))
                        .unwrap_or("")
                        .to_owned();

                    let penalty = record.penalty;

                    csv_writer.serialize(AssignmentRecord {
                        query_id: record.query_id.clone(),
                        assigned_taxid: node.0,
                        assigned_name: assigned_name.to_owned(),
                        assigned_rank: rank,
                        penalty,
                    })?;
                }

                count += 1;
                eprint!("\rProcessed {count} query sequences");
            }
        });
    });

    Ok(())
}

pub fn assign(args: AssignArgs) -> anyhow::Result<()> {
    let taxonomy = PreprocessedTaxonomy::deserialize_with_format(
        &args.preprocessed_taxonomy,
        args.taxonomy_format,
    )?;
    eprintln!("Loaded taxonomy");

    let taxonomy = &taxonomy;
    let q = Rational::from(args.q);

    let (sender, receiver) = crossbeam::channel::unbounded();

    let (rl, rr) = rayon::join(
        move || load_reads_and_produce_assignments(&args.preprocessed_reads, taxonomy, q, sender),
        move || write_assignments(&args.output, taxonomy, receiver),
    );

    rl?;
    rr?;

    Ok(())
}
