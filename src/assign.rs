use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io;

use clap::Args;
use crossbeam::channel::Sender;
use csv::{StringRecord, WriterBuilder};
use extended_rational::Rational;
use itertools::Itertools;
use rayon::iter::{ParallelBridge, ParallelIterator};
use serde::{Serialize, Deserialize};

use crate::preprocess_blastout;
use crate::preprocess_taxonomy::{PreprocessedTaxonomy, PreprocessedTaxonomyFormat, SomeTaxonomy};
use crate::taxonomy::formats::ncbi::NcbiTaxonomy;
use crate::taxonomy::{LabelledTaxonomy, NodeId};


#[derive(Clone, Debug, Default)]
struct TaxonAnnotations {
    descendants: usize,
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

    for desc_id in taxonomy.postorder_descendants(reads_lca_id) {
        let (descendants, matches) = if taxonomy.is_leaf(desc_id) {
            let mut anns = TaxonAnnotations::new();
            if reads_node_ids.contains(&desc_id) {
                anns.matches = 1;
            } else {
                anns.matches = 0;
            }
            let matches = anns.matches;
            anns.nonmatches = 1 - matches;

            assert!(ann_map.insert(desc_id, anns).is_none());

            (0, matches)
        } else {
            let anns = ann_map
                .get_mut(&desc_id)
                .expect("Node annotations should have been created by its children");

            anns.nonmatches = 1 + anns.descendants - anns.matches;

            (anns.descendants, anns.matches)
        };

        if desc_id != taxonomy.get_root() && desc_id != reads_lca_id {
            let parent_anns = ann_map
                .entry(taxonomy.find_parent(desc_id).unwrap())
                .or_insert_with(TaxonAnnotations::new);

            parent_anns.descendants += 1 + descendants;
            parent_anns.matches += matches;
        }
    }

    Ok((reads_lca_id, ann_map))
}

fn annotate_precision_recall(ann_map: &mut TaxonomyAnnotations, reads_lca_id: NodeId) {
    let lca_ann = ann_map.get(&reads_lca_id)
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
    #[clap(long, arg_enum, default_value_t = PreprocessedTaxonomyFormat::CBOR)]
    taxonomy_format: PreprocessedTaxonomyFormat,

    /// Path to the preprocessed taxonomy (presumably from preprocess-taxonomy)
    preprocessed_taxonomy: String,

    /// Path to the preprocessed reads file (presumably from preprocess-blastout). It should be a
    /// file containing the parsed output of a mapping program, in the format
    ///
    /// [read_id]	[species_id_1];...;[species_id_n]
    ///          (tab)
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
    query_id_col: String,
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
    let query_id_col = record[0].to_owned();
    let subjects_col = &record[1];

    if let Some((subjects_id_col, weights_id)) = subjects_col.split_once('/') {
        Ok(ReadsCsvHeader {
            query_id_col,
            subjects_id_col: subjects_id_col.to_owned(),
            weights_col: Some(weights_id.to_owned())
        })
    } else {
        Ok(ReadsCsvHeader {
            query_id_col,
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

        let subject_node_ids: HashSet<NodeId> = subject_ids.iter()
            .filter_map(|&name| {
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

        if subject_node_ids.len() < 2 { return Ok(()); }

        let (lca_id, mut ann_map) = annotate_match_count(tax, &subject_node_ids)?;

        annotate_precision_recall(&mut ann_map, lca_id);
        let (assigned_nodes, penalty) = assign_reads(&ann_map, q);

        sender.send(QueryAssignments {
            query_id: query_id.to_owned(),
            assigned_nodes,
            penalty: f32::from(penalty),
        })?;

        anyhow::Ok(())
    })
}

fn ncbi_taxonomy_lookup_taxid<'a, Names>(
    tax: &'a NcbiTaxonomy<Names>,
) -> impl Fn(&str) -> Vec<NodeId> + 'a
{
    |name: &str| {
        name.parse()
            .into_iter()
            .map(NodeId)
            .filter_map(|node| tax.fixup_node(node))
            .collect::<HashSet<_>>()
            .into_iter()
            .collect_vec()
    }
}

fn labelled_taxonomy_lookup_taxid<'a, Tax>(tax: &'a Tax) -> impl Fn(&str) -> Vec<NodeId> + 'a
where
    Tax: LabelledTaxonomy,
{
    |name: &str| tax.nodes_with_label(name).collect_vec()
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

    match &taxonomy.tree {
        SomeTaxonomy::NcbiTaxonomyWithSingleClassNames(tax) => {
            if header.subjects_id_col == preprocess_blastout::fields::STAXID.0 {
                produce_assignments(
                    csv_reader,
                    &header,
                    tax,
                    ncbi_taxonomy_lookup_taxid(tax),
                    q,
                    sender,
                )?;
            } else if header.subjects_id_col == preprocess_blastout::fields::SSCINAME.0 {
                produce_assignments(
                    csv_reader,
                    &header,
                    tax,
                    labelled_taxonomy_lookup_taxid(tax),
                    q,
                    sender,
                )?;
            } else {
                anyhow::bail!("Unable to map read subject ID's to NCBI Taxonomy ID's")
            }
        }
        SomeTaxonomy::NcbiTaxonomyWithManyNames(tax) => {
            if header.subjects_id_col == preprocess_blastout::fields::STAXID.0 {
                produce_assignments(
                    csv_reader,
                    &header,
                    tax,
                    ncbi_taxonomy_lookup_taxid(tax),
                    q,
                    sender,
                )?;
            } else if header.subjects_id_col == preprocess_blastout::fields::SSCINAME.0 {
                produce_assignments(
                    csv_reader,
                    &header,
                    tax,
                    labelled_taxonomy_lookup_taxid(tax),
                    q,
                    sender,
                )?;
            } else {
                anyhow::bail!("Unable to map read subject ID's to NCBI Taxonomy ID's")
            }
        },
        SomeTaxonomy::NewickTaxonomy(tax) => {
            eprintln!("Warning: Unable to verify Newick taxonomy labels, check that subject ID's are correctly matched");

            produce_assignments(
                csv_reader,
                &header,
                tax,
                labelled_taxonomy_lookup_taxid(tax),
                q,
                sender,
            )?;
        }
    }

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

fn write_assignments<W, Tax, I>(writer: W, tax: &Tax, records: I) -> anyhow::Result<()>
where
    W: io::Write,
    Tax: LabelledTaxonomy,
    I: IntoIterator<Item = QueryAssignments>,
{
    let mut csv_writer = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(writer);

    eprintln!("");

    let mut count = 0;

    for record in records {
        for assigned_node in record.assigned_nodes {
            let query_id = record.query_id.clone();

            let assigned_taxid = assigned_node.0;

            let assigned_name = match tax.labels_of(assigned_node).exactly_one() {
                Ok(label) => {
                    if label.is_empty() {
                        eprintln!("\nWarning: Empty label found for taxid {assigned_node}");
                    }
                    label
                }
                Err(mut err) => {
                    if let Some(label) = err.next() {
                        eprintln!("\nWarning: Found more than one label for taxid {assigned_node}");
                        label
                    } else {
                        eprintln!("\nWarning: No label found for taxid {assigned_node}");
                        ""
                    }
                }
            };
            let assigned_name = assigned_name.to_owned();

            let assigned_rank = tax
                .find_rank(assigned_node)
                .and_then(|rank_sym| tax.rank_sym_str(rank_sym))
                .unwrap_or("")
                .to_owned();

            let penalty = record.penalty;

            csv_writer.serialize(AssignmentRecord {
                query_id,
                assigned_taxid,
                assigned_name,
                assigned_rank,
                penalty,
            })?;
        }

        count += 1;
        eprint!("\rProcessed {count} query sequences");
    }

    Ok(())
}

fn open_output_and_write_assignments<I>(
    output: &str,
    taxonomy: &PreprocessedTaxonomy,
    records: I,
) -> anyhow::Result<()>
where
    I: IntoIterator<Item = QueryAssignments>,
{
    match &taxonomy.tree {
        SomeTaxonomy::NcbiTaxonomyWithSingleClassNames(tax) => {
            if output == "-" {
                write_assignments(io::stdout(), tax, records)?;
            } else {
                write_assignments(File::create(output)?, tax, records)?;
            };
        }
        SomeTaxonomy::NcbiTaxonomyWithManyNames(tax) => {
            if output == "-" {
                write_assignments(io::stdout(), tax, records)?;
            } else {
                write_assignments(File::create(output)?, tax, records)?;
            };
        }
        SomeTaxonomy::NewickTaxonomy(tax) => {
            if output == "-" {
                write_assignments(io::stdout(), tax, records)?;
            } else {
                write_assignments(File::create(output)?, tax, records)?;
            };
        }
    }

    Ok(())
}

pub fn assign(args: AssignArgs) -> anyhow::Result<()> {
    let taxonomy = PreprocessedTaxonomy::deserialize_with_format(
        &args.preprocessed_taxonomy,
        args.taxonomy_format,
    )?;
    let taxonomy = &taxonomy;
    let q = Rational::from(args.q);

    let (sender, receiver) = crossbeam::channel::unbounded();

    let (rl, rr) = rayon::join(
        move || load_reads_and_produce_assignments(&args.preprocessed_reads, taxonomy, q, sender),
        move || open_output_and_write_assignments(&args.output, taxonomy, receiver),
    );

    rl?;
    rr?;

    Ok(())
}
