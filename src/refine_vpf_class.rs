use std::{collections::{HashSet, HashMap}, iter, path::Path, io, error};

use anyhow::{Result, Context};
use clap::Args;
use itertools::Itertools;
use serde::{Serialize, Deserialize, Serializer};

use crate::{
    preprocess_taxonomy::{PreprocessedTaxonomy, PreprocessedTaxonomyFormat, SomeTaxonomy, with_some_taxonomy},
    taxonomy::{NodeId, LabelledTaxonomy, Taxonomy}, assign::AssignmentRecord, util::writing_new_file_or_stdout,
};

#[derive(Args)]
pub struct RefineVpfClassArgs {
    #[clap(long, arg_enum, default_value_t)]
    taxonomy_format: PreprocessedTaxonomyFormat,

    /// Path to the preprocessed taxonomy (presumably from preprocess-taxonomy)
    preprocessed_taxonomy: String,

    /// TANGO3 metagenomic assignment output
    metagenomic_assignment: String,

    /// VPF-Class output file (presumably host_[rank].tsv)
    vpf_class_prediction: String,

    /// Does the presence of a higher-ranked node in the metagenomic assignment imply the presence
    /// of any of its descendants?
    #[clap(long)]
    include_descendants: bool,

    /// When mapping VPF classifications to taxids, ignore their rank.
    #[clap(long)]
    allow_different_ranks: bool,

    /// VPF-Class output's taxonomic rank. If not specified, it will be inferred from
    /// `VPF_CLASS_PREDICTION`'s file name
    #[clap(short, long)]
    rank: Option<String>,

    /// Refined VPF-Class output file (use '-' for STDOUT)
    #[clap(short, long)]
    output: String
}

impl RefineVpfClassArgs {
    fn get_or_infer_rank(&self) -> Option<String> {
        self.rank.to_owned().or_else(|| {
            let path: &Path = self.vpf_class_prediction.as_ref();

            path.file_name().and_then(|s| {
                s.to_string_lossy()
                    .strip_suffix(".tsv")
                    .and_then(|name| name.strip_prefix("host_"))
                    .map(ToOwned::to_owned)
            })
        })
    }
}

#[derive(Default)]
struct RankAssignments {
    rank_assignments: HashMap<NodeId, HashSet<NodeId>>,
    rank_contigs: HashMap<NodeId, HashSet<String>>,
    dropped_records: usize,
    total_records: usize,
}

impl RankAssignments {
    fn add_descendant_from_rank(&mut self, rank_node: NodeId, descendant: NodeId) {
        self.rank_assignments.entry(rank_node).or_default().insert(descendant);
    }

    fn add_ascendant_at_rank(&mut self, rank_node: NodeId, ascendant: NodeId) {
        self.rank_assignments.entry(rank_node).or_default().insert(ascendant);
    }

    fn add_contig_for_rank(&mut self, rank_node: NodeId, contig: String) {
        self.rank_contigs.entry(rank_node).or_default().insert(contig);
    }

    fn add_contig_for_rank_cloned(&mut self, rank_node: NodeId, contig: &str) {
        self.rank_contigs.entry(rank_node)
            .or_default()
            .get_or_insert_owned(contig);
    }

    fn lookup_rank_node(&self, rank_node: NodeId) -> Option<(&HashSet<NodeId>, &HashSet<String>)> {
        let taxids = self.rank_assignments.get(&rank_node);
        let contigs = self.rank_contigs.get(&rank_node);

        assert_eq!(taxids.is_some(), contigs.is_some());

        Option::zip(taxids, contigs)
    }
}

fn load_rank_assignments<Tax: Taxonomy, P: AsRef<Path>>(
    tax: &Tax,
    rank_sym: Tax::RankSym,
    include_descendants: bool,
    assignments_path: P,
) -> Result<RankAssignments> {
    let csv_reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(assignments_path)?;

    let mut assignments = RankAssignments::default();

    for record in csv_reader.into_deserialize() {
        let record: AssignmentRecord = record?;
        let node = NodeId(record.assigned_taxid);

        let valid_ancestor = iter::once(node)
            .chain(tax.ancestors(node))
            .find(|&ancestor| Some(rank_sym) == tax.find_rank(ancestor));

        if let Some(valid_ancestor) = valid_ancestor {
            assignments.add_descendant_from_rank(valid_ancestor, node);
            assignments.add_contig_for_rank(valid_ancestor, record.assigned_name);

        } else if include_descendants {
            let mut valid_descendants = tax
                .postorder_descendants(node)
                .filter(|&desc| Some(rank_sym) == tax.find_rank(desc))
                .peekable();

            if valid_descendants.peek().is_some() {
                for valid_descendant in valid_descendants {
                    assignments.add_ascendant_at_rank(valid_descendant, node);
                    assignments.add_contig_for_rank_cloned(valid_descendant, &record.assigned_name);
                }

            } else {
                assignments.dropped_records += 1;

                eprintln!(
                    "Warning: No valid descendants found for {:?} (taxid:{}) (rank: {})",
                    record.assigned_name,
                    node,
                    tax.find_rank(node).and_then(|sym| tax.rank_sym_str(sym)).unwrap_or(""),
                );
            }
        } else {
            assignments.dropped_records += 1;

            eprintln!(
                "Warning: No valid ancestors found for {:?} (taxid:{}) (rank: {})",
                record.assigned_name,
                node,
                tax.find_rank(node).and_then(|sym| tax.rank_sym_str(sym)).unwrap_or(""),
            );
        }

        assignments.total_records += 1;
    }

    Ok(assignments)
}

#[derive(Serialize, Deserialize)]
struct VpfClassRecord {
    virus_name: String,
    class_name: String,
    membership_ratio: f64,
    virus_hit_score: f64,
    confidence_score: f64,
}

#[derive(Serialize)]
struct EnrichedVpfClassRecord<'a> {
    #[serde(flatten)]
    vpf_class_record: VpfClassRecord,

    #[serde(serialize_with = "serialize_assigned_taxids")]
    assigned_taxids: HashSet<NodeId>,
    #[serde(serialize_with = "serialize_assigned_contigs")]
    assigned_contigs: HashSet<&'a str>,
}

impl<'a> EnrichedVpfClassRecord<'a> {
    fn new(vpf_class_record: VpfClassRecord) -> Self {
        EnrichedVpfClassRecord {
            vpf_class_record,
            assigned_taxids: HashSet::new(),
            assigned_contigs: HashSet::new(),
        }
    }

    fn has_assignments(&self) -> bool {
        !self.assigned_taxids.is_empty() || !self.assigned_contigs.is_empty()
    }
}

fn serialize_assigned_taxids<S: Serializer>(
    taxids: &HashSet<NodeId>,
    serializer: S,
) -> Result<S::Ok, S::Error> {
    serializer.serialize_str(&taxids.iter().map(|id| id.to_string()).join(";"))
}

fn serialize_assigned_contigs<S: Serializer>(
    contigs: &HashSet<&str>,
    serializer: S,
) -> Result<S::Ok, S::Error> {
    serializer.serialize_str(&contigs.iter().join(";"))
}

fn open_and_refine_vpf_class<'a, Tax: LabelledTaxonomy, P: AsRef<Path>>(
    tax: &'a Tax,
    rank_sym: Tax::RankSym,
    rank_assignments: &'a RankAssignments,
    allow_different_ranks: bool,
    path: P,
) -> Result<impl Iterator<Item = csv::Result<EnrichedVpfClassRecord<'a>>> + 'a> {
    let csv_reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(path)?;

    let result = csv_reader.into_deserialize().filter_map(move |record| {
        match record {
            Ok(record) => {
                let mut enriched = EnrichedVpfClassRecord::new(record);

                let mut found_mismatched_ranks = false;
                let mut found_taxids = false;

                for node in tax.nodes_with_label(&enriched.vpf_class_record.class_name) {
                    found_taxids = true;

                    if let Some((taxids, contigs)) = rank_assignments.lookup_rank_node(node) {
                        if allow_different_ranks || tax.find_rank(node) == Some(rank_sym) {
                            enriched.assigned_taxids.extend(taxids);
                            enriched.assigned_contigs.extend(contigs.iter().map(|s| s.as_str()));
                        } else {
                            found_mismatched_ranks = true;
                        }
                    }
                }

                if !found_taxids {
                    eprintln!("Warning: No taxids found for classification {}", &enriched.vpf_class_record.class_name);
                }

                if found_mismatched_ranks {
                    eprintln!(
                        "Warning: Found taxids for classification {} but had mismatching ranks, skipping (use --allow-different-ranks to disable the check)",
                        &enriched.vpf_class_record.class_name
                    );
                }

                if enriched.has_assignments() {
                    Some(Ok(enriched))
                } else {
                    None
                }
            }
            Err(err) => Some(Err(err)),
        }
    });

    Ok(result)
}

fn write_refined_output<'a, W, I, E: 'static>(writer: W, refinement: I) -> Result<()>
where
    W: io::Write,
    I: IntoIterator<Item = Result<EnrichedVpfClassRecord<'a>, E>>,
    E: error::Error + Send + Sync,
{
    let mut csv_writer = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(writer);

    for record in refinement {
        csv_writer.serialize(record?)?;
    }

    Ok(())
}

pub fn refine_vpf_class(args: RefineVpfClassArgs) -> Result<()> {
    let rank = args.get_or_infer_rank().ok_or_else(|| {
        anyhow::anyhow!(
            "--rank not specified and could not be deduced from {:?}",
            &args.vpf_class_prediction
        )
    })?;

    let taxonomy = PreprocessedTaxonomy::deserialize_with_format(
        &args.preprocessed_taxonomy,
        args.taxonomy_format
    )?;

    with_some_taxonomy!(&taxonomy.tree, tax => {
        let rank_sym = tax
            .lookup_rank_sym(&rank)
            .ok_or_else(|| anyhow::anyhow!("Rank not found in taxonomy: {rank}"))?;

        let rank_assignments = load_rank_assignments(
            tax,
            rank_sym,
            args.include_descendants,
            &args.metagenomic_assignment,
        )
        .context("Error collecting metagenomic assignment taxids")?;

        eprintln!(
            "Loaded assignments (dropped {} records out of {})",
            rank_assignments.dropped_records, rank_assignments.total_records
        );

        let refinement = open_and_refine_vpf_class(
            tax,
            rank_sym,
            &rank_assignments,
            args.allow_different_ranks,
            &args.vpf_class_prediction,
        )
        .context("Error reading VPF-Class output")?;

        writing_new_file_or_stdout!(&args.output, writer => {
            write_refined_output(writer, refinement)?;
        });
    });

    Ok(())
}
