use std::{
    collections::{HashMap, HashSet},
    io, iter,
    num::ParseFloatError,
    path::Path,
};

use anyhow::{Context, Result};
use clap::{Args, ArgEnum};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::{
    tango_assign::AssignmentRecord,
    filter::{self, FromStrFilter},
    preprocess_taxonomy::{
        with_some_taxonomy, PreprocessedTaxonomy, PreprocessedTaxonomyFormat, SomeTaxonomy,
    },
    taxonomy::{LabelledTaxonomy, NodeId, Taxonomy},
    util::writing_new_file_or_stdout,
};

/// Refine a VPF-Class host prediction to include only ascendants from taxonomic nodes found in a
/// metagenomic assignment.
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
    output: String,

    /// Filters to apply before processing. Example: --filter 'membership_ratio>=0.2'
    #[clap(long)]
    filter: Vec<String>,

    /// Print a summary of class frequencies and their average membership ratio (incompatible with
    /// '--output -').
    #[clap(long)]
    print_summary: bool,

    /// Column to sort --print-summary by
    #[clap(long, arg_enum, default_value_t = SummarySortBy::VirusCount)]
    summary_sort_by: SummarySortBy,
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
        self.rank_assignments
            .entry(rank_node)
            .or_default()
            .insert(descendant);
    }

    fn add_ascendant_at_rank(&mut self, rank_node: NodeId, ascendant: NodeId) {
        self.rank_assignments
            .entry(rank_node)
            .or_default()
            .insert(ascendant);
    }

    fn add_contig_for_rank(&mut self, rank_node: NodeId, contig: String) {
        self.rank_contigs
            .entry(rank_node)
            .or_default()
            .insert(contig);
    }

    fn add_contig_for_rank_cloned(&mut self, rank_node: NodeId, contig: &str) {
        self.rank_contigs
            .entry(rank_node)
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

fn load_rank_assignments<Tax: Taxonomy>(
    tax: &Tax,
    rank_sym: Tax::RankSym,
    include_descendants: bool,
    assignments_path: &Path,
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
            assignments.add_contig_for_rank(valid_ancestor, record.query_id);

        } else if include_descendants {
            let mut valid_descendants = tax
                .postorder_descendants(node)
                .filter(|&desc| Some(rank_sym) == tax.find_rank(desc))
                .peekable();

            if valid_descendants.peek().is_some() {
                for valid_descendant in valid_descendants {
                    assignments.add_ascendant_at_rank(valid_descendant, node);
                    assignments.add_contig_for_rank_cloned(valid_descendant, &record.query_id);
                }
            } else {
                assignments.dropped_records += 1;

                eprintln!(
                    "Warning: No valid descendants found for {:?} (taxid:{}) (rank: {})",
                    record.assigned_name,
                    node,
                    tax.find_rank(node)
                        .and_then(|sym| tax.rank_sym_str(sym))
                        .unwrap_or(""),
                );
            }
        } else {
            assignments.dropped_records += 1;

            eprintln!(
                "Warning: No valid ancestors found for {:?} (taxid:{}) (rank: {})",
                record.assigned_name,
                node,
                tax.find_rank(node)
                    .and_then(|sym| tax.rank_sym_str(sym))
                    .unwrap_or(""),
            );
        }

        assignments.total_records += 1;
    }

    Ok(assignments)
}

#[derive(Debug, Serialize, Deserialize)]
struct VpfClassRecord {
    pub virus_name: String,
    pub class_name: String,
    pub membership_ratio: f64,
    pub virus_hit_score: f64,
    pub confidence_score: f64,
}

enum VpfClassRecordFieldValue {
    VirusName(String),
    ClassName(String),
    MembershipRatio(f64),
    VirusHitScore(f64),
    ConfidenceScore(f64),
}

struct VpfClassRecordFilter {
    field_value: VpfClassRecordFieldValue,
    op: filter::Op,
}

impl VpfClassRecordFilter {
    fn apply(&self, record: &VpfClassRecord) -> bool {
        use VpfClassRecordFieldValue::*;

        match &self.field_value {
            VirusName(value) => self.op.apply(&record.virus_name, value),
            ClassName(value) => self.op.apply(&record.class_name, value),
            MembershipRatio(value) => self.op.apply(&record.membership_ratio, value),
            VirusHitScore(value) => self.op.apply(&record.virus_hit_score, value),
            ConfidenceScore(value) => self.op.apply(&record.confidence_score, value),
        }
    }
}

#[derive(Debug, Error)]
enum VpfClassRecordFilterParseError {
    #[error("Error parsing floating point number")]
    FloatParseError(#[from] ParseFloatError),

    #[error("Unknown field {:?}, Valid names are virus_name, class_name, membership_ratio, virus_hit_score and confidence_score", .0)]
    UnknownField(String),
}

impl FromStrFilter for VpfClassRecordFilter {
    type Err = VpfClassRecordFilterParseError;

    fn try_from_parts(key: &str, op: filter::Op, value: &str) -> Result<Self, Self::Err> {
        use VpfClassRecordFieldValue::*;
        use VpfClassRecordFilterParseError::*;

        let field_value = match key {
            "virus_name" => VirusName(value.to_owned()),
            "class_name" => ClassName(value.to_owned()),
            "membership_ratio" => MembershipRatio(value.parse()?),
            "virus_hit_score" => VirusHitScore(value.parse()?),
            "confidence_score" => ConfidenceScore(value.parse()?),
            _ => return Err(UnknownField(key.to_owned())),
        };

        Ok(VpfClassRecordFilter { field_value, op })
    }
}

struct EnrichedVpfClassRecord<'a> {
    vpf_class_record: VpfClassRecord,
    assigned_taxids: HashSet<NodeId>,
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

#[derive(Debug, Serialize)]
struct CsvEnrichedVpfClassRecord {
    virus_name: String,
    class_name: String,
    membership_ratio: f64,
    virus_hit_score: f64,
    confidence_score: f64,
    assigned_taxids: String,
    num_assigned_contigs: usize,
}

impl<'a> From<EnrichedVpfClassRecord<'a>> for CsvEnrichedVpfClassRecord {
    fn from(other: EnrichedVpfClassRecord<'a>) -> Self {
        CsvEnrichedVpfClassRecord {
            virus_name: other.vpf_class_record.virus_name,
            class_name: other.vpf_class_record.class_name,
            membership_ratio: other.vpf_class_record.membership_ratio,
            virus_hit_score: other.vpf_class_record.virus_hit_score,
            confidence_score: other.vpf_class_record.confidence_score,
            assigned_taxids: other.assigned_taxids.iter().join(";"),
            num_assigned_contigs: other.assigned_contigs.len(),
        }
    }
}


fn enrich_vpf_class_record<'a, Tax: LabelledTaxonomy>(
    tax: &'a Tax,
    rank_sym: Tax::RankSym,
    rank_assignments: &'a RankAssignments,
    allow_different_ranks: bool,
    record: VpfClassRecord,
) -> Option<EnrichedVpfClassRecord<'a>> {
    let mut enriched = EnrichedVpfClassRecord::new(record);

    let mut found_mismatched_ranks = false;
    let mut found_taxids = false;

    for node in tax.nodes_with_label(&enriched.vpf_class_record.class_name) {
        found_taxids = true;

        if let Some((taxids, contigs)) = rank_assignments.lookup_rank_node(node) {
            if allow_different_ranks || tax.find_rank(node) == Some(rank_sym) {
                enriched.assigned_taxids.extend(taxids);
                enriched
                    .assigned_contigs
                    .extend(contigs.iter().map(|s| s.as_str()));
            } else {
                found_mismatched_ranks = true;
            }
        }
    }

    if !found_taxids {
        eprintln!(
            "Warning: No taxids found for classification {}",
            &enriched.vpf_class_record.class_name
        );
    }

    if found_mismatched_ranks {
        eprintln!(
            "Warning: Found taxids for classification {} but had mismatching ranks, skipping (use --allow-different-ranks to disable the check)",
            &enriched.vpf_class_record.class_name
        );
    }

    if enriched.has_assignments() {
        Some(enriched)
    } else {
        None
    }
}

#[derive(Default, Serialize)]
struct ClassStats<'a> {
    virus_count: u32,
    assigned_taxids: HashSet<NodeId>,
    assigned_contigs: HashSet<&'a str>
}

impl<'a> ClassStats<'a> {
    fn add(&mut self, record: &EnrichedVpfClassRecord<'a>) {
        self.virus_count += 1;
        self.assigned_taxids.extend(&record.assigned_taxids);
        self.assigned_contigs.extend(&record.assigned_contigs);
    }
}

#[derive(ArgEnum, Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum SummarySortBy {
    ClassName,
    VirusCount,
    AssignedTaxIds,
    AssignedContigs,
}

#[derive(Default)]
struct EnrichmentSummary<'a> {
    class_stats: HashMap<String, ClassStats<'a>>,
}

impl<'a> EnrichmentSummary<'a> {
    fn account(&mut self, record: &EnrichedVpfClassRecord<'a>) {
        let class_name = &*record.vpf_class_record.class_name;

        self
            .class_stats
            .raw_entry_mut()
            .from_key(class_name)
            .or_insert_with(|| (class_name.to_owned(), ClassStats::default()))
            .1
            .add(record);
    }

    fn write<W: io::Write>(self, writer: W, sort: SummarySortBy) -> csv::Result<()> {
        #[derive(Serialize)]
        struct Record {
            class_name: String,
            virus_count: i32,
            assigned_taxids: i32,
            assigned_contigs: i32,
        }

        let mut records = self
            .class_stats
            .into_iter()
            .map(|(class_name, stats)| Record {
                class_name,
                virus_count: stats.virus_count as i32,
                assigned_taxids: stats.assigned_taxids.len() as i32,
                assigned_contigs: stats.assigned_contigs.len() as i32,
            })
            .collect_vec();

        match sort {
            SummarySortBy::ClassName => records.sort_by(|r1, r2| r1.class_name.cmp(&r2.class_name)),
            SummarySortBy::VirusCount => records.sort_by_key(|r| -r.virus_count),
            SummarySortBy::AssignedTaxIds => records.sort_by_key(|r| -r.assigned_taxids),
            SummarySortBy::AssignedContigs => records.sort_by_key(|r| -r.assigned_contigs),
        }

        let mut csv_writer = csv::WriterBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_writer(writer);

        for row in records {
            csv_writer.serialize(&row)?;
        }

        csv_writer.flush()?;
        Ok(())
    }
}

fn write_refined_output<'a, W, I>(writer: W, refinement: I) -> Result<()>
where
    W: io::Write,
    I: IntoIterator<Item = EnrichedVpfClassRecord<'a>>,
{
    let mut csv_writer = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(writer);

    for record in refinement {
        let serializable_record = CsvEnrichedVpfClassRecord::from(record);
        csv_writer.serialize(&serializable_record)?;
    }

    Ok(())
}

pub fn refine_vpf_class(args: RefineVpfClassArgs) -> Result<()> {
    if args.print_summary && args.output == "-" {
        anyhow::bail!(
            "--output - is incompatible with --print-summary. Use --output /path/to/file instead."
        );
    }

    let filters = VpfClassRecordFilter::parse_filters(&args.filter)
        .context("Error parsing --filter")?;

    let rank = args.get_or_infer_rank().ok_or_else(|| {
        anyhow::anyhow!(
            "--rank not specified and could not be deduced from {:?}",
            &args.vpf_class_prediction
        )
    })?;

    let taxonomy = PreprocessedTaxonomy::deserialize_with_format(
        &args.preprocessed_taxonomy,
        args.taxonomy_format,
    )?;

    with_some_taxonomy!(&taxonomy.tree, tax => {
        let rank_sym = tax
            .lookup_rank_sym(&rank)
            .ok_or_else(|| anyhow::anyhow!("Rank not found in taxonomy: {rank}"))?;

        let rank_assignments = load_rank_assignments(
            tax,
            rank_sym,
            args.include_descendants,
            args.metagenomic_assignment.as_ref(),
        )
        .context("Error collecting metagenomic assignment taxids")?;

        eprintln!(
            "Loaded assignments (dropped {} records out of {})",
            rank_assignments.dropped_records, rank_assignments.total_records
        );

        let csv_reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_path(&args.vpf_class_prediction)
            .context("Could not open VPF-Class prediction")?;

        itertools::process_results(csv_reader.into_deserialize(), |records| {
            let refinement = records
                .filter(|record| filters.iter().all(|f| f.apply(record)))
                .filter_map(|record| enrich_vpf_class_record(tax, rank_sym, &rank_assignments, args.allow_different_ranks, record));

            if args.print_summary {
                let mut summary = EnrichmentSummary::default();

                let refinement = refinement.inspect(|record| summary.account(record));

                writing_new_file_or_stdout!(&args.output, writer => {
                    write_refined_output(writer, refinement)
                        .context("Error writing refined output")?;
                });

                summary.write(io::stdout(), args.summary_sort_by)
                    .context("Error writing summary")?;

            } else {
                writing_new_file_or_stdout!(&args.output, writer => {
                    write_refined_output(writer, refinement)
                        .context("Error writing refined output")?;
                });
            }

            anyhow::Ok(())
        })??;
    });

    Ok(())
}
