use std::{
    collections::{HashMap, HashSet},
    io, iter,
    path::Path,
};

use anyhow::{Context, Result};
use clap::Args;
use lending_iterator::prelude::*;

use crate::{
    crispr_match::CrisprMatchData,
    preprocessed_taxonomy::{with_some_taxonomy, PreprocessedTaxonomyArgs},
    refine_vpf_class::enrichment::{CrisprEnrichment, NoEnrichment, NoEnrichmentContext},
    tango_assign::AssignmentRecord,
    taxonomy::{LabelledTaxonomy, NodeId, Taxonomy},
    tool::vpf_class::{self, VpfClassRecord, VpfClassRecordFilter, VpfClassRecordHKT},
    util::{self, csv_stream::CsvReaderIterExt, filter::FromStrFilter, writing_new_file_or_stdout},
};

use self::{
    enrichment::{CsvEnrichedVpfClassRecord, EnrichedVpfClassRecord, Enrichment},
    summary::{EnrichmentSummary, RecordDropReason, SummarySortBy},
};

pub(crate) mod enrichment;
pub(crate) mod summary;

/// Refine a VPF-Class host prediction to include only ascendants from taxonomic nodes found in a
/// metagenomic assignment.
#[derive(Args)]
pub struct RefineVpfClassArgs {
    #[clap(flatten)]
    taxonomy: PreprocessedTaxonomyArgs,

    /// TANGO metagenomic assignment output (see `meteor tango-assign`)
    metagenomic_assignment: String,

    /// Use CRISPR match data from this file (see `meteor crispr-match`) to restrict the refinement.
    #[clap(short, long)]
    crispr_matches: Option<String>,

    /// VPF-Class output file (presumably host_[rank].tsv)
    vpf_class_prediction: String,

    /// Does the presence of a higher-ranked node in the metagenomic assignment allow predictions
    /// about any of its descendants by itself? If given, all of its RANK descendants are assumed
    /// to be present in the sample.
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
    #[clap(short, long, default_value = "-")]
    output: String,

    /// Filters to apply before processing. Example: --filter 'membership_ratio>=0.2'
    #[clap(long)]
    filter: Vec<String>,

    /// Log every dropped record and the reason.
    #[clap(short, long)]
    verbose: bool,

    /// Print a summary of class frequencies and their average membership ratio (use '-' for
    /// STDOUT)
    #[clap(long, value_name = "SUMMARY_FILE")]
    write_class_summary: Option<String>,

    ///// Column to sort --write-summary by
    #[clap(long, value_enum, default_value_t = SummarySortBy::VirusCount)]
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
    // node at rank -> related nodes in assignment (depending on --include-descendants)
    rank_assignments: HashMap<NodeId, HashSet<NodeId>>,
    // node at rank -> related contigs (depending on --include-descendants)
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
    verbose: bool,
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
            .chain(tax.strict_ancestors(node))
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

                if verbose {
                    eprintln!(
                        "Warning: No valid descendants found for assignment {:?} (taxid:{}) (rank: {})",
                        record.assigned_name,
                        node,
                        tax.find_rank(node)
                            .and_then(|sym| tax.rank_sym_str(sym))
                            .unwrap_or(""),
                    );
                }
            }
        } else {
            assignments.dropped_records += 1;

            if verbose {
                eprintln!(
                    "Warning: No valid ancestors found for assignment {:?} (taxid:{}) (rank: {})",
                    record.assigned_name,
                    node,
                    tax.find_rank(node)
                        .and_then(|sym| tax.rank_sym_str(sym))
                        .unwrap_or(""),
                );
            }
        }

        assignments.total_records += 1;
    }

    Ok(assignments)
}

struct RefinementContext<Tax: Taxonomy> {
    taxonomy: Tax,
    verbose: bool,
    rank_sym: Tax::RankSym,
    rank_assignments: RankAssignments,
    filters: Vec<VpfClassRecordFilter>,
}

impl<Tax: LabelledTaxonomy> RefinementContext<Tax> {
    fn enrich_vpf_class_record<'a, S: AsRef<str>, CE: Enrichment>(
        &'a self,
        allow_different_ranks: bool,
        crispr_context: &CE::Context,
        record: VpfClassRecord<S>,
    ) -> Result<EnrichedVpfClassRecord<'a, S, CE>, RecordDropReason> {
        let mut enriched = EnrichedVpfClassRecord::<'a, S, CE>::new(record);

        let mut found_mismatched_ranks = false;
        let mut found_taxids = false;

        for node in self.taxonomy.nodes_with_label(enriched.vpf_class_record.class_name.as_ref()) {
            found_taxids = true;

            if let Some((taxids, contigs)) = self.rank_assignments.lookup_rank_node(node) {
                if allow_different_ranks || self.taxonomy.find_rank(node) == Some(self.rank_sym) {
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
                "Warning: No taxids known for host prediction {}",
                enriched.vpf_class_record.class_name.as_ref(),
            );
            return Err(RecordDropReason::UnknownHostTaxId);
        }

        if found_mismatched_ranks {
            eprintln!(
                "Warning: Found taxids for host prediction {} but had mismatching ranks, skipping (use --allow-different-ranks to disable this check)",
                enriched.vpf_class_record.class_name.as_ref(),
            );
        }

        if !enriched.has_assignments() {
            return Err(RecordDropReason::NotInAssignment);
        }

        let crisprs_found = enriched.crispr_enrichment.enrich(
            crispr_context,
            &enriched.vpf_class_record,
            &enriched.assigned_taxids,
            &enriched.assigned_contigs,
        );

        if !crisprs_found {
            Err(RecordDropReason::NoCrisprInfo)
        } else {
            Ok(enriched)
        }
    }

    #[apply(Gat!)]
    fn process_vpf_class_records<CE: Enrichment, Iter, W: io::Write>(
        &self,
        args: &RefineVpfClassArgs,
        crispr_context: &CE::Context,
        mut records: Iter,
        output_writer: W,
    ) -> anyhow::Result<EnrichmentSummary<'_, CE>>
    where
        Iter: for<'n> LendingIterator<Item<'n> = csv::Result<VpfClassRecord<&'n str>>>,
    {
        let mut csv_writer = csv::WriterBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_writer(output_writer);

        let filters = &self.filters;
        let mut summary = EnrichmentSummary::<CE>::default();

        while let Some(record) = records.next() {
            let record = record?;

            if !self.verbose && summary.num_records() % 100 == 0 {
                eprint!("\rProcessing record #{}", summary.num_records());
            }

            if !filters.iter().all(|f| f.apply(&record)) {
                summary.account_dropped(RecordDropReason::Filtered, self.verbose);
                continue;
            }

            let enriched = match self.enrich_vpf_class_record::<&str, CE>(args.allow_different_ranks, crispr_context, record) {
                Ok(enriched) => {
                    summary.account_kept(self.verbose);
                    enriched
                },
                Err(drop_reason) => {
                    summary.account_dropped(drop_reason, self.verbose);
                    continue;
                }
            };

            if args.write_class_summary.is_some() {
                summary.account_classes(&enriched);
            }

            if let Err(e) = csv_writer.serialize(CsvEnrichedVpfClassRecord::<&str, CE>::from(enriched)) {
                if util::is_broken_pipe(&e) {
                    return Ok(summary);
                } else {
                    return Err(e.into());
                }
            }
        }

        util::ignore_broken_pipe(csv_writer.flush())?;

        if !self.verbose {
            eprintln!();
        }

        Ok(summary)
    }
}

pub fn refine_vpf_class(args: RefineVpfClassArgs) -> Result<()> {
    if args.write_class_summary.as_deref() == Some("-") && args.output == "-" {
        anyhow::bail!(
            "Both --output and --write-summary can't be stdout."
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

    let taxonomy = args.taxonomy.deserialize()?;

    with_some_taxonomy!(&taxonomy.tree, tax => {
        let rank_sym = tax
            .lookup_rank_sym(&rank)
            .ok_or_else(|| anyhow::anyhow!("Rank not found in taxonomy: {rank}"))?;

        let rank_assignments = load_rank_assignments(
            tax,
            rank_sym,
            args.include_descendants,
            args.verbose,
            args.metagenomic_assignment.as_ref(),
        )
        .context("Error collecting metagenomic assignment taxids")?;

        eprintln!(
            "Loaded assignments (dropped {} records out of {})",
            rank_assignments.dropped_records, rank_assignments.total_records
        );

        let context = RefinementContext {
            taxonomy: tax,
            rank_sym,
            verbose: args.verbose,
            rank_assignments, filters
        };

        let crispr_matches = if let Some(crispr_match_path) = &args.crispr_matches {
            let data = CrisprMatchData::load(crispr_match_path.as_ref())
                .context("Loading CRISPR match data")?;

            eprintln!(
                "Loaded CRISPR match data ({} viruses and {} hosts present)",
                data.crispr_matches.0.key_interner().len(),
                data.crispr_matches.0.value_interner().len()
            );

            Some(data)
        } else {
            None
        };

        let csv_records = vpf_class::vpf_class_records_reader(args.vpf_class_prediction.as_ref())
            .context("Could not open VPF-Class prediction")?
            .into_lending_iter()
            .into_deserialize::<VpfClassRecordHKT>(None);

        macro_rules! write_summaries {
            ($summary:ident) => {{
                if let Some(summary_file) = &args.write_class_summary {
                    writing_new_file_or_stdout!(summary_file, writer => {
                        let writer = writer.context("Error creating summary file")?;

                        util::ignore_broken_pipe_anyhow(
                            $summary.write_classes(writer, args.summary_sort_by)
                                .context("Error writing summary")
                        )?;
                    });
                }

                $summary.write_drop_stats(io::stderr())?;
            }}
        }

        writing_new_file_or_stdout!(&args.output, writer => {
            let writer = writer.context("Error creating output file")?;

            match crispr_matches {
                Some(crispr_matches) => {
                    let summary = context.process_vpf_class_records::<CrisprEnrichment, _, _>(&args, &&crispr_matches, csv_records, writer)?;
                    write_summaries!(summary);
                }
                None => {
                    let summary = context.process_vpf_class_records::<NoEnrichment, _, _>(&args, &NoEnrichmentContext, csv_records, writer)?;
                    write_summaries!(summary);
                }
            }
        });
    });

    Ok(())
}
