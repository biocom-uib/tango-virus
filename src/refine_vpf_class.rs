use std::{collections::HashSet, iter, path::Path, fs::File, io, error};

use anyhow::{Result, Context};
use clap::Args;
use itertools::Either;
use serde::{Serialize, Deserialize};

use crate::{
    preprocess_taxonomy::{PreprocessedTaxonomy, PreprocessedTaxonomyFormat, SomeTaxonomy},
    taxonomy::{NodeId, LabelledTaxonomy, Taxonomy}, assign::AssignmentRecord,
};

#[derive(Args)]
pub struct RefineVpfClassArgs {
    #[clap(long, arg_enum, default_value_t = PreprocessedTaxonomyFormat::CBOR)]
    taxonomy_format: PreprocessedTaxonomyFormat,

    /// Path to the preprocessed taxonomy (presumably from preprocess-taxonomy)
    preprocessed_taxonomy: String,

    /// TANGO3 metagenomic assignment output
    metagenomic_assignment: String,

    /// VPF-Class output file (presumably host_[rank].tsv)
    vpf_class_prediction: String,

    /// Does the presence of a higher-ranked node in the metagenomic assignment imply the presence
    /// of any lower node?
    #[clap(long)]
    include_descendants: bool,

    /// When mapping VPF classifications to taxids, ignore their rank.
    #[clap(long)]
    allow_different_ranks: bool,

    /// VPF-Class output's taxonomic rank. If not specified, will be inferred from
    /// `VPF_CLASS_PREDICTION`'s file name
    #[clap(short, long)]
    rank: Option<String>,

    /// Refined VPF-Class output file (use '-' for STDOUT)
    #[clap(short, long)]
    output: String
}

impl RefineVpfClassArgs {
    fn get_or_infer_rank(&self) -> Option<&str> {
        self.rank.as_deref().or_else(|| {
            self.vpf_class_prediction
                .strip_prefix("host_")
                .and_then(|s| s.strip_suffix(".tsv"))
        })
    }
}

fn load_assignment_found_taxids<Tax: Taxonomy, P: AsRef<Path>>(
    tax: &Tax,
    rank: &str,
    include_descendants: bool,
    assignments_path: P,
) -> Result<(Tax::RankSym, HashSet<NodeId>)> {
    let rank_sym = tax
        .lookup_rank_sym(rank)
        .ok_or_else(|| anyhow::anyhow!("Given rank not found in taxonomy"))?;

    let csv_reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(assignments_path)?;

    let mut present = HashSet::new();

    for record in csv_reader.into_deserialize() {
        let record: AssignmentRecord = record?;
        let node = NodeId(record.assigned_taxid);

        let valid_ancestor = iter::once(node)
            .chain(tax.ancestors(node))
            .find(|&ancestor| rank_sym == tax.get_rank(ancestor));

        if let Some(valid_ancestor) = valid_ancestor {
            present.insert(valid_ancestor);
        } else if include_descendants {
            let mut valid_descendants = tax
                .postorder_descendants(node)
                .filter(|&desc| rank_sym == tax.get_rank(desc))
                .peekable();

            if valid_descendants.peek().is_some() {
                present.extend(valid_descendants);
            } else {
                eprintln!(
                    "Warning: No valid descendants found for {:?} (taxid:{})",
                    record.assigned_name, node
                );
            }
        } else {
            eprintln!(
                "Warning: No valid ancestors found for {:?} (taxid:{})",
                record.assigned_name, node
            );
        }
    }

    Ok((rank_sym, present))
}

#[derive(Serialize, Deserialize)]
struct VpfClassRecord {
    virus_name: String,
    class_name: String,
    membership_ratio: f64,
    virus_hit_score: f64,
    confidence_score: f64,
}

fn open_and_refine_vpf_class<'a, Tax: LabelledTaxonomy<RankSym: 'static>, P: AsRef<Path>>(
    tax: &'a Tax,
    rank_sym: Tax::RankSym,
    present: &'a HashSet<NodeId>,
    allow_different_ranks: bool,
    path: P,
) -> Result<impl Iterator<Item = csv::Result<VpfClassRecord>> + 'a> {
    let csv_reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(path)?;

    let result = csv_reader.into_deserialize().filter(move |record| {
        if let Ok(record) = record {
            let record: &VpfClassRecord = &record;

            let mut had_candidate: u32 = 0;

            for node in tax.nodes_with_label(&record.class_name) {
                if present.contains(&node) {
                    if allow_different_ranks || tax.get_rank(node) == rank_sym {
                        return true;
                    } else {
                        had_candidate |= 0b10;
                    }
                } else {
                    had_candidate |= 0b1;
                }
            }

            if had_candidate & 0b10 != 0 {
                eprintln!(
                    "Warning: Found taxids for classification {} but had mismatching ranks, skipping (use --allow-different-ranks to disable the check)",
                    &record.class_name
                );
            } else if had_candidate & 0b01 == 0 {
                eprintln!("Warning: No taxids found for classification {}", &record.class_name);
            }

            false

        } else {
            true
        }
    });

    Ok(result)
}

fn write_refined_output<W, I, E: 'static>(writer: W, refinement: I) -> Result<()>
where
    W: io::Write,
    I: Iterator<Item = Result<VpfClassRecord, E>>,
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
    let taxonomy = PreprocessedTaxonomy::deserialize_with_format(
        &args.preprocessed_taxonomy,
        args.taxonomy_format
    )?;

    let rank = args.get_or_infer_rank().ok_or_else(|| {
        anyhow::anyhow!(
            "--rank not specified and could not be deduced from {:?}",
            &args.vpf_class_prediction
        )
    })?;

    let present;

    let refinement = match &taxonomy.tree {
        SomeTaxonomy::NcbiTaxonomyWithSingleClassNames(tax) => {
            let rank_sym;

            (rank_sym, present) = load_assignment_found_taxids(
                tax,
                rank,
                args.include_descendants,
                &args.metagenomic_assignment,
            )
            .context("Error collecting metagenomic assignment taxids")?;

            Either::Left(
                open_and_refine_vpf_class(
                    tax,
                    rank_sym,
                    &present,
                    args.allow_different_ranks,
                    &args.vpf_class_prediction,
                )
                .context("Error reading VPF-Class output")?,
            )
        }
        SomeTaxonomy::NewickTaxonomy(tax) => {
            let rank_sym;

            (rank_sym, present) = load_assignment_found_taxids(
                tax,
                rank,
                args.include_descendants,
                &args.metagenomic_assignment,
            )
            .context("Error collecting metagenomic assignment taxids")?;

            Either::Right(
                open_and_refine_vpf_class(
                    tax,
                    rank_sym,
                    &present,
                    args.allow_different_ranks,
                    &args.vpf_class_prediction,
                )
                .context("Error reading VPF-Class output")?,
            )
        }
    };

    if args.output == "-" {
        write_refined_output(std::io::stdout(), refinement)?;
    } else {
        write_refined_output(File::create(&args.output)?, refinement)?;
    }

    Ok(())
}
