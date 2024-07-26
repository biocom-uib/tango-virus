use std::path::Path;

use anyhow::{anyhow, Context};
use clap::{Args, ValueEnum};

use itertools::{izip, Itertools};
use polars::datatypes::{DataType, ListChunked};
use polars::lazy::dsl::{col, Expr, GetOutput};
use polars::prelude::{JoinType, JoinValidation, LazyFrame};
use polars::series::{IntoSeries, Series};
use crate::preprocessed_taxonomy::PreprocessedTaxonomyArgs;
use crate::taxonomy::formats::ncbi::{NcbiTaxonomy, TaxId};
use crate::taxonomy::{Taxonomy, LabelledTaxonomy, NodeId};
use crate::tool::vpf_class;

#[derive(Copy, Clone, Eq, PartialEq, PartialOrd, Ord, Hash, ValueEnum)]
pub enum VpfClassRank {
    Family,
    Genus
}

#[derive(Args)]
pub struct JoinVpfClassArgs {
    /// Path to the VPF-Class output directory
    #[clap(long, help_heading = "VPF-Class")]
    vpf_class: String,

    /// Minimum required membership_ratio of a prediction
    #[clap(long, help_heading = "VPF-Class", requires = "vpf_class")]
    vpf_class_membership_ratio: Option<f32>,

    /// Minimum required confidence_score of a prediction
    #[clap(long, help_heading = "VPF-Class", requires = "vpf_class")]
    vpf_class_confidence_score: Option<f32>,

    /// Ranks to use from VPF-Class (family and/or genus)
    #[clap(long, value_enum, value_delimiter = ',', help_heading = "VPF-Class")]
    vpf_class_rank: Vec<VpfClassRank>,
}

/// Combine virus taxonomy annotations from multiple tools
#[derive(Args)]
pub struct JoinVirusTaxonomyArgs {
    #[clap(flatten)]
    taxonomy: PreprocessedTaxonomyArgs,

    /// Output path (use '-' for STDOUT)
    #[clap(short, long, default_value = "-")]
    output: String,

    #[clap(flatten)]
    vpf_class_args: Option<JoinVpfClassArgs>,
}

fn vpf_class_filter_expr(args: &JoinVpfClassArgs) -> Option<Expr> {
    let mut filter = None;

    if let Some(membership_ratio) = args.vpf_class_membership_ratio {
        filter = Some(col("membership_ratio").gt_eq(membership_ratio));
    }

    if let Some(confidence_score) = args.vpf_class_confidence_score {
        let cf_filter = col("confidence_score").gt_eq(confidence_score);

        filter = if let Some(mr_filter) = filter {
            Some(mr_filter.and(cf_filter))
        } else {
            Some(cf_filter)
        };
    }

    filter
}

fn lookup_taxid<Tax>(tax: &Tax, rank_sym: Tax::RankSym, name: &str) -> Option<NodeId>
where
    Tax: LabelledTaxonomy,
{
    let mut candidate = None;

    for node in tax.nodes_with_label(name) {
        if tax.find_rank(node) == Some(rank_sym) {
            return Some(node)
        } else if candidate.is_none() {
            candidate = Some(node);
        } else {
            eprintln!("Found more than one matching taxid for {name}");
        }
    }

    candidate
}

fn join_vpf_class_ranks<Tax>(tax: &Tax, family: LazyFrame, genus: LazyFrame) -> LazyFrame
where
    Tax: LabelledTaxonomy,
{
    let joined = family
        .join_builder()
        .with(genus)
        .how(JoinType::Full)
        .validate(JoinValidation::ManyToMany)
        .on([col("virus_name")])
        .coalesce(polars::prelude::JoinCoalesce::CoalesceColumns)
        .finish();

    joined
        .with_column(
            col("family_taxids")
                .map_many(
                    |params| {
                        let family_taxidss = params[0].list()?;
                        let genus_taxidss = params[1].list()?;

                        let result = izip!(family_taxidss, genus_taxidss)
                            .map(|(family_taxids, _genus_taxids)| {
                                family_taxids.unwrap()
                            })
                            .collect::<ListChunked>()
                            .into_series();

                        Ok(Some(result))
                    },
                    &[col("genus_taxids")],
                    GetOutput::from_type(DataType::List(Box::new(DataType::UInt64))),
                )
                .alias("taxids"),
        )
        .select([col("virus_name"), col("taxids")])
}

fn load_vpf_class<Tax>(args: &JoinVpfClassArgs, tax: &Tax) -> anyhow::Result<()>
where
    Tax: LabelledTaxonomy,
{
    let path: &Path = args.vpf_class.as_ref();

    let filter = vpf_class_filter_expr(args);

    let ranks = &args.vpf_class_rank;

    let load_rank = |rank| {
        let mut prediction = vpf_class::load_as_lazyframe(&path.join(format!("{}.tsv", rank)))
            .with_context(|| format!("loading VPF-Class {} prediction", rank))?;

        if let Some(filter) = &filter {
            prediction = prediction.filter(filter.clone());
        }

        let rank_sym = tax
            .lookup_rank_sym(rank)
            .ok_or_else(|| anyhow!("Could not find rank {rank} in the taxonomy"))?;

        let taxid_col = format!("{}_taxid", rank);

        prediction = prediction
            .select([col("virus_name"), col("class_name")])
            .with_column(
                col("class_name")
                    .map(
                        |class_name| {
                            let result = class_name
                                .str()?
                                .iter()
                                .map(|class| Some(u64::from(lookup_taxid(tax, rank_sym, class?)?)))
                                .collect();

                            Ok(Some(result))
                        },
                        GetOutput::from_type(DataType::UInt64),
                    )
            )
            .rename(["class_name"], [&taxid_col])
            .group_by([col("virus_name")])
            .agg([col(&taxid_col).implode()]);

        anyhow::Ok(prediction)
    };

    let family = if ranks.iter().contains(&VpfClassRank::Family) {
        Some(load_rank("family")?)
    } else {
        None
    };

    let genus = if ranks.iter().contains(&VpfClassRank::Genus) {
        Some(load_rank("genus")?)
    } else {
        None
    };

    match (family, genus) {
        (Some(family), None) => family,
        (None, Some(genus)) => genus,
        (Some(family), Some(genus)) => join_vpf_class_ranks(family, genus),
    }

    Ok(())
}

fn join_virus_taxonomy(args: JoinVirusTaxonomyArgs) {

}
