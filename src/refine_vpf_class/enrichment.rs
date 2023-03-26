use std::collections::HashSet;

use itertools::Itertools;
use serde::Serialize;

use crate::{util::vpf_class_record::VpfClassRecord, taxonomy::NodeId};

use super::CrisprMatchData;


pub trait Enrichment: Default {
    type Context;
    type CsvFields: Serialize + From<Self>;

    type SummaryClassData: Default;
    type SummaryClassStats: Serialize + From<Self::SummaryClassData>;

    fn enrich<S: AsRef<str>>(
        &mut self,
        context: &Self::Context,
        vpf_class_record: &VpfClassRecord<S>,
        assigned_taxids: &HashSet<NodeId>,
        assigned_contigs: &HashSet<&str>,
    ) -> bool;

    fn add_class_data(&self, class_data: &mut Self::SummaryClassData);

    fn into_csv_fields(self) -> Self::CsvFields {
        self.into()
    }
}

pub struct NoEnrichmentContext;

#[derive(Default, Serialize)]
pub struct NoEnrichmentData;

#[derive(Default, Serialize)]
pub struct NoEnrichment {}

impl Enrichment for NoEnrichment {
    type Context = NoEnrichmentContext;
    type CsvFields = NoEnrichment;
    type SummaryClassData = NoEnrichmentData;
    type SummaryClassStats = NoEnrichmentData;

    #[must_use]
    fn enrich<S: AsRef<str>>(
        &mut self,
        _context: &Self::Context,
        _vpf_class_record: &VpfClassRecord<S>,
        _assigned_taxids: &HashSet<NodeId>,
        _assigned_contigs: &HashSet<&str>,
    ) -> bool {
        true
    }

    fn add_class_data(&self, _class_data: &mut Self::SummaryClassData) {
    }
}

#[derive(Default, Serialize)]
pub struct CrisprEnrichment<'a> {
    pub crispr_matches: HashSet<&'a str>,
}

#[derive(Serialize)]
pub struct CrisprMatchesCsvFields {
    pub crispr_matches: String,
}

impl<'a> From<CrisprEnrichment<'a>> for CrisprMatchesCsvFields {
    fn from(value: CrisprEnrichment) -> Self {
        Self {
            crispr_matches: value.crispr_matches.iter().join(";"),
        }
    }
}

#[derive(Default)]
pub struct CrisprEnrichmentSummaryClassData<'a> {
    pub crispr_match_contigs: HashSet<&'a str>,
}

#[derive(Default, Serialize)]
pub struct CrisprEnrichmentSummaryClassStats {
    pub crispr_match_contigs: String,
}

impl<'a> From<CrisprEnrichmentSummaryClassData<'a>> for CrisprEnrichmentSummaryClassStats {
    fn from(value: CrisprEnrichmentSummaryClassData<'a>) -> Self {
        Self {
            crispr_match_contigs: value.crispr_match_contigs.iter().join(";"),
        }
    }
}

impl<'a> Enrichment for CrisprEnrichment<'a> {
    type Context = &'a CrisprMatchData;
    type CsvFields = CrisprMatchesCsvFields;

    fn enrich<S: AsRef<str>>(
        &mut self,
        context: &&'a CrisprMatchData,
        vpf_class_record: &VpfClassRecord<S>,
        _assigned_taxids: &HashSet<NodeId>,
        assigned_contigs: &HashSet<&str>,
    ) -> bool {
        self.crispr_matches.extend(
            context.crispr_matches
                .0
                .lookup(vpf_class_record.virus_name.as_ref())
                .filter(|contig| assigned_contigs.contains(contig))
        );

        !self.crispr_matches.is_empty()
    }

    type SummaryClassData = CrisprEnrichmentSummaryClassData<'a>;
    type SummaryClassStats = CrisprEnrichmentSummaryClassStats;

    fn add_class_data(&self, class_data: &mut Self::SummaryClassData) {
        class_data.crispr_match_contigs.extend(&self.crispr_matches);
    }
}


pub struct EnrichedVpfClassRecord<'a, S, CE> {
    pub vpf_class_record: VpfClassRecord<S>,
    pub assigned_taxids: HashSet<NodeId>,
    pub assigned_contigs: HashSet<&'a str>,
    pub crispr_enrichment: CE,
}

impl<'a, S, CE: Default> EnrichedVpfClassRecord<'a, S, CE> {
    pub fn new(vpf_class_record: VpfClassRecord<S>) -> Self {
        EnrichedVpfClassRecord {
            vpf_class_record,
            assigned_taxids: HashSet::new(),
            assigned_contigs: HashSet::new(),
            crispr_enrichment: CE::default(),
        }
    }

    pub fn has_assignments(&self) -> bool {
        !self.assigned_taxids.is_empty() || !self.assigned_contigs.is_empty()
    }
}

#[derive(Debug, Serialize)]
pub struct CsvEnrichedVpfClassRecord<S, CrisprMatchFields> {
    #[serde(flatten)]
    pub vpf_class_record: VpfClassRecord<S>,

    pub assigned_taxids: String,
    pub num_assigned_contigs: usize,

    #[serde(flatten)]
    pub crispr_match_fields: CrisprMatchFields,
}

impl<'a, S, CE> From<EnrichedVpfClassRecord<'a, S, CE>>
    for CsvEnrichedVpfClassRecord<S, CE::CsvFields>
where
    CE: Enrichment,
{
    fn from(other: EnrichedVpfClassRecord<'a, S, CE>) -> Self {
        CsvEnrichedVpfClassRecord {
            vpf_class_record: other.vpf_class_record,
            assigned_taxids: other.assigned_taxids.iter().join(";"),
            num_assigned_contigs: other.assigned_contigs.len(),
            crispr_match_fields: other.crispr_enrichment.into_csv_fields(),
        }
    }
}
