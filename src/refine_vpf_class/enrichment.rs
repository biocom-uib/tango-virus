use std::collections::HashSet;

use itertools::Itertools;
use serde::{ser::SerializeStruct, Serialize, Serializer};

use crate::{
    crispr_match::CrisprMatchData,
    taxonomy::NodeId,
    util::csv_flatten_fix::{serialize_flat_struct, SerializeFlat},
    tool::vpf_class::VpfClassRecord,
};

pub trait Enrichment: Default {
    type Context;
    type CsvFields: SerializeFlat + From<Self>;

    type SummaryClassData: Default + Clone;
    type SummaryClassStats: SerializeFlat + From<Self::SummaryClassData>;

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

#[derive(Default, Clone)]
pub struct NoEnrichmentData;

impl SerializeFlat for NoEnrichmentData {
    const FIELD_COUNT: usize = 1;

    fn serialize_flat<Ser: SerializeStruct>(&self, _row: &mut Ser) -> Result<(), Ser::Error> {
        Ok(())
    }
}

#[derive(Default, Clone)]
pub struct NoEnrichment;

impl SerializeFlat for NoEnrichment {
    const FIELD_COUNT: usize = 1;

    fn serialize_flat<Ser: SerializeStruct>(&self, _row: &mut Ser) -> Result<(), Ser::Error> {
        Ok(())
    }
}

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

#[derive(Default, Clone, Serialize)]
pub struct CrisprEnrichment<'a> {
    pub crispr_matches: HashSet<&'a str>,
}

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

impl SerializeFlat for CrisprMatchesCsvFields {
    const FIELD_COUNT: usize = 1;

    fn serialize_flat<Ser: SerializeStruct>(&self, row: &mut Ser) -> Result<(), Ser::Error> {
        row.serialize_field("crispr_matches", &self.crispr_matches)
    }
}

#[derive(Default, Clone)]
pub struct CrisprEnrichmentSummaryClassData<'a> {
    pub crispr_match_contigs: HashSet<&'a str>,
}

#[derive(Default, Copy, Clone)]
pub struct CrisprEnrichmentSummaryClassStats {
    pub num_crispr_matched_contigs: usize,
}

impl<'a> From<CrisprEnrichmentSummaryClassData<'a>> for CrisprEnrichmentSummaryClassStats {
    fn from(value: CrisprEnrichmentSummaryClassData<'a>) -> Self {
        Self {
            num_crispr_matched_contigs: value.crispr_match_contigs.len(),
        }
    }
}

impl SerializeFlat for CrisprEnrichmentSummaryClassStats {
    const FIELD_COUNT: usize = 1;

    fn serialize_flat<Ser: SerializeStruct>(&self, row: &mut Ser) -> Result<(), Ser::Error> {
        row.serialize_field("num_crispr_matched_contigs", &self.num_crispr_matched_contigs)
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

#[derive(Debug)]
pub struct CsvEnrichedVpfClassRecord<S, CE: Enrichment> {
    pub vpf_class_record: VpfClassRecord<S>,

    pub assigned_taxids: String,
    pub num_assigned_contigs: usize,

    pub crispr_match_fields: CE::CsvFields,
}

impl<S: AsRef<str>, CE: Enrichment> Serialize for CsvEnrichedVpfClassRecord<S, CE> {
    fn serialize<Ser: Serializer>(&self, serializer: Ser) -> Result<Ser::Ok, Ser::Error> {
        serialize_flat_struct(serializer, "CsvEnrichedVpfClassRecord", self)
    }
}

impl<S, CE> SerializeFlat for CsvEnrichedVpfClassRecord<S, CE>
where
    S: AsRef<str>,
    CE: Enrichment,
{
    const FIELD_COUNT: usize = VpfClassRecord::<S>::FIELD_COUNT + 2 + CE::CsvFields::FIELD_COUNT;

    fn serialize_flat<Ser: SerializeStruct>(&self, row: &mut Ser) -> Result<(), Ser::Error> {
        self.vpf_class_record.serialize_flat(row)?;

        row.serialize_field("assigned_taxids", &self.assigned_taxids)?;
        row.serialize_field("num_assigned_contigs", &self.num_assigned_contigs)?;

        self.crispr_match_fields.serialize_flat(row)?;

        Ok(())
    }
}

impl<'a, S, CE> From<EnrichedVpfClassRecord<'a, S, CE>> for CsvEnrichedVpfClassRecord<S, CE>
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
