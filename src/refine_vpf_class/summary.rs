use std::{collections::{HashSet, HashMap}, io};

use clap::ValueEnum;
use itertools::Itertools;
use serde::{Serialize, ser::SerializeStruct, Serializer};

use crate::{taxonomy::NodeId, util::csv_flatten_fix::{SerializeFlat, serialize_flat_struct}};

use super::enrichment::{EnrichedVpfClassRecord, Enrichment};


struct ClassData<'a, CE: Enrichment> {
    virus_count: u32,
    assigned_taxids: HashSet<NodeId>,
    assigned_contigs: HashSet<&'a str>,
    crispr_enrichment_data: CE::SummaryClassData,
}

struct ClassStats<CE: Enrichment> {
    class_name: String,
    virus_count: i32,
    num_assigned_taxids: i32,
    num_assigned_contigs: i32,
    crispr_enrichment_stats: CE::SummaryClassStats,
}

impl<'a, CE: Enrichment> Default for ClassData<'a, CE> {
    fn default() -> Self {
        Self {
            virus_count: Default::default(),
            assigned_taxids: Default::default(),
            assigned_contigs: Default::default(),
            crispr_enrichment_data: Default::default(),
        }
    }
}

impl<'a, CE: Enrichment> ClassData<'a, CE> {
    fn add<S>(&mut self, record: &EnrichedVpfClassRecord<'a, S, CE>) {
        self.virus_count += 1;
        self.assigned_taxids.extend(&record.assigned_taxids);
        self.assigned_contigs.extend(&record.assigned_contigs);

        record.crispr_enrichment.add_class_data(&mut self.crispr_enrichment_data);
    }
}

impl<CE: Enrichment> SerializeFlat for ClassStats<CE> {
    const FIELD_COUNT: usize = 4 + CE::SummaryClassStats::FIELD_COUNT;

    fn serialize_flat<Ser: SerializeStruct>(&self, row: &mut Ser) -> Result<(), Ser::Error> {
        row.serialize_field("class_name", &self.class_name)?;
        row.serialize_field("virus_count", &self.virus_count)?;
        row.serialize_field("num_assigned_taxids", &self.num_assigned_taxids)?;
        row.serialize_field("num_assigned_contigs", &self.num_assigned_contigs)?;
        self.crispr_enrichment_stats.serialize_flat(row)?;
        Ok(())
    }
}

impl<CE: Enrichment> Serialize for ClassStats<CE> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serialize_flat_struct(serializer, "ClassStats", self)
    }
}


#[derive(ValueEnum, Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum SummarySortBy {
    ClassName,
    VirusCount,
    AssignedTaxIds,
    AssignedContigs,
}

pub struct EnrichmentSummary<'a, CE: Enrichment> {
    class_data: HashMap<String, ClassData<'a, CE>>,
}

impl<'a, CE: Enrichment> Default for EnrichmentSummary<'a, CE> {
    fn default() -> Self {
        Self { class_data: Default::default() }
    }
}

impl<'a, CE: Enrichment> EnrichmentSummary<'a, CE> {
    pub fn account<S: AsRef<str>>(&mut self, record: &EnrichedVpfClassRecord<'a, S, CE>) {
        let class_name = record.vpf_class_record.class_name.as_ref();

        self
            .class_data
            .raw_entry_mut()
            .from_key(class_name)
            .or_insert_with(|| (class_name.to_owned(), ClassData::default()))
            .1
            .add(record);
    }

    pub fn write<W: io::Write>(self, writer: W, sort: SummarySortBy) -> csv::Result<()> {
        let mut class_stats = self
            .class_data
            .into_iter()
            .map(|(class_name, data)| ClassStats::<CE> {
                class_name,
                virus_count: data.virus_count as i32,
                num_assigned_taxids: data.assigned_taxids.len() as i32,
                num_assigned_contigs: data.assigned_contigs.len() as i32,
                crispr_enrichment_stats: data.crispr_enrichment_data.into(),
            })
            .collect_vec();

        match sort {
            SummarySortBy::ClassName => class_stats.sort_by(|r1, r2| r1.class_name.cmp(&r2.class_name)),
            SummarySortBy::VirusCount => class_stats.sort_by_key(|r| -r.virus_count),
            SummarySortBy::AssignedTaxIds => class_stats.sort_by_key(|r| -r.num_assigned_taxids),
            SummarySortBy::AssignedContigs => class_stats.sort_by_key(|r| -r.num_assigned_contigs),
        }

        let mut csv_writer = csv::WriterBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_writer(writer);

        for row in class_stats {
            csv_writer.serialize(&row)?;
        }

        csv_writer.flush()?;
        Ok(())
    }
}
