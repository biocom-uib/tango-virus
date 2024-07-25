use std::{num::ParseFloatError, path::Path, fs::File};

use lending_iterator::HKT;
use serde::{Serialize, Deserialize, ser::SerializeStruct};
use thiserror::Error;

use crate::util::{filter::{self, FromStrFilter}, csv_flatten_fix::SerializeFlat};

#[derive(Debug, Serialize, Deserialize)]
pub struct VpfClassRecord<S = String> {
    pub virus_name: S,
    pub class_name: S,
    pub membership_ratio: f64,
    pub virus_hit_score: f64,
    pub confidence_score: f64,
}

impl<S: AsRef<str>> SerializeFlat for VpfClassRecord<S> {
    const FIELD_COUNT: usize = 5;

    fn serialize_flat<Ser>(&self, row: &mut Ser) -> Result<(), Ser::Error>
    where
        Ser: SerializeStruct,
    {
        row.serialize_field("virus_name", self.virus_name.as_ref())?;
        row.serialize_field("class_name", self.class_name.as_ref())?;
        row.serialize_field("membership_ratio", &self.membership_ratio)?;
        row.serialize_field("virus_hit_score", &self.virus_hit_score)?;
        row.serialize_field("confidence_score", &self.confidence_score)?;
        Ok(())
    }
}

pub type VpfClassRecordHKT = HKT!(VpfClassRecord<&str>);

pub enum VpfClassRecordFieldValue {
    VirusName(String),
    ClassName(String),
    MembershipRatio(f64),
    VirusHitScore(f64),
    ConfidenceScore(f64),
}

pub struct VpfClassRecordFilter {
    field_value: VpfClassRecordFieldValue,
    op: filter::Op,
}

impl VpfClassRecordFilter {
    pub fn apply<S: AsRef<str>>(&self, record: &VpfClassRecord<S>) -> bool {
        use VpfClassRecordFieldValue::*;

        match &self.field_value {
            VirusName(value) => self.op.apply(record.virus_name.as_ref(), value),
            ClassName(value) => self.op.apply(record.class_name.as_ref(), value),
            MembershipRatio(value) => self.op.apply(&record.membership_ratio, value),
            VirusHitScore(value) => self.op.apply(&record.virus_hit_score, value),
            ConfidenceScore(value) => self.op.apply(&record.confidence_score, value),
        }
    }
}

#[derive(Debug, Error)]
pub enum VpfClassRecordFilterParseError {
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

pub fn vpf_class_records_reader(path: &Path) -> csv::Result<csv::Reader<File>> {
    csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(path)
}

pub fn vpf_class_record_schema() -> polars::prelude::Schema {
    use polars::prelude::{DataType, Schema};

    let mut schema = Schema::new();

    schema.with_column("virus_name".into(), DataType::String);
    schema.with_column("class_name".into(), DataType::String);
    schema.with_column("membership_ratio".into(), DataType::Float32);
    schema.with_column("virus_hit_score".into(), DataType::Float32);
    schema.with_column("confidence_score".into(), DataType::Float32);

    schema
}
