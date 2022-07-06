use std::{fs::File, path::Path};

use anyhow::Context;
use clap::{Args, ValueEnum};
use serde::{de::Error, Deserialize, Serialize, Deserializer};

use crate::taxonomy::formats::{newick::NewickTaxonomy, ncbi::{self, NcbiTaxonomy}};

#[derive(Deserialize, Serialize)]
pub enum SomeTaxonomy {
    NcbiTaxonomyWithSingleClassNames(NcbiTaxonomy<ncbi::SingleClassNames>),
    NcbiTaxonomyWithManyNames(NcbiTaxonomy<ncbi::AllNames>),
    NewickTaxonomy(NewickTaxonomy),
}

#[derive(Deserialize, Serialize)]
pub struct PreprocessedTaxonomy {
    #[serde(deserialize_with = "deserialize_check_version")]
    version: u64,
    pub tree: SomeTaxonomy,
    pub ordered_ranks: Option<Vec<String>>,
}

impl Default for PreprocessedTaxonomyFormat {
    fn default() -> Self {
        Self::Bincode
    }
}

macro_rules! with_some_taxonomy {
    ($some_taxonomy:expr, $tax:pat => $body:expr $(,)?) => {
        match $some_taxonomy {
            crate::preprocessed_taxonomy::SomeTaxonomy::NcbiTaxonomyWithSingleClassNames($tax) => $body,
            crate::preprocessed_taxonomy::SomeTaxonomy::NcbiTaxonomyWithManyNames($tax) => $body,
            crate::preprocessed_taxonomy::SomeTaxonomy::NewickTaxonomy($tax) => $body,
        }
    };
}

macro_rules! with_some_ncbi_or_newick_taxonomy {
    (
        $some_taxonomy:expr,
        ncbi: $ncbi_tax:pat => $ncbi_body:expr,
        newick: $newick_tax:pat => $newick_body:expr
        $(,)?
    ) => {
        match $some_taxonomy {
            SomeTaxonomy::NcbiTaxonomyWithSingleClassNames($ncbi_tax) => $ncbi_body,
            SomeTaxonomy::NcbiTaxonomyWithManyNames($ncbi_tax) => $ncbi_body,
            SomeTaxonomy::NewickTaxonomy($newick_tax) => $newick_body,
        }
    };
}

pub(crate) use with_some_ncbi_or_newick_taxonomy;
pub(crate) use with_some_taxonomy;

#[derive(ValueEnum, Debug, Copy, Clone)]
pub enum PreprocessedTaxonomyFormat {
    Bincode,
    #[cfg(feature = "taxonomy-serialize-cbor")]
    Cbor,
    #[cfg(feature = "taxonomy-serialize-json")]
    Json,
}

fn deserialize_check_version<'de, D: Deserializer<'de>>(d: D) -> Result<u64, D::Error> {
    let version = u64::deserialize(d)?;

    if version != PreprocessedTaxonomy::FORMAT_VERSION {
        let msg = format!("Preprocessed taxonomy version mismatch (found: {}, current: {}). Re-generate it using preprocess-taxonomy.", version, PreprocessedTaxonomy::FORMAT_VERSION);
        return Err(D::Error::custom(msg));
    }

    Ok(version)
}

impl PreprocessedTaxonomy {
    pub const FORMAT_VERSION: u64 = 1
        + NcbiTaxonomy::<ncbi::SingleClassNames>::FORMAT_VERSION as u64
        + NcbiTaxonomy::<ncbi::AllNames>::FORMAT_VERSION as u64
        + NewickTaxonomy::FORMAT_VERSION as u64;

    pub fn new(tree: SomeTaxonomy, ordered_ranks: Option<Vec<String>>) -> Self {
        Self {
            version: Self::FORMAT_VERSION,
            tree,
            ordered_ranks,
        }
    }

    pub fn deserialize_with_format<P: AsRef<Path>>(
        path: P,
        format: PreprocessedTaxonomyFormat,
    ) -> anyhow::Result<Self> {
        match format {
            PreprocessedTaxonomyFormat::Bincode => Ok(bincode::deserialize_from(File::open(path)?)?),

            #[cfg(feature ="taxonomy-serialize-cbor")]
            PreprocessedTaxonomyFormat::Cbor => Ok(serde_cbor::from_reader(File::open(path)?)?),

            #[cfg(feature ="taxonomy-serialize-json")]
            PreprocessedTaxonomyFormat::Json => Ok(serde_json::from_reader(File::open(path)?)?),
        }
    }

    pub fn serialize_with_format<P: AsRef<Path>>(
        &self,
        path: P,
        format: PreprocessedTaxonomyFormat,
    ) -> anyhow::Result<()> {
        match format {
            PreprocessedTaxonomyFormat::Bincode => {
                Ok(bincode::serialize_into(File::create(path)?, self)?)
            }
            #[cfg(feature ="taxonomy-serialize-cbor")]
            PreprocessedTaxonomyFormat::Cbor => {
                Ok(serde_cbor::to_writer(File::create(path)?, self)?)
            }
            #[cfg(feature ="taxonomy-serialize-json")]
            PreprocessedTaxonomyFormat::Json => {
                Ok(serde_json::to_writer(File::create(path)?, self)?)
            }
        }
    }
}

#[derive(Args)]
pub struct PreprocessedTaxonomyArgs {
    /// Format in which the preprocessed taxonomy was serialized.
    #[clap(long, value_enum, default_value_t)]
    taxonomy_format: PreprocessedTaxonomyFormat,

    /// Path to the preprocessed taxonomy (from preprocess-taxonomy)
    preprocessed_taxonomy: String,
}

impl PreprocessedTaxonomyArgs {
    pub fn deserialize(&self) -> anyhow::Result<PreprocessedTaxonomy> {
        PreprocessedTaxonomy::deserialize_with_format(
            &self.preprocessed_taxonomy,
            self.taxonomy_format
        )
        .context("Could not deserialize the preprocessed taxonomy")
    }
}
