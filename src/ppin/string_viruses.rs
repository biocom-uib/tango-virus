use std::{path::Path, collections::{HashMap, HashSet}, fs::File, io::{self, BufReader}};

use clap::Args;
use flate2::bufread::GzDecoder;
use polars::prelude::LazyFrame;
use string_interner::{StringInterner, DefaultSymbol};

use crate::{
    preprocessed_taxonomy::{with_some_ncbi_or_newick_taxonomy, PreprocessedTaxonomyArgs},
    tango_assign::AssignmentRecord,
    taxonomy::Taxonomy,
    util::blastout::{self, BlastOutFilter, BlastOutFmt},
    util::filter::FromStrFilter,
};


struct StringProteinAliases<Entry> {
    string_ids: StringInterner,
    aliases: StringInterner,
    sources: StringInterner,

    mapping: HashMap<DefaultSymbol, Vec<Entry>>,
}

#[derive(Copy, Clone)]
struct ProteinStringId<S = DefaultSymbol> {
    string_id: S,
    source: S,
}

impl ProteinStringId<DefaultSymbol> {
    pub fn resolve(self, interners: &StringProteinAliases<Self>) -> Option<ProteinStringId<&str>> {
        Some(ProteinStringId {
            string_id: interners.string_ids.resolve(self.string_id)?,
            source: interners.sources.resolve(self.source)?
        })
    }
}

#[derive(Copy, Clone)]
struct ProteinAlias<S = DefaultSymbol> {
    alias: S,
    source: S,
}

impl ProteinAlias<DefaultSymbol> {
    pub fn resolve(self, interners: &StringProteinAliases<Self>) -> Option<ProteinAlias<&str>> {
        Some(ProteinAlias {
            alias: interners.aliases.resolve(self.alias)?,
            source: interners.sources.resolve(self.source)?
        })
    }
}

impl StringProteinAliases<ProteinStringId> {
    pub fn load_csv<R: io::Read>(mut reader: csv::Reader<R>) -> anyhow::Result<Self> {
        let mut result = StringProteinAliases {
            aliases: StringInterner::new(),
            string_ids: StringInterner::new(),
            sources: StringInterner::new(),
            mapping: HashMap::new(),
        };

        let mut current = None;

        for record in reader.records() {
            let record = record?;

            assert!(record.len() == 3);

            let string_id = result.string_ids.get_or_intern(&record[0]);
            let alias = result.aliases.get_or_intern(&record[1]);

            let entry = match current {
                Some((current_alias, current_ids)) if current_alias == alias => current_ids,
                _ => result.mapping.entry(alias).or_default(),
            };

            for source in record[2].split_whitespace() {
                let source = result.sources.get_or_intern(source);
                entry.push(ProteinStringId { string_id, source });
            }

            current = Some((alias, entry));
        }

        Ok(result)
    }

    pub fn find_string_ids(&self, alias: &str) -> Option<&[ProteinStringId]> {
        let alias = self.aliases.get(alias)?;
        let entry = self.mapping.get(&alias)?;
        Some(entry)
    }
}

fn read_string_protein_aliases(
    path: &Path,
) -> anyhow::Result<StringProteinAliases<ProteinStringId>> {
    let gz_decoder = GzDecoder::new(BufReader::new(File::open(path)?));

    let csv_reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(gz_decoder);

    StringProteinAliases::load_csv(csv_reader)
}

fn read_uniprot_blastout(
    path: impl AsRef<Path>,
    blast_outfmt: &BlastOutFmt,
    filters: &[BlastOutFilter],
) -> anyhow::Result<LazyFrame> {
    use blastout::fields;

    let mut wanted_columns = Vec::from([fields::QSEQID.0, fields::SACCVER.0, fields::STAXID.0]);

    for filter in filters {
        let column = filter.get_column();

        if !wanted_columns.contains(&column) {
            wanted_columns.push(column);
        }
    }

    let df = blastout::load_blastout(path, blast_outfmt, Some(wanted_columns))?;

    Ok(blastout::apply_filters(
        df,
        filters.iter().cloned().collect(),
    ))
}

#[derive(Args)]
pub struct StringVirusInteractions {
    #[clap(flatten)]
    taxonomy: PreprocessedTaxonomyArgs,

    /// TANGO3 metagenomic assignment output
    metagenomic_assignment: String,

    /// Does the presence of a higher-ranked node in the metagenomic assignment imply the presence
    /// of any of its descendants?
    #[clap(long)]
    include_descendants: bool,

    /// Path to the viral UniProt BLAST search. At least, the fields qseqid and saccver must be
    /// present.
    #[clap(long)]
    blastout_path: String,

    /// BLAST+ output format. Only 6 and 7 are supported at the moment. Columns can be specified
    /// just like in BLAST+, like '7 qseqid saccver pident bitscore'. Irrelevant columns can be
    /// ignored using _ instead, e.g. '7 qseqid saccver _ bitscore'.
    #[clap(long)]
    blast_outfmt: String,

    /// Filters to apply to the viral UniProt BLAST search. Example: --filter 'pident>=90'
    #[clap(long)]
    filter: Vec<String>,

    /// Directory for STRING-Viruses downloads.
    #[clap(long)]
    stringdb_dir: String,

    /// Output file
    #[clap(short, long)]
    output: String,
}

fn find_stringdb_interactions(args: StringVirusInteractions) -> anyhow::Result<()> {
    with_some_ncbi_or_newick_taxonomy!(args.taxonomy.deserialize()?.tree,
        ncbi: taxonomy => {
            let postorder_taxids: Vec<_> = taxonomy.postorder_descendants(taxonomy.get_root()).collect();

            let taxids: HashSet<_> = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(true)
                .from_path(&args.metagenomic_assignment)?
                .into_deserialize()
                .map(|record| record.map(|a: AssignmentRecord| a.assigned_taxid))
                .collect::<Result<_, _>>()?;

            let blast_outfmt = args.blast_outfmt.parse()?;
            let blast_filters = BlastOutFilter::parse_filters(&args.filter)?;
            let string_virus_hits = read_uniprot_blastout(&args.blastout_path, &blast_outfmt, &blast_filters)?;

            //let viral_proteins_uniprot = read_uniprot_blastout(
                //viral_proteins_uniprot_blastout,
                //viral_proteins_uniprot_blast_outfmt,
                //viral_proteins_uniprot_blast_filter,
            //)
            //.context("Reading BLAST search of viral proteins against UniProt")?;
        },
        newick: _ => {
            println!("Newick taxonomies are not supported in this utility");
        },
    );

    Ok(())
}
