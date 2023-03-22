use std::{path::Path, io, marker::PhantomData, collections::HashSet, fs::File};

use anyhow::Context;
use clap::Args;
use lending_iterator::{HKT, LendingIterator};
use serde::{Serialize, Deserialize};
use string_interner::DefaultSymbol;

use crate::{
    fetch::{self, string_viruses::ProteinLinksRecord},
    preprocessed_taxonomy::{with_some_ncbi_or_newick_taxonomy, PreprocessedTaxonomyArgs},
    tango_assign::AssignmentRecord,
    taxonomy::{
        formats::ncbi::{NamesAssoc, NcbiTaxonomy},
        NodeId, Taxonomy,
    },
    util::{
        csv_stream::CsvReaderIterExt,
        interned_mapping::{InternedMultiMapping, Interner, PairInterner},
        maybe_gzdecoder, writing_new_file_or_stdout,
    },
};

use super::ProteinVirusTaxidMapping;

struct StringProteinAliases<Entry> {
    mapping: InternedMultiMapping<PairInterner>,
    phantom: PhantomData<Entry>,
}

#[derive(Copy, Clone)]
struct ProteinStringId<S = DefaultSymbol> {
    string_id: S,
    source: S,
}

#[derive(Copy, Clone)]
struct ProteinAlias<S = DefaultSymbol> {
    alias: S,
    source: S,
}

impl<Entry> StringProteinAliases<Entry> {
    pub fn get_source_sym(&self, source: &str) -> Option<DefaultSymbol> {
        self.mapping.value_interner().1.get(source)
    }
}

#[allow(dead_code)]
impl StringProteinAliases<ProteinStringId> {
    pub fn read_tsv<R: io::Read>(reader: R) -> anyhow::Result<Self> {
        let mapping = InternedMultiMapping::read_grouped_tsv_with::<SourceRepeaterHKT, _, _>(
            reader,
            false,
            |record| {
                // (alias, (string_id, source))
                Ok((
                    &record[1],
                    SourceRepeater {
                        value: &record[0],
                        source: record[2].split_whitespace(),
                    },
                ))
            },
        )?;

        Ok(StringProteinAliases {
            mapping,
            phantom: PhantomData,
        })
    }

    pub fn find_string_ids(&self, alias: &str) -> impl Iterator<Item = ProteinStringId> + '_ {
        self.mapping
            .lookup_syms(alias)
            .map(|(string_id, source)| ProteinStringId { string_id, source })
    }

    pub fn resolve_string_id(&self, string_id: ProteinStringId) -> Option<ProteinStringId<&str>> {
        self
            .mapping
            .value_interner()
            .resolve((string_id.string_id, string_id.source))
            .map(|(string_id, source)| ProteinStringId { string_id, source })
    }

    pub fn resolve_alias_name(&self, alias: DefaultSymbol) -> Option<&str> {
        self.mapping.key_interner().resolve(alias)
    }

    pub fn resolve_alias(&self, ProteinAlias { alias, source }: ProteinAlias) -> Option<ProteinAlias<&str>> {
        let alias = self.resolve_alias_name(alias)?;

        let source = self
            .mapping
            .value_interner()
            .1
            .resolve(source)?;

        Some(ProteinAlias { alias, source })
    }
}

struct SourceRepeater<'a> {
    value: &'a str,
    source: std::str::SplitWhitespace<'a>,
}

type SourceRepeaterHKT = HKT!(SourceRepeater<'_>);

impl<'a> Iterator for SourceRepeater<'a> {
    type Item = (&'a str, &'a str);

    fn next(&mut self) -> Option<Self::Item> {
        let source = self.source.next()?;
        Some((self.value, source))
    }
}

#[allow(dead_code)]
impl StringProteinAliases<ProteinAlias> {
    pub fn read_tsv<R: io::Read>(reader: R) -> anyhow::Result<Self> {
        let mapping = InternedMultiMapping::read_grouped_tsv_with::<SourceRepeaterHKT, _, _>(
            reader,
            false,
            |record| {
                // (string_id, (alias, source))
                Ok((
                    &record[0],
                    SourceRepeater {
                        value: &record[1],
                        source: record[2].split_whitespace(),
                    },
                ))
            },
        )?;

        Ok(StringProteinAliases {
            mapping,
            phantom: PhantomData,
        })
    }

    pub fn find_aliases(&self, string_id: &str) -> impl Iterator<Item = ProteinAlias> + '_ {
        self.mapping
            .lookup_syms(string_id)
            .map(|(alias, source)| ProteinAlias { alias, source })
    }

    pub fn resolve_alias_name(&self, alias: DefaultSymbol) -> Option<&str> {
        self.mapping.value_interner().0.resolve(alias)
    }

    pub fn resolve_alias(&self, alias: ProteinAlias) -> Option<ProteinAlias<&str>> {
        self
            .mapping
            .value_interner()
            .resolve((alias.alias, alias.source))
            .map(|(alias, source)| ProteinAlias { alias, source })
    }

    pub fn resolve_string_id(&self, ProteinStringId { string_id, source }: ProteinStringId) -> Option<ProteinStringId<&str>> {
        let string_id = self.mapping.key_interner().resolve(string_id)?;

        let source = self
            .mapping
            .value_interner()
            .1
            .resolve(source)?;

        Some(ProteinStringId { string_id, source })
    }
}

#[derive(Args)]
pub struct StringVirusesInteractionsArgs {
    #[clap(flatten)]
    taxonomy: PreprocessedTaxonomyArgs,

    /// TANGO metagenomic assignment output (see `meteor tango-assign --help`)
    metagenomic_assignment: String,

    /// Does the presence of a higher-ranked node in the metagenomic assignment imply the presence
    /// of any of its descendants?
    #[clap(long)]
    include_descendants: bool,

    /// Path to the UniProt - Viral contig name mapping. See `meteor ppi match-viral-proteins`
    /// present.
    #[clap(long)]
    uniprot_virus_mapping: String,

    /// Directory for STRING-Viruses downloads (see `meteor fetch string-viruses --help`)
    #[clap(long)]
    stringdb_dir: String,

    /// Output path. The result is a tab-separated values file containing three columns:
    /// protein_1, protein_2 and combined_score. Use '-' to print to standard output.
    #[clap(short, long, default_value = "-")]
    output: String,
}

fn collect_assignment_taxids<Names: NamesAssoc + 'static>(
    assignment: &Path,
    taxonomy: &NcbiTaxonomy<Names>,
    include_descendants: bool,
) -> anyhow::Result<HashSet<NodeId>> {
    let mut taxids = HashSet::new();

    let mut csv_reader_iter = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(assignment)?
        .into_lending_iter()
        .into_deserialize::<HKT!(AssignmentRecord<&str>)>(None);

    while let Some(record) = csv_reader_iter.next() {
        let taxid = NodeId(record?.assigned_taxid);

        taxids.extend(taxonomy.ancestors(taxid));

        if include_descendants {
            taxids.extend(taxonomy.postorder_descendants(taxid));
        } else {
            taxids.insert(taxid);
        }
    }

    Ok(taxids)
}

fn read_string_uniprot_mapping(stringdb_dir: &Path) -> anyhow::Result<StringProteinAliases<ProteinAlias>> {
    use fetch::string_viruses::FileName;

    maybe_gzdecoder!(FileName::ProteinAliases.open_read(stringdb_dir)?, rdr => {
        StringProteinAliases::<ProteinAlias>::read_tsv(rdr)
    })
}

fn read_uniprot_virus_mapping(path: &Path) -> anyhow::Result<ProteinVirusTaxidMapping> {
    let file = File::open(path)?;

    ProteinVirusTaxidMapping::read_tsv(file)
}

struct VirusAliasHitInfo<'a> {
    accession: &'a str,
    hit_virus_name: &'a str,
}

enum NodeKind<'a, F>
where
    F: FnOnce() -> Option<VirusAliasHitInfo<'a>>
{
    Host(&'a str),
    MaybeVirus(F),
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum EdgeKind {
    VirusVirus,
    VirusHost,
    HostVirus,
    HostHost,
}

#[derive(Serialize)]
pub struct PredictedInteraction<S = String> {
    pub protein_accession_a: S,
    pub protein_accession_b: S,
    pub virus_name: Option<S>,
    pub virus_taxid: Option<usize>,
    pub host_taxid: Option<usize>,
    pub edge_kind: EdgeKind,
}

fn inspect_node_from_string_id<'b, 'a: 'b>(
    assigned_taxids: &HashSet<NodeId>,
    uniprot_virus_mapping: &'a ProteinVirusTaxidMapping,
    string_uniprot_mapping: &'a StringProteinAliases<ProteinAlias>,
    preferred_source_sym: DefaultSymbol,
    string_id: &'b str,
) -> Option<(NodeId, NodeKind<'a, impl FnOnce() -> Option<VirusAliasHitInfo<'a>> + 'b>)> {
    let Some((taxid, _)) = fetch::string_viruses::parse_string_id(string_id) else {
        eprintln!("Warning: Could not parse StringDB identifier: {string_id}");
        return None;
    };

    let mut aliases = string_uniprot_mapping.find_aliases(string_id);

    let result = if assigned_taxids.contains(&NodeId(taxid)) {
        NodeKind::Host(
            string_uniprot_mapping.resolve_alias_name(
                aliases
                .find(|alias| alias.source == preferred_source_sym)?.alias
            )?
        )
    } else {
        NodeKind::MaybeVirus(move || {
            let mut candidate = None;

            for alias in aliases {
                let Some(accession) = string_uniprot_mapping.resolve_alias_name(alias.alias) else {
                    continue;
                };

                let is_preferred_source = alias.source == preferred_source_sym;

                for (hit_virus_name, hit_taxid) in uniprot_virus_mapping.0.lookup(accession) {
                    if hit_taxid != taxid {
                        eprintln!("Warning: Virus TaxID mismatch between StringDB and UniProt BLAST search: {string_id}");
                        continue;
                    }

                    if is_preferred_source {
                        candidate = Some(VirusAliasHitInfo { accession, hit_virus_name });
                        break;
                    } else if candidate.is_none() {
                        candidate = Some(VirusAliasHitInfo { accession, hit_virus_name });
                    }
                }
            }

            candidate
        })
    };

    Some((NodeId(taxid), result))
}

fn process_protein_links_record<'a>(
    assigned_taxids: &HashSet<NodeId>,
    uniprot_virus_mapping: &'a ProteinVirusTaxidMapping,
    string_uniprot_mapping: &'a StringProteinAliases<ProteinAlias>,
    preferred_source_sym: DefaultSymbol,
    record: ProteinLinksRecord<&str>,
) -> Option<PredictedInteraction<&'a str>> {
    let (taxid_a, node_kind_a) = inspect_node_from_string_id(
        assigned_taxids,
        uniprot_virus_mapping,
        string_uniprot_mapping,
        preferred_source_sym,
        record.protein_id_a,
    )?;

    let (taxid_b, node_kind_b) = inspect_node_from_string_id(
        assigned_taxids,
        uniprot_virus_mapping,
        string_uniprot_mapping,
        preferred_source_sym,
        record.protein_id_b,
    )?;

    use NodeKind::*;
    use EdgeKind::*;

    match (node_kind_a, node_kind_b) {
        (Host(acc_a), Host(acc_b)) => {
            if taxid_a != taxid_b {
                eprintln!("Warning: Host-host interaction TaxID mismatch: {} - {}", record.protein_id_a, record.protein_id_b);
                None
            } else {
                Some(PredictedInteraction {
                    protein_accession_a: acc_a,
                    protein_accession_b: acc_b,
                    virus_name: None,
                    virus_taxid: None,
                    host_taxid: Some(taxid_a.0),
                    edge_kind: HostHost,
                })
            }
        }
        (Host(acc_a), MaybeVirus(best_hit_b)) => {
            let best_hit_b = best_hit_b()?;

            Some(PredictedInteraction {
                protein_accession_a: acc_a,
                protein_accession_b: best_hit_b.accession,
                virus_name: Some(best_hit_b.hit_virus_name),
                virus_taxid: Some(taxid_b.0),
                host_taxid: Some(taxid_a.0),
                edge_kind: HostVirus,
            })
        }
        (MaybeVirus(best_hit_a), Host(acc_b)) => {
            let best_hit_a = best_hit_a()?;

            Some(PredictedInteraction {
                protein_accession_a: best_hit_a.accession,
                protein_accession_b: acc_b,
                virus_name: Some(best_hit_a.hit_virus_name),
                virus_taxid: Some(taxid_a.0),
                host_taxid: Some(taxid_b.0),
                edge_kind: VirusHost,
            })
        }
        (MaybeVirus(best_hit_a), MaybeVirus(best_hit_b)) => {
            let best_hit_a = best_hit_a()?;
            let best_hit_b = best_hit_b()?;

            if taxid_a != taxid_b {
                eprintln!("Warning: Virus-virus interaction TaxID mismatch: {} - {}", record.protein_id_a, record.protein_id_b);
                None
            } else {
                Some(PredictedInteraction {
                    protein_accession_a: best_hit_a.accession,
                    protein_accession_b: best_hit_b.accession,
                    virus_name: Some(best_hit_a.hit_virus_name),
                    virus_taxid: Some(taxid_a.0),
                    host_taxid: None,
                    edge_kind: VirusVirus,
                })
            }
        }
    }
}

pub fn string_viruses_interactions(args: StringVirusesInteractionsArgs) -> anyhow::Result<()> {
    let uniprot_virus_mapping = read_uniprot_virus_mapping(args.uniprot_virus_mapping.as_ref())
        .context(format!("Reading the UniProt - Viral contig name mapping from {}", &args.uniprot_virus_mapping))?;

    let assigned_taxids = with_some_ncbi_or_newick_taxonomy!(args.taxonomy.deserialize()?.tree,
        ncbi: taxonomy => {
            collect_assignment_taxids(
                args.metagenomic_assignment.as_ref(),
                &taxonomy,
                args.include_descendants
            )
            .context(format!("Reading metagenomic assigment at {}", &args.metagenomic_assignment))?
        },
        newick: _ => {
            anyhow::bail!("Newick taxonomies are not supported in this utility");
        },
    );

    let stringdb_dir = args.stringdb_dir.as_ref();

    let string_uniprot_mapping = read_string_uniprot_mapping(stringdb_dir)
        .context(format!("Reading the StringDB - UniProt mapping in {}", &args.stringdb_dir))?;

    let preferred_source_sym = string_uniprot_mapping.get_source_sym("UniProtKB-EI")
        .expect("Could not find StringDB alias source: UniProtKB-EI");

    maybe_gzdecoder!(fetch::string_viruses::FileName::ProteinLinks.open_read(stringdb_dir)?, protein_links_reader => {
        writing_new_file_or_stdout!(&args.output, writer => {
            let mut csv_writer = csv::WriterBuilder::new()
                .delimiter(b'\t')
                .has_headers(true)
                .from_writer(writer?);

            let mut records = csv::ReaderBuilder::new()
                .has_headers(false)
                .delimiter(b'\t')
                .from_reader(protein_links_reader)
                .into_lending_iter()
                .into_deserialize::<HKT!(fetch::string_viruses::ProteinLinksRecord<&str>)>(None);

            while let Some(record) = records.next() {
                let record = record?;

                let interaction = process_protein_links_record(&assigned_taxids, &uniprot_virus_mapping, &string_uniprot_mapping, preferred_source_sym, record);

                if let Some(interaction) = interaction {
                    csv_writer.serialize(interaction)?;
                }
            }

            csv_writer.flush()?;
        });
    });

    Ok(())
}
