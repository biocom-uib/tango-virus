use std::{fs::{self, File}, path::{Path, PathBuf}, io::{self, BufReader, Write, ErrorKind}, sync::mpsc::Receiver, process};

use anyhow::{Context, anyhow};
use clap::Args;
use flate2::bufread::GzDecoder;
use ftp::FtpStream;
use itertools::{Itertools, process_results};

use super::uniprot_xml_parser::{UniProtXmlReader, EntryBuilder, dbref_types};

/// Prepare a BLAST database using UniProt protein sequences
#[derive(Args)]
pub struct PrepareUniProtBlastDBArgs {
    /// UniProt divisions to include. Example: viruses. See https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/README
    #[clap(num_args = 1..)]
    uniprot_divisions: Vec<String>,

    /// Do not include SwissProt proteins
    #[clap(long)]
    skip_swissprot: bool,

    /// Do not include TrEMBL proteins
    #[clap(long)]
    skip_trembl: bool,

    /// Always re-download exiting files
    #[clap(long)]
    force_download: bool,

    /// Directory to use to download UniProt files and prepare the inputs for makeblastdb.
    #[clap(long, short = 'd')]
    work_dir: String,

    /// Output BLAST DB path
    #[clap(short, long)]
    output: String,
}


const FTP_SERVER: &str = "ftp.uniprot.org:21";
const REMOTE_DIVISIONS_DIRECTORY: &str = "/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions";
const SWISSPROT_FILE_NAME_PART: &str = "sprot";
const TREMBL_FILE_NAME_PART: &str = "trembl";

macro_rules! format_uniprot_division_file_name {
    (db: $db:ident, division: $division:ident) => {
        format!("uniprot_{}_{}.xml.gz", $db, $division)
    }
}


fn ensure_makeblastdb_is_executable() -> anyhow::Result<bool> {
    let result = process::Command::new("makeblastdb")
        .arg("-version")
        .stdin(process::Stdio::null())
        .stdout(process::Stdio::null())
        .stderr(process::Stdio::null())
        .status();

    match &result {
        Err(e) if e.kind() == ErrorKind::NotFound => Ok(false),
        _ => Ok(result?.success()),
    }
}

fn ftp_connect() -> anyhow::Result<FtpStream> {
    let mut conn = FtpStream::connect(FTP_SERVER)?;
    conn.login("anonymous", "anonymous")?;
    conn.cwd(REMOTE_DIVISIONS_DIRECTORY)?;
    Ok(conn)
}

fn fetch(args: &PrepareUniProtBlastDBArgs) -> anyhow::Result<Vec<PathBuf>> {
    let work_dir = Path::new(&args.work_dir);

    fs::create_dir_all(work_dir)?;

    let dbs = {
        let mut dbs = Vec::new();

        if !args.skip_swissprot {
            dbs.push(SWISSPROT_FILE_NAME_PART);
        }

        if !args.skip_trembl {
            dbs.push(TREMBL_FILE_NAME_PART);
        }

        dbs
    };

    let mut ftp_conn = None;

    let mut paths = Vec::new();

    for db in dbs {
        for division in &args.uniprot_divisions {
            let file_name = format_uniprot_division_file_name!(db: db, division: division);
            let file_path = work_dir.join(&file_name);

            if args.force_download || !file_path.exists() {
                let conn = if let Some(conn) = &mut ftp_conn {
                    conn
                } else {
                    ftp_conn.insert(ftp_connect()?)
                };

                println!("Downloading {file_name:?} to {file_path:?}");

                conn.retr(&file_name, |reader| {
                    let r =
                        File::create(&file_path).and_then(|mut file| io::copy(reader, &mut file));

                    Ok(r)
                })??;

                paths.push(file_path);
            }
        }
    }

    if let Some(conn) = &mut ftp_conn {
        conn.quit()?;
    }

    Ok(paths)
}

fn load_compressed_uniprot_division(path: &Path) -> anyhow::Result<UniProtXmlReader<BufReader<GzDecoder<BufReader<File>>>>> {
    let file = File::open(path)?;
    let decoder = GzDecoder::new(BufReader::new(file));
    let reader = UniProtXmlReader::new(BufReader::new(decoder));

    Ok(reader)
}

//struct LazyIntoIterator<T, F: FnOnce() -> T>(F);

//impl<T: IntoIterator, F: FnOnce() -> T> IntoIterator for LazyIntoIterator<T, F> {
//    type Item = T::Item;
//    type IntoIter = T::IntoIter;

//    fn into_iter(self) -> Self::IntoIter {
//        self.0().into_iter()
//    }
//}

#[derive(Default)]
struct UniProtEntryBuilder {
    primary_accession: Option<String>,
    taxid: Option<i64>,
    sequence: Option<String>,
}

struct UniProtEntry {
    primary_accession: String,
    taxid: i64,
    sequence: String,
}

impl EntryBuilder for UniProtEntryBuilder {
    type Entry = UniProtEntry;

    fn accession(&mut self, acc: &str) {
        self.primary_accession.get_or_insert_with(|| acc.to_owned());
    }

    fn name(&mut self, _name: &str) {}

    fn sequence(&mut self, sequence: &str) {
        self.sequence.get_or_insert_with(|| sequence.to_owned());
    }

    fn begin_dbreference(&mut self, _ref_type: &str, _ref_id: &str) {}
    fn dbreference_property(&mut self, _prop_type: &str, _prop_value: &str) {}
    fn finish_dbreference(&mut self) {}

    fn begin_organism(&mut self) {}
    fn organism_scientific_name(&mut self, _name: &str) {}
    fn organism_common_name(&mut self, _name: &str) {}

    fn begin_organism_dbreference(&mut self, ref_type: &str, ref_id: &str) {
        if self.taxid.is_none() && ref_type == dbref_types::NCBI_TAXONOMY {
            self.taxid = ref_id.parse().ok();
        }
    }

    fn organism_dbreference_property(&mut self, _prop_type: &str, _prop_value: &str) {}
    fn finish_organism_dbreference(&mut self) {}
    fn finish_organism(&mut self) {}

    fn finish(self) -> Option<Self::Entry> {
        Some(UniProtEntry {
            primary_accession: self.primary_accession?,
            taxid: self.taxid?,
            sequence: self.sequence?,
        })
    }
}


const FASTA_FILE_NAME: &str = "sequences.fa";
const TAXID_MAP_FILE_NAME: &str = "taxid.tab";

fn write_sequences(fasta_path: &Path, rseq: Receiver<(String, String)>) -> anyhow::Result<()> {
    let mut fasta_file = File::create(fasta_path)?;

    for (primary_accession, sequence) in rseq {
        writeln!(&mut fasta_file, ">{primary_accession}")?;
        writeln!(&mut fasta_file, "{sequence}")?;
    }

    Ok(())
}

fn write_taxid_mapping(taxid_map_path: &Path, rmap: Receiver<(String, i64)>) -> anyhow::Result<()> {
    let mut taxid_map_writer = csv::WriterBuilder::new()
        .has_headers(false)
        .delimiter(b' ')
        .from_path(taxid_map_path)?;

    for (primary_accession, taxid) in rmap {
        taxid_map_writer.serialize((&primary_accession, taxid))?;
    }

    taxid_map_writer.flush()?;
    Ok(())
}

fn prepare_files_for_makeblastdb(args: &PrepareUniProtBlastDBArgs) -> anyhow::Result<()> {
    let paths = fetch(args).context("Downloading UniProt division files")?;

    let work_dir = Path::new(&args.work_dir);
    let fasta_path = work_dir.join(FASTA_FILE_NAME);
    let taxid_map_path = work_dir.join(TAXID_MAP_FILE_NAME);

    let all_entries = paths
        .iter()
        .map(|path| {
            let reader = load_compressed_uniprot_division(path).context("Loading {path:?}")?;
            Ok(reader
                .into_entries_iter::<UniProtEntryBuilder>()
                .map(|entry| entry.context("Processing entries in {path:?}")))
        })
        .flatten_ok()
        .map(|entry| entry.flatten());

    let (sseq, rseq) = std::sync::mpsc::sync_channel(100);
    let (smap, rmap) = std::sync::mpsc::sync_channel(1000);

    let (processing_errors, (fasta_errors, taxid_map_errors)) = rayon::join(
        move || {
            process_results(all_entries, |all_entries| {
                for entry in all_entries {
                    sseq.send((entry.primary_accession.clone(), entry.sequence))?;
                    smap.send((entry.primary_accession, entry.taxid))?;
                }
                Ok(())
            })
            .and_then(|r| r)
        },
        move || {
            rayon::join(
                move || write_sequences(&fasta_path, rseq),
                move || write_taxid_mapping(&taxid_map_path, rmap),
            )
        },
    );

    processing_errors.context("Processing XML entries")?;
    fasta_errors.context("Writing the FASTA file")?;
    taxid_map_errors.context("Writing the TAXID mapping file")?;

    Ok(())
}

pub fn prepare_uniprot_blastdb(args: PrepareUniProtBlastDBArgs) -> anyhow::Result<()> {
    if !ensure_makeblastdb_is_executable()? {
        anyhow::bail!("makeblastdb executable not found, please make sure that BLAST+ is installed");
    }

    prepare_files_for_makeblastdb(&args)
        .context("Preparing files for makeblastdb")?;

    let work_dir = Path::new(&args.work_dir);
    let fasta_path = work_dir.join(FASTA_FILE_NAME);
    let taxid_map_path = work_dir.join(TAXID_MAP_FILE_NAME);

    let db_path = Path::new(&args.output);

    if let Some(parent) = db_path.parent() {
        fs::create_dir_all(parent)?;
    }

    let fasta_path = fasta_path.to_str().ok_or(anyhow!("FASTA path contains invalid unicode"))?;
    let taxid_map_path = taxid_map_path.to_str().ok_or(anyhow!("TAXID mapping path contains invalid unicode"))?;

    let makeblastdb_status = process::Command::new("makeblastdb")
        .arg("-parse_seqids")
        .args(["-input_type", "fasta"])
        .args(["-in", fasta_path])
        .args(["-taxid_map", taxid_map_path])
        .args(["-out", &args.output])
        .stdin(process::Stdio::null())
        .status()
        .context("Executing makeblastdb")?;

    if makeblastdb_status.success() {
        println!("makeblastdb finished successfully");
    } else if let Some(code) = makeblastdb_status.code() {
        println!("makeblastdb finished with exit code {code}");
    } else {
        println!("makeblastdb was interrupted");
    }

    Ok(())
}
