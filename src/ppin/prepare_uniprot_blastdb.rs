use std::{fs::{self, File}, path::{Path, PathBuf}, io::{BufReader, Write, ErrorKind}, sync::mpsc::{Receiver, SyncSender}, process, time::Duration};

use anyhow::{Context, anyhow};
use clap::{Args, ArgGroup};
use flate2::bufread::GzDecoder;
use ftp::FtpStream;

use crate::util::progress_monitor;

use super::uniprot_xml_parser::{UniProtXmlReader, EntryBuilder, dbref_types, SubfieldsAction, IdentityBuilder};

/// Prepare a BLAST database using UniProt protein sequences
#[derive(Args)]
#[clap(group(
        ArgGroup::new("makebastdb")
            .required(true)
            .multiple(true)
            .args(&["skip_makeblastdb", "output"])
))]
pub struct PrepareUniProtBlastDBArgs {
    /// UniProt divisions to include. Example: viruses. See https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/README
    #[clap(num_args = 1..)]
    uniprot_divisions: Vec<String>,

    /// Do not include SwissProt proteins
    #[clap(long)]
    exclude_swissprot: bool,

    /// Do not include TrEMBL proteins
    #[clap(long)]
    exclude_trembl: bool,

    /// Always re-download existing files
    #[clap(long)]
    force_download: bool,

    /// Directory to use to download UniProt files and prepare the inputs for makeblastdb.
    #[clap(long, short = 'd')]
    work_dir: String,

    /// Do not run makeblastdb, just set up input files
    #[clap(long)]
    skip_makeblastdb: bool,

    /// Output BLAST DB path. Default: <WORK_DIR>/blastdb
    #[clap(short, long)]
    output: Option<String>,
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

        if !args.exclude_swissprot {
            dbs.push(SWISSPROT_FILE_NAME_PART);
        }

        if !args.exclude_trembl {
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

                let file_size = conn.size(&file_name)?.ok_or(anyhow!(
                    "File {file_name:?} does not exist in the FTP server"
                ))?;

                conn.retr(&file_name, |reader| {
                    let io_res = File::create(&file_path).and_then(|file_writer| {
                        progress_monitor::monitor_copy(
                            reader,
                            file_writer,
                            Duration::from_secs(1),
                            progress_monitor::default_progress_callback(file_size as u64),
                            progress_monitor::default_finish_callback(),
                        )
                    });

                    println!();

                    Ok(io_res)
                })??;
            } else {
                println!("File {file_name:?} is already present as {file_path:?}, skipping (delete it or use --force-download to re-download)");
            }

            paths.push(file_path);
        }
    }

    if let Some(conn) = &mut ftp_conn {
        conn.quit()?;
    }

    Ok(paths)
}

#[derive(Clone, Debug, Default)]
struct UniProtEntryBuilder {
    primary_accession: Option<String>,
    taxid: Option<i64>,
    sequence: Option<String>,
}

#[derive(Debug)]
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

    fn sequence(&mut self, sequence: &str) {
        self.sequence.get_or_insert_with(|| sequence.to_owned());
    }

    fn begin_organism(&mut self) -> SubfieldsAction {
        SubfieldsAction::Populate
    }

    fn begin_organism_dbreference(&mut self, ref_type: &str, ref_id: &str) -> SubfieldsAction {
        if self.taxid.is_none() && ref_type == dbref_types::NCBI_TAXONOMY {
            self.taxid = ref_id.parse().ok();
        }

        SubfieldsAction::Skip
    }

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

fn read_all_entries<P: AsRef<Path>>(
    paths: &[P],
    sseq: SyncSender<(String, String)>,
    smap: SyncSender<(String, i64)>,
) -> anyhow::Result<()> {
    struct UniProtEntryReader<F: FnMut(UniProtEntry) -> anyhow::Result<()>>(F);

    impl<F> progress_monitor::ReaderConsumer for UniProtEntryReader<F>
    where
        F: FnMut(UniProtEntry) -> anyhow::Result<()>,
    {
        type Output = anyhow::Result<()>;

        fn read_from<R: std::io::Read>(mut self, reader: R) -> Self::Output {
            let decoder = GzDecoder::new(BufReader::new(reader));

            let uniprot_reader = UniProtXmlReader::new(BufReader::new(decoder));

            for entry in uniprot_reader.into_entries::<UniProtEntryBuilder>() {
                (self.0)(entry?)?;
            }

            Ok(())
        }
    }

    for path in paths {
        let path = path.as_ref();

        let path_str = path.to_string_lossy();
        println!("Loading {}", &path_str);

        let file_size = fs::metadata(path)
            .context(format!("Querying the size of {path_str}"))?
            .len();

        progress_monitor::monitor_reader(
            File::open(path).context(format!("Opening UniProt XML {path:?}"))?,
            Duration::from_secs(1),
            progress_monitor::default_progress_callback(file_size),
            progress_monitor::default_finish_callback(),
            UniProtEntryReader(|entry| {
                sseq.send((entry.primary_accession.clone(), entry.sequence))?;
                smap.send((entry.primary_accession, entry.taxid))?;
                Ok(())
            }),
        )
        .context(format!("Reading UniProt entries from {path_str}"))?;
    }

    anyhow::Ok(())
}

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
    let work_dir = Path::new(&args.work_dir);
    let fasta_path = work_dir.join(FASTA_FILE_NAME);
    let taxid_map_path = work_dir.join(TAXID_MAP_FILE_NAME);

    let paths = fetch(args).context("Downloading UniProt division files")?;

    let (sseq, rseq) = std::sync::mpsc::sync_channel(100);
    let (smap, rmap) = std::sync::mpsc::sync_channel(1000);

    let (processing_errors, (fasta_errors, taxid_map_errors)) = rayon::join(
        move || read_all_entries(&paths, sseq, smap),
        move || {
            rayon::join(
                move || write_sequences(&fasta_path, rseq),
                move || write_taxid_mapping(&taxid_map_path, rmap),
            )
        },
    );

    processing_errors.context("Processing UniProt XML files")?;
    fasta_errors.context("Writing the FASTA file")?;
    taxid_map_errors.context("Writing the TAXID mapping file")?;

    Ok(())
}

pub fn prepare_uniprot_blastdb(args: PrepareUniProtBlastDBArgs) -> anyhow::Result<()> {
    if !args.skip_makeblastdb && !ensure_makeblastdb_is_executable()? {
        anyhow::bail!(
            "makeblastdb executable not found, please make sure that BLAST+ is installed"
        );
    }

    println!("Preparing makeblastdb input files from UniProt divisions");

    prepare_files_for_makeblastdb(&args).context("Preparing files for makeblastdb")?;

    let work_dir = Path::new(&args.work_dir);
    let fasta_path = work_dir.join(FASTA_FILE_NAME);
    let taxid_map_path = work_dir.join(TAXID_MAP_FILE_NAME);

    let fasta_path = fasta_path
        .to_str()
        .ok_or(anyhow!("FASTA path contains invalid unicode"))?;

    let taxid_map_path = taxid_map_path
        .to_str()
        .ok_or(anyhow!("TAXID mapping path contains invalid unicode"))?;

    let db_path = args
        .output
        .as_ref()
        .map(PathBuf::from)
        .unwrap_or(work_dir.join("blastdb"));

    if let Some(parent) = db_path.parent() {
        fs::create_dir_all(parent)?;
    }

    let output = &*db_path.to_string_lossy();

    let mut makeblastdb_command = process::Command::new("makeblastdb");

    makeblastdb_command
        .arg("-parse_seqids")
        .args(["-input_type", "fasta"])
        .args(["-in", fasta_path])
        .args(["-taxid_map", taxid_map_path])
        .args(["-out", output])
        .stdin(process::Stdio::null());

    if args.skip_makeblastdb {
        println!("TAXID mapping file written to {taxid_map_path}");
        println!("Sequences written to {fasta_path}");

        println!("Recommended command: {makeblastdb_command:?}");
    } else {
        println!("Executing: {makeblastdb_command:?}");

        let makeblastdb_status = makeblastdb_command
            .status()
            .context("Executing makeblastdb")?;

        if makeblastdb_status.success() {
            println!("makeblastdb finished successfully");
        } else if let Some(code) = makeblastdb_status.code() {
            println!("makeblastdb finished with exit code {code}");
        } else {
            println!("makeblastdb was interrupted");
        }
    }

    Ok(())
}
