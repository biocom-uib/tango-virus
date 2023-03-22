use std::{fs::{self, File}, path::{Path, PathBuf}, io::{BufReader, Write}, sync::mpsc::{Receiver, SyncSender}, process, time::Duration};

use anyhow::Context;
use clap::{Args, ArgGroup};
use flate2::read::GzDecoder;
use tempdir::TempDir;

use crate::{util::{progress_monitor, cli_tools::{BlastTool, CliTool}}, fetch};

use super::uniprot_xml_parser::{UniProtXmlReader, EntryBuilder, dbref_types, SubfieldsAction};

/// Prepare a BLAST database using UniProt protein sequences
#[derive(Args, Clone)]
#[clap(group(ArgGroup::new("makebastdb").multiple(false)))]
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

    /// Directory with UniProt files. See `meteor fetch uniprot --help`.
    #[clap(long)]
    uniprot_dir: String,

    /// Directory to use to prepare the inputs for makeblastdb. A temporary directory is used by
    /// default.
    #[clap(long)]
    work_dir: Option<String>,

    /// Do not run makeblastdb, just set up input files
    #[clap(long, requires = "work_dir", group = "makeblastdb")]
    skip_makeblastdb: bool,

    /// Installation prefix of the NCBI BLAST+ toolkit to use.
    #[clap(long)]
    blast_prefix: Option<String>,

    /// Output BLAST DB path. Default: <WORK_DIR>/blastdb
    #[clap(short, long, group = "makeblastdb")]
    output: Option<String>,
}

impl From<PrepareUniProtBlastDBArgs> for fetch::uniprot::UniProtFetchArgs {
    fn from(value: PrepareUniProtBlastDBArgs) -> Self {
        fetch::uniprot::UniProtFetchArgs {
            download_dir: value.uniprot_dir,
            uniprot_divisions: value.uniprot_divisions,
            exclude_swissprot: value.exclude_swissprot,
            exclude_trembl: value.exclude_trembl,
            force_download: false,
        }
    }
}

#[derive(Clone, Debug, Default)]
struct UniProtEntryBuilder {
    primary_accession: Option<String>,
    taxid: Option<usize>,
    sequence: Option<String>,
}

#[derive(Debug)]
struct UniProtEntry {
    primary_accession: String,
    taxid: usize,
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
    smap: SyncSender<(String, usize)>,
) -> anyhow::Result<()> {
    struct UniProtEntryReader<F: FnMut(UniProtEntry) -> anyhow::Result<()>>(F);

    impl<F> progress_monitor::ReaderConsumer for UniProtEntryReader<F>
    where
        F: FnMut(UniProtEntry) -> anyhow::Result<()>,
    {
        type Output = anyhow::Result<()>;

        fn read_from<R: std::io::Read>(mut self, reader: R) -> Self::Output {
            let decoder = GzDecoder::new(reader);

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

fn write_taxid_mapping(taxid_map_path: &Path, rmap: Receiver<(String, usize)>) -> anyhow::Result<()> {
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

fn prepare_files_for_makeblastdb(work_dir: &Path, uniprot_files: &[PathBuf]) -> anyhow::Result<()> {
    let fasta_path = work_dir.join(FASTA_FILE_NAME);
    let taxid_map_path = work_dir.join(TAXID_MAP_FILE_NAME);

    let (sseq, rseq) = std::sync::mpsc::sync_channel(100);
    let (smap, rmap) = std::sync::mpsc::sync_channel(1000);

    let (processing_errors, (fasta_errors, taxid_map_errors)) = rayon::join(
        move || read_all_entries(&uniprot_files, sseq, smap),
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

/*pub fn build_sample_blastp_command(
    args: &PrepareUniProtBlastDBArgs,
    blast_prefix: Option<&Path>,
    db_path: &Path,
) -> anyhow::Result<process::Command> {
    let (found, mut cmd) = if let Some(blastp_tool) = BlastTool::resolve(blast_prefix, "blastp")? {
        (true, blastp_tool.new_command())
    } else {
        (false, process::Command::new("blastp"))
    };

    cmd.arg("-db")
        .arg(db_path)
        .arg("-query")
        .arg("PATH_TO_VIRAL_PROTEINS_FASTA")
        .arg("-outfmt")
        .arg("7 qseqid sseqid staxid evalue pident")
        .arg("-out")
        .arg("PATH_TO_BLAST_OUTPUT");

    Ok(cmd)
}*/

pub fn prepare_uniprot_blastdb(args: PrepareUniProtBlastDBArgs) -> anyhow::Result<()> {
    let blast_prefix = args.blast_prefix.as_ref().map(|s| s.as_ref());

    let mut makeblastdb_command = if args.skip_makeblastdb {
        process::Command::new("makeblastdb") // for display purposes only
    } else if let Some(tool) = BlastTool::resolve(blast_prefix, "makeblastdb")? {
        tool.new_command()
    } else {
        anyhow::bail!(
            "makeblastdb executable not found, please make sure that BLAST+ is installed"
        );
    };

    let mut temp_dir = None;

    let work_dir = if let Some(work_dir) = &args.work_dir {
        Path::new(work_dir)
    } else {
        temp_dir.insert(TempDir::new("prepare-uniprot-blastdb")?).path()
    };

    let fasta_path = work_dir.join(FASTA_FILE_NAME);
    let taxid_map_path = work_dir.join(TAXID_MAP_FILE_NAME);

    if !fasta_path.exists() || !taxid_map_path.exists() {
        println!("Preparing makeblastdb input files from UniProt divisions");

        let uniprot_files = fetch::uniprot::list_files(&args.clone().into())
            .context("Searching UniProt division files.")?;

        prepare_files_for_makeblastdb(work_dir, &uniprot_files)
            .context("Preparing files for makeblastdb")?;

        println!("TAXID mapping file written to {taxid_map_path:?}");
        println!("Sequences written to {fasta_path:?}");
    } else {
        println!("makeblastdb input files already exist, proceeding to makeblastdb");
    }

    let db_path = args
        .output
        .as_ref()
        .map(PathBuf::from)
        .unwrap_or(work_dir.join("blastdb"));

    if let Some(parent) = db_path.parent() {
        fs::create_dir_all(parent)?;
    }

    makeblastdb_command
        .arg("-parse_seqids")
        .args(["-input_type", "fasta"])
        .args(["-dbtype", "prot"])
        .arg("-in")
        .arg(&fasta_path)
        .arg("-taxid_map")
        .arg(&taxid_map_path)
        .arg("-out")
        .arg(&db_path)
        .stdin(process::Stdio::null());

    if args.skip_makeblastdb {
        println!("Recommended makeblastdb command:");
        println!("  $ {makeblastdb_command:?}");

        //println!("Next, the recommended blastp command is:");
        //println!("  $ {:?}", build_sample_blastp_command(&args, blast_prefix, &db_path)?);

    } else {
        println!("Running {makeblastdb_command:?}");

        let makeblastdb_status = makeblastdb_command
            .status()
            .context("Executing makeblastdb")?;

        if makeblastdb_status.success() {
            println!("makeblastdb finished successfully");

            //println!("Next, the recommended blastp command is:");
            //println!("  $ {:?}", build_sample_blastp_command(&args, blast_prefix, &db_path)?);

        } else {
            if let Some(code) = makeblastdb_status.code() {
                println!("makeblastdb finished with exit code {code}");
            } else {
                println!("Unknown error running makeblastdb");
            }

            if let Some(temp_dir) = temp_dir.take() {
                let path = temp_dir.into_path();
                println!("Intentionally not deleting temporary work directory for inspection due to error: {path:?}");
            }
        }
    }

    Ok(())
}
