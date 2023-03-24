use std::{path::{PathBuf, Path}, fs::{self, File}, time::Duration};

use anyhow::anyhow;
use clap::Args;
use ftp::FtpStream;

use crate::util::progress_monitor;

const FTP_SERVER: &str = "ftp.uniprot.org:21";
const REMOTE_DIVISIONS_DIRECTORY: &str = "/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions";
const SWISSPROT_FILE_NAME_PART: &str = "sprot";
const TREMBL_FILE_NAME_PART: &str = "trembl";

macro_rules! format_uniprot_division_file_name {
    (db: $db:ident, division: $division:ident) => {
        format!("uniprot_{}_{}.xml.gz", $db, $division)
    }
}

/// Fetch UniProt division sequences and metadata.
#[derive(Args)]
pub struct UniProtFetchArgs {
    /// UniProt divisions to include. Example: viruses. See https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/README
    #[clap(num_args = 1..)]
    pub uniprot_divisions: Vec<String>,

    /// Do not include SwissProt proteins
    #[clap(long)]
    pub exclude_swissprot: bool,

    /// Do not include TrEMBL proteins
    #[clap(long)]
    pub exclude_trembl: bool,

    /// Always re-download existing files
    #[clap(long)]
    pub force_download: bool,

    /// Directory to download UniProt files into
    #[clap(long, short = 'd')]
    pub download_dir: String,
}

impl UniProtFetchArgs {
    fn db_filename_parts(&self) -> Vec<&str> {
        let mut dbs = Vec::new();

        if !self.exclude_swissprot {
            dbs.push(SWISSPROT_FILE_NAME_PART);
        }

        if !self.exclude_trembl {
            dbs.push(TREMBL_FILE_NAME_PART);
        }

        dbs
    }
}

fn ftp_connect() -> anyhow::Result<FtpStream> {
    let mut conn = FtpStream::connect(FTP_SERVER)?;
    conn.login("anonymous", "anonymous")?;
    conn.cwd(REMOTE_DIVISIONS_DIRECTORY)?;
    Ok(conn)
}

pub fn fetch(args: &UniProtFetchArgs) -> anyhow::Result<Vec<PathBuf>> {
    let work_dir = Path::new(&args.download_dir);

    fs::create_dir_all(work_dir)?;

    let dbs = args.db_filename_parts();

    let mut ftp_conn = None;

    let mut paths = Vec::new();

    for db in dbs {
        for division in &args.uniprot_divisions {
            let file_name = format_uniprot_division_file_name!(db: db, division: division);
            let remote_path = format!("{REMOTE_DIVISIONS_DIRECTORY}/{file_name}");
            let local_partial_path = work_dir.join(format!("{file_name}.part"));
            let local_path = work_dir.join(&file_name);

            if args.force_download || !local_path.exists() {
                let conn = if let Some(conn) = &mut ftp_conn {
                    conn
                } else {
                    ftp_conn.insert(ftp_connect()?)
                };

                eprintln!("Downloading {remote_path:?} to {local_path:?}");

                let file_size = conn.size(&file_name)?.ok_or(anyhow!(
                    "File {file_name:?} does not exist in the FTP server"
                ))?;

                let r = conn.retr(&file_name, |reader| {
                    let io_res = File::create(&local_partial_path).and_then(|file_writer| {
                        progress_monitor::monitor_copy(
                            reader,
                            file_writer,
                            Duration::from_secs(1),
                            progress_monitor::default_progress_callback(file_size as u64),
                            progress_monitor::default_finish_callback(),
                        )
                    });

                    eprintln!();

                    Ok(io_res)
                });

                if matches!(&r, Err(_) | Ok(Err(_))) {
                    fs::remove_file(&local_partial_path)?;
                }

                r??;

                fs::rename(&local_partial_path, &local_path)?;
            } else {
                println!("File {file_name:?} is already present as {local_path:?}, skipping (delete it or use --force-download to re-download)");
            }

            paths.push(local_path);
        }
    }

    if let Some(conn) = &mut ftp_conn {
        conn.quit()?;
    }

    Ok(paths)
}

pub fn list_files(args: &UniProtFetchArgs) -> anyhow::Result<Vec<PathBuf>> {
    let mut paths = Vec::new();

    for db in args.db_filename_parts() {
        for division in &args.uniprot_divisions {
            let file_name = format_uniprot_division_file_name!(db: db, division: division);

            let path = Path::new(&args.download_dir).join(&file_name);

            if path.exists() {
                paths.push(path);
            } else {
                anyhow::bail!("{} does not exist in {}. Run `meteor fetch uniprot` first?", &file_name, &args.download_dir);
            }

        }
    }

    Ok(paths)
}
