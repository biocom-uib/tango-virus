use std::{path::Path, fs::{self, File}, time::Duration};

use anyhow::anyhow;
use clap::Args;
use flate2::read::GzDecoder;
use ftp::FtpStream;

use crate::util::progress_monitor;

const FTP_SERVER: &str = "ftp.ncbi.nlm.nih.gov:21";
const REMOTE_TAXONOMY_DIRECTORY: &str = "/pub/taxonomy";
const TAXDUMP_FILE_NAME: &str = "taxdump.tar.gz";

/// Fetch the NCBI Taxonomy
#[derive(Args)]
pub struct NcbiTaxonomyFetchArgs {
    /// Always re-download existing files
    #[clap(long, short)]
    pub force_download: bool,

    /// Directory to download the NCBI Taxonomy into
    #[clap(long, short = 'd')]
    pub download_dir: String,

    /// Do not unpack taxdump
    #[clap(long)]
    pub no_extract: bool,
}

fn ftp_connect() -> anyhow::Result<FtpStream> {
    let mut conn = FtpStream::connect(FTP_SERVER)?;
    conn.login("anonymous", "anonymous")?;
    conn.cwd(REMOTE_TAXONOMY_DIRECTORY)?;
    Ok(conn)
}

pub fn fetch(args: &NcbiTaxonomyFetchArgs) -> anyhow::Result<()> {
    let work_dir = Path::new(&args.download_dir);

    fs::create_dir_all(work_dir)?;

    let mut ftp_conn = None;

    let file_name = TAXDUMP_FILE_NAME;
    let remote_path = file_name;
    let local_partial_path = work_dir.join(format!("{file_name}.part"));
    let local_path = work_dir.join(file_name);

    if args.force_download || !local_path.exists() {
        let conn = if let Some(conn) = &mut ftp_conn {
            conn
        } else {
            ftp_conn.insert(ftp_connect()?)
        };

        println!("Downloading {remote_path:?} to {local_path:?}");

        let file_size = conn.size(file_name)?.ok_or(anyhow!(
            "File {file_name:?} does not exist in the FTP server"
        ))?;

        let r = conn.retr(file_name, |reader| {
            let io_res = File::create(&local_partial_path).and_then(|file_writer| {
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
        });

        if matches!(&r, Err(_) | Ok(Err(_))) {
            fs::remove_file(&local_partial_path)?;
        }

        r??;

        fs::rename(&local_partial_path, &local_path)?;

        if !args.no_extract {
            let decoder = GzDecoder::new(File::open(&local_path)?);
            let mut archive = tar::Archive::new(decoder);
            archive.unpack(work_dir)?;

            println!("{file_name} extracted into {work_dir:?}");
        }
    } else {
        println!("File {file_name:?} is already present as {local_path:?}, skipping (delete it or use --force-download to re-download)");
    }

    if let Some(conn) = &mut ftp_conn {
        conn.quit()?;
    }

    Ok(())
}
