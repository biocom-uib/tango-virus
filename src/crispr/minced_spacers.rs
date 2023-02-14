use std::{path::{Path, PathBuf}, env, process};

const MINCED_JAR_ENV_VAR: &str = "METEOR_MINCED";

const CRISPRS_FILE_NAME: &str = "crisprs.txt";
const SPACERS_FILE_NAME: &str = "crisprs_spacers.fa"; // hardcoded in minced since v0.1.5

const BLASTOUT_FILE_NAME: &str = "spacer_search.blastout";

struct MincedArgs {
    work_dir: PathBuf,
    minced_jar: PathBuf,
}

impl MincedArgs {
    fn new(work_dir: &Path, minced_jar: Option<&Path>) -> Option<Self> {
        let minced_jar = minced_jar.map(PathBuf::from).or_else(|| {
            let path_str = env::var(MINCED_JAR_ENV_VAR).ok()?;
            Some(PathBuf::from(path_str))
        });

        Some(MincedArgs {
            work_dir: PathBuf::from(work_dir),
            minced_jar: minced_jar?,
        })
    }

    fn minced_cmd(&self) -> process::Command {
        let mut cmd = process::Command::new("java");

        cmd.args(["-jar", &self.work_dir.to_string_lossy()]);
        cmd.current_dir(&self.work_dir);

        cmd
    }

    fn blast_cmd(&self, perc_identity: i32, subject: &Path, output: &Path) -> process::Command {
        let mut cmd = process::Command::new("blastn");

        let query_path = self.work_dir.join(SPACERS_FILE_NAME);
        let blastout_path = self.work_dir.join(BLASTOUT_FILE_NAME);

        cmd
            .args(["-task", "blastn-short"])
            .args(["-perc_identity", &perc_identity.to_string()])
            .args(["-query", &query_path.to_string_lossy()])
            .args(["-subject", &subject.to_string_lossy()])
            .args(["-out", &blastout_path.to_string_lossy()]);

        cmd
    }
}

fn produce_spacers(args: &MincedArgs, metagenomic_seqs: &Path) {
}
