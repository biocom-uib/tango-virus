use std::{hash::Hash, path::{PathBuf, Path, self}, process::Command, io, marker::PhantomData, fs::File};

use anyhow::Context;
use regex::Regex;
use serde::Deserialize;
use string_interner::StringInterner;

use crate::util::{cli_tools::{self, CliTool}, interned_mapping::{Interner, InternedMultiMapping}};

const FAA_FILE_NAME: &str = "predicted_viral_proteins.faa";

const PROTEIN_VIRUS_NAME_REGEX: &str = r#"_\d+$"#;
const DEFAULT_PRODIGAL_PROCEDURE: &str = "meta";

const DEFAULT_PERC_IDENTITY: i32 = 95;

const BLASTOUT_FMT: &str = "6 qseqid saccver staxid";
const BLASTOUT_FILE_NAME: &str = "uniprot_viral_search.blastout";


#[derive(Deserialize)]
struct ProteinProteinTaxidRecord<'a> {
    predicted_name: &'a str,
    uniprot_name: &'a str,
    taxid: i64,
}

pub struct ProdigalBlastpPipeline {
    work_dir: PathBuf,
    prodigal: cli_tools::Prodigal,
    prodigal_procedure: String,
    blastp: cli_tools::BlastTool,
    blastp_num_threads: Option<i32>,
    dry_run: bool,
}

impl ProdigalBlastpPipeline {
    pub fn new(work_dir: PathBuf, prodigal: cli_tools::Prodigal, blastp: cli_tools::BlastTool) -> Self {
        Self {
            work_dir,
            prodigal,
            prodigal_procedure: DEFAULT_PRODIGAL_PROCEDURE.to_owned(),
            blastp,
            blastp_num_threads: None,
            dry_run: false,
        }
    }

    pub fn set_prodigal_procedure(&mut self, procedure: &str) {
        self.prodigal_procedure = procedure.to_owned();
    }

    pub fn set_blastp_num_threads(&mut self, num_threads: Option<i32>) {
        self.blastp_num_threads = num_threads;
    }

    pub fn set_dry_run(&mut self, dry_run: bool) {
        self.dry_run = dry_run;
    }

    fn prodigal_cmd(&self, viral_seqs: &Path) -> io::Result<Command> {
        let mut cmd = self.prodigal.new_command();

        cmd.current_dir(&self.work_dir)
            .arg("-a")
            .arg(FAA_FILE_NAME)
            .arg("-i")
            .arg(path::absolute(viral_seqs)?)
            .args(["-p", &self.prodigal_procedure]);

        Ok(cmd)
    }

    fn blastp_cmd(&self, uniprot_blastdb: &Path, perc_identity: i32) -> io::Result<Command> {
        let mut cmd = self.blastp.new_command();

        if let Some(num_threads) = self.blastp_num_threads {
            cmd.args(["-num_threads", &num_threads.to_string()]);
        }

        cmd.current_dir(&self.work_dir)
            .arg("-db")
            .arg(path::absolute(uniprot_blastdb)?)
            .arg("-perc_identity")
            .arg(perc_identity.to_string())
            .arg("-query")
            .arg(FAA_FILE_NAME)
            .args(["-outfmt", BLASTOUT_FMT])
            .arg("-out")
            .arg(BLASTOUT_FILE_NAME);

        Ok(cmd)
    }

    pub fn find_and_match_proteins(
        &self,
        viral_seqs: &Path,
        uniprot_blastdb: &Path,
        perc_identity: Option<i32>,
    ) -> anyhow::Result<()> {

        if self.work_dir.join(FAA_FILE_NAME).exists() {
            eprintln!("{FAA_FILE_NAME} already exists in the working directory, skipping MinCED.
                Delete the file to re-generate it");

        } else {
            let mut prodigal_cmd = self.prodigal_cmd(viral_seqs)?;

            eprintln!("Running {prodigal_cmd:?}");

            if !self.dry_run {
                let prodigal_status = prodigal_cmd.status().context("Running Prodigal")?;

                if !prodigal_status.success() {
                    anyhow::bail!("Prodigal finished with non-zero exit code: {}", prodigal_status.code().unwrap_or(0));
                }
            }
        }

        if self.work_dir.join(BLASTOUT_FILE_NAME).exists() {
            eprintln!("{BLASTOUT_FILE_NAME} already exists in the working directory, skipping BLAST
                search. Delete the file to re-generate it");

        } else {
            let perc_identity = perc_identity.unwrap_or(DEFAULT_PERC_IDENTITY);
            let mut blastp_cmd = self.blastp_cmd(uniprot_blastdb, perc_identity)?;

            eprintln!("Running {blastp_cmd:?}");

            if !self.dry_run {
                let blastp_status = blastp_cmd.status().context("Running blastp")?;

                if !blastp_status.success() {
                    anyhow::bail!("blastp finished with non-zero exit code: {}", blastp_status.code().unwrap_or(0));
                }
            }
        }

        Ok(())
    }
}


#[derive(Default)]
struct PrjInterner<T: Default + Copy + Eq + Ord + Hash, I: Interner>(I, PhantomData<T>);

impl<T: Copy + Default + Eq + Ord + Hash, I: Interner> Interner for PrjInterner<T, I> {
    type Value<'a> = (I::Value<'a>, T);

    type Symbol = (I::Symbol, T);

    fn get_or_intern(&mut self, (value, x): Self::Value<'_>) -> Self::Symbol {
        (self.0.get_or_intern(value), x)
    }

    fn resolve(&self, (sym, x): Self::Symbol) -> Option<Self::Value<'_>> {
        Some((self.0.resolve(sym)?, x))
    }
}

pub struct VirusProteinTaxidMapping(InternedMultiMapping::<PrjInterner<i64, StringInterner>>);

impl ProdigalBlastpPipeline {
    pub fn collect_viral_protein_mapping(&self) -> anyhow::Result<VirusProteinTaxidMapping> {
        let protein_virus_name_regex =
            Regex::new(PROTEIN_VIRUS_NAME_REGEX).expect("Error compiling PROTEIN_VIRUS_NAME_REGEX");

        let file_reader = File::open(self.work_dir.join(BLASTOUT_FILE_NAME))?;

        let mapping = InternedMultiMapping::read_tsv_with(file_reader, false, |record| {
            let record = record.deserialize::<ProteinProteinTaxidRecord>(None)?;

            if let Some(suffix_match) = protein_virus_name_regex.find(record.predicted_name) {
                let virus_name = &record.predicted_name[..suffix_match.start()];

                Ok((virus_name, (record.uniprot_name, record.taxid)))
            } else {
                anyhow::bail!(
                    "Could not parse virus name from blast output: {:?}",
                    record.predicted_name
                );
            }
        })
        .context("Reading BLAST+ blastp output")?;

        Ok(VirusProteinTaxidMapping(mapping))
    }
}
