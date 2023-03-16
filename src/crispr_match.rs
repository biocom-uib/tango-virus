use std::{io, fs::File};

use clap::Args;

use crate::util::{cli_tools, interned_mapping::InternedMultiMapping};

#[derive(Default)]
pub struct VirusHostMapping(pub InternedMultiMapping);

impl VirusHostMapping {
    pub fn read_tsv<R: io::Read>(&self, reader: R) -> anyhow::Result<Self> {
        Ok(Self(InternedMultiMapping::read_tsv(reader, true)?))
    }

    pub fn write_tsv<W: io::Write>(&self, writer: W) -> io::Result<()> {
        self.0.write_tsv(writer, &["virus_name", "host_name"])
    }
}

mod minced_spacers;

#[derive(Args)]
pub struct CrisprMatchArgs {
    /// Directory to place intermediate files in, instead of a temporary directory. Existing files
    /// will not be re-created.
    #[clap(long, short = 'd')]
    work_dir: Option<String>,

    /// Path to the MinCED .jar file.
    #[clap(long = "minced", value_name = "MINCED_JAR")]
    minced_jar: String,

    /// Installation prefix of the NCBI BLAST+ toolkit to use.
    #[clap(long)]
    blast_prefix: Option<String>,

    /// Number of threads (blastn option -num_threads).
    #[clap(long)]
    num_threads: Option<i32>,

    /// Minimum hit identity % for the blastn search (0-100).
    #[clap(long)]
    perc_identity: Option<i32>,

    #[clap(long)]
    dry_run: bool,

    /// FASTA file of the metaviromic sample.
    viral_seqs: String,

    /// FASTA file of the metagenomic sample.
    metagenomic_seqs: String,

    /// Output file (tab-separated values)
    #[clap(long, short)]
    output: String,
}

pub fn crispr_match(args: CrisprMatchArgs) -> anyhow::Result<()> {
    let Some(java) = cli_tools::Java::resolve()? else {
        anyhow::bail!("Could not locate a Java executable in PATH");
    };

    let Some(minced) = minced_spacers::MincedTool::resolve(java, args.minced_jar.as_ref()) else {
        anyhow::bail!("Could not locate the provided MinCED .jar file");
    };

    let blast_prefix = args.blast_prefix.as_ref().map(|s| s.as_ref());

    let Some(blastn) = cli_tools::BlastTool::resolve(blast_prefix, "blastn")? else {
        anyhow::bail!("Could not find a blastn executable");
    };

    let viral_seqs = args.viral_seqs.as_ref();
    let metagenomic_seqs = args.metagenomic_seqs.as_ref();

    let run_pipeline = |work_dir| {
        let mut pipeline = minced_spacers::MincedSpacersPipeline::new(minced, blastn, work_dir);
        pipeline.set_blastn_num_threads(args.num_threads);
        pipeline.set_dry_run(args.dry_run);
        pipeline.match_spacers(viral_seqs, metagenomic_seqs, args.perc_identity)?;
        pipeline.collect_virus_host_mapping()
    };

    let mapping = if let Some(work_dir) = args.work_dir {
        let work_dir = work_dir.as_ref();
        run_pipeline(work_dir)?
    } else {
        let temp = tempdir::TempDir::new("crispr_match")?;
        run_pipeline(temp.path())?
    };

    let output_file = File::create(&args.output)?;
    mapping.write_tsv(output_file)?;

    Ok(())
}
