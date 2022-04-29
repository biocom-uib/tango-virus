# TANGO-Virus

## Build instructions

First, you will need `cargo`. To install it, we recommend using `rustup`, which
can be installed easily as follows (this includes `cargo`):

```
curl https://sh.rustup.rs -sSf | sh -s -- --profile minimal --default-toolchain nightly
```

Then you can install `tango-virus` as follows (by default, it is installed in
`~/.cargo/bin`):

```
cargo install --git https://github.com/biocom-uib/tango-virus
```

## Usage

TANGO-Virus has four main subcommands, which are usually run the same order as
explained. To obtain more information on each command, use `tango-virus
<subcommand> --help`.

### `preprocess-taxonomy`

This subcommand loads, preprocesses and serializes the provided taxonomy (i.e.
NCBI's). To use the NCBI taxonomy, download and extract it (it can be found
[here](https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz)). Assuming that
the extracted taxonomy is in `/path/to/taxdump/`, you could run

```
tango-virus preprocess-taxonomy /path/to/taxdump/ --contract preprocessed_taxonomy.cbor
```

_Note that this step only needs to be run once per taxonomy update._

### `preprocess-blastout`

This subcommand loads the output of NCBI's BLAST+ (`blastn`) to produce a
suitable input file for the `assign` subcommand. Assuming that `blastn` was
executed using `-outfmt '6 qseqid staxid'` or `-outfmt '7 qseqid staxid'`
and, for instance, `-db nt`:

```
tango-virus preprocess-blastout /path/to/blast/hits.out \
    --blast-outfmt '7 qseqid staxid' \
    --query-id-col qseqid \
    --subject-id-col staxid \
    -o preprocessed_hits.tsv
```

The argument to `--blast-outfmt` should be exactly the same as `blastn`'s
`-outfmt`, but only supports formats 6 and 7. Columns `staxid` and `ssciname`
are not included by default but at least one of those is required by `assign`
(`staxid` is preferred), see `blastn -help` for more information about
specifying output columns).

### `assign`

This subcommand assigns hits from `preprocess-blastout` to nodes in the
taxonomy. This is a reimplementation of
[TANGO](https://www.cs.upc.edu/~valiente/tango/) using the serialized taxonomy
from the `preprocess-taxonomy` step and the output from `preprocess-blastout`.
To run it, execute

```
tango-virus assign preprocessed_taxonomy.cbor preprocessed_hits.tsv \
    -q 0.5 -o assignment.tsv
```

### `refine-vpf-class`

This step takes as input the results from `assign` and the host prediction
output from `vpf-class` ([GitHub](https://github.com/biocom-uib/vpf-tools)) and
refines it to include only ascendants from taxonomic nodes found in the
`assign` step. If the output is located at `/path/to/vpf-class/output/`, then
the refinement of the `host_genus` prediction can be done as follows:

```
tango-virus refine-vpf-class preprocessed_taxonomy.cbor assignment.tsv \
    /path/to/vpf-class/output/host_genus.tsv \
    -o refined_host_genus.tsv
```

(*TODO*) The output contains two additional columns: for each VPF-Class
prediction, `assignments` is a semicolon-separated list of `taxid`s found in
`assignment.tsv` which are descendants of this prediction and
`assigned_contigs` is a semicolon-separated list of the metagenomic contigs
that were grouped as those descendants.
