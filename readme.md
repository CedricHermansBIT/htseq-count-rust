# TallySeq

TallySeq is a fast Rust implementation of `htseq-count` semantics, designed as a practical drop-in counting replacement while retaining HTSeq-compatible results and file formats.

The `main` branch is tested against HTSeq 2.1.2 and supports the main `htseq-count` command-line surface, including multiple input files, single-end and paired-end data, SAM/BAM, stdin, compressed annotations, metadata columns, matrix output formats and annotated SAM/BAM output. CRAM input is supported on Linux and macOS.

## Installation

Build from source with a recent stable Rust toolchain:

```bash
git clone https://github.com/CedricHermansBIT/TallySeq
cd TallySeq
cargo build --release --locked
```

The primary binary is created at `target/release/tallyseq`. Release archives also include an `htseq-count` compatibility copy of the same executable.

## Usage

```bash
# Single alignment file
tallyseq -s no sample.bam genes.gtf

# Multiple samples in one count table
tallyseq -n 4 -s no sample1.bam sample2.bam sample3.cram genes.gtf.gz

# Read a single alignment stream from stdin
samtools view -h sample.bam | tallyseq -s no - genes.gtf

# Paired-end, coordinate-sorted input
tallyseq -r pos -s reverse paired.sorted.bam genes.gtf

# Include annotation metadata and a header
tallyseq \
  --additional-attr gene_name \
  --add-chromosome-info \
  --with-header \
  sample1.bam sample2.bam genes.gtf
```

The final positional file is the GTF/GFF annotation. All preceding positional files are alignment inputs.

### HTSeq-compatible options

- `-m, --mode`: `union`, `intersection-strict`, or `intersection-nonempty`. Default: `union`.
- `-s, --stranded`: `yes`, `no`, or `reverse`. Default: `yes`.
- `-t, --type`: feature type. Default: `exon`. May be supplied repeatedly.
- `-i, --idattr`: feature ID attribute. Default: `gene_id`. Repeated values are joined with `:`.
- `--additional-attr`: add an annotation attribute column to the output. May be repeated.
- `--add-chromosome-info`: add chromosome as an output metadata column.
- `--feature-query`: filter annotation features, for example `gene_name == "ACTB"`.
- `-a, --minaqual`: minimum MAPQ. Default: `10`.
- `--nonunique`: `none`, `all`, `fraction`, or `random`.
- `--secondary-alignments`: `ignore` or `score`.
- `--supplementary-alignments`: `ignore` or `score`.
- `-r, --order`: paired-end input order, `name` or `pos`.
- `--max-reads-in-buffer`: maximum unmatched-mate buffer for `--order pos`.
- `-n, --nprocesses`: number of alignment files processed concurrently.
- `-f, --format`: deprecated compatibility option. Accepted but ignored, matching modern HTSeq; input type is auto-detected.
- `-q, --quiet`: suppress progress output.
- `-d, --delimiter`: tabular output delimiter.
- `-c, --counts_output`: write counts to a file.
- `--with-header`: add input filenames as column headers.
- `--append-output`: append tabular output to an existing file.
- `--counts-output-sparse`: use sparse storage for supported matrix output formats.
- `-o, --samout`: write annotated alignments with the `XF` assignment tag. Supply once per alignment input.
- `-p, --samout-format`: `SAM` or `BAM`.
- `--version`: print the program version.

A Rust-specific `--threads` option controls BAM/CRAM decoding threads per input file. This is deliberately separate from HTSeq's `-n/--nprocesses`.

### Input compatibility

Normal local SAM and BAM files use the fast pure-Rust reader directly. On Linux and macOS, HTSlib is used as a compatibility bridge when needed for:

- CRAM
- stdin
- alignment files without a recognized extension
- format autodetection outside the direct SAM/BAM path

The Windows build keeps the dependency stack fully native and uses the pure-Rust SAM/BAM reader for files, stdin and format sniffing. CRAM input is not currently available in the Windows package.

GTF/GFF annotations may be plain text or gzip-compressed (`.gz` / `.gzip`). The annotation parser accepts the common GTF/GFF2 and GFF3 attribute separators and handles quoted semicolons correctly.

### Count output formats

Without `-c`, counts are written as a tabular table to stdout. With `-c`, the suffix selects the output format:

- `.tsv`, `.csv`, `.txt`: tabular output
- `.mtx`: Matrix Market plus `_features.tsv` and `_samples.tsv`
- `.h5ad`: AnnData/H5AD
- `.loom`: Loom

Matrix outputs use float32 values, matching HTSeq. `--counts-output-sparse` writes coordinate Matrix Market output and CSR-backed H5AD data.

### Annotated SAM/BAM output

`--samout` works for single-end and paired-end data, including name-grouped and position-sorted pairs. One output filename is required for each input alignment file.

```bash
tallyseq \
  -r pos \
  -o sample1.annotated.bam \
  -o sample2.annotated.bam \
  -p BAM \
  sample1.bam sample2.bam genes.gtf
```

The normal counting path does not create the SAM/BAM annotation machinery unless `--samout` is requested.

## Paired-end implementation

Two streaming strategies are used rather than loading every read name into memory:

- With `--order name`, only the current query-name group is retained and reciprocal mates are paired within that group.
- With `--order pos`, unmatched records are retained in a bounded hash map until their reciprocal mate is encountered.

Missing mates are counted using the alignment that is available, matching HTSeq. The second mate's strand interpretation is inverted in stranded counting, also matching HTSeq.

## HTSeq compatibility testing

The branch contains three complementary test layers:

1. generated differential tests against HTSeq 2.1.2;
2. adversarial interval-boundary tests and official upstream HTSeq fixtures;
3. an end-to-end compatibility suite for the broader CLI and file formats.

The compatibility suite covers:

- multiple alignment inputs and `-n/--nprocesses`
- headers, append mode and annotation metadata columns
- `--feature-query`
- gzipped annotations
- stdin and extension-independent alignment detection
- CRAM
- single-end and paired-end `--samout`
- SAM and BAM annotation output
- coordinate-sorted paired `samout`
- Matrix Market, sparse H5AD and Loom

Core counting differential coverage includes CIGAR `M`, `=`, `X`, insertions, deletions, skipped regions and clipping; all overlap and strandedness modes; MAPQ/NH behavior; secondary/supplementary alignments; repeated feature types and ID attributes; and both paired-end ordering modes.

One unusual behavior is deliberately retained for HTSeq 2.1.2 parity: if paired mate 1 exists but lacks an `NH` tag, current HTSeq does not inspect an `NH` tag present only on mate 2. The Rust implementation mirrors that behavior.

## Dependencies and reproducible builds

`Cargo.lock` is committed and CI builds with `cargo build --release --locked`. Release packaging is validated natively on Linux x86_64/ARM64, macOS x86_64/ARM64, and Windows x86_64.

The implementation uses the lightweight pure-Rust `bam` crate on the normal SAM/BAM hot path. On Unix platforms, HTSlib is included for CRAM and compatibility conversion; the Windows package omits HTSlib so it can build natively. `rust-hdf5` is used for native H5AD/Loom output.

### Optional feature-tree export

The Graphviz feature-tree export is disabled by default and is not part of normal counting.

```bash
tallyseq --export-feature-tree feature_tree.dot reads.bam genes.gtf
```

The older `--export_feature_map` spelling remains as an alias. There is intentionally no `-f` short form because `-f` belongs to HTSeq's deprecated `--format` option.

## Real-data benchmark

A reproducible benchmark is available in `benchmarks/benchmark_real_data.py`. It downloads real paired-end Pasilla RNA-seq chromosome 4 BAM files and the matching Drosophila BDGP5.78 GTF from Zenodo record 61771 and verifies their published MD5 checksums.

Every benchmark scenario first requires exact count parity with HTSeq. A differing feature or special count invalidates that benchmark run.

The `full` profile covers single-end and paired-end data, all three overlap modes, all three strandedness modes, multimapper options, MAPQ thresholds, secondary/supplementary scoring, repeated feature types/IDs, and both paired input orderings.

```bash
python -m pip install HTSeq
python benchmarks/benchmark_real_data.py \
  --build \
  --profile full \
  --rust-threads 1 \
  --nprocesses 1 \
  --repeats 3
```

Downloads are cached under `benchmarks/data`; JSON and Markdown results are written under `benchmarks/results`.

The benchmark records the TallySeq thread/process settings, CPU model, logical CPU count, total RAM, platform/kernel, Rust/Cargo versions, Git commit, and filesystem provenance.

A second benchmark, `benchmarks/benchmark_scaling.py`, measures scaling with alignment count. It deterministically cycles the real Pasilla-derived single-end records to 100,000, 1,000,000, 5,000,000, and 10,000,000 alignments by default. Exact HTSeq/TallySeq count equality is required at every target and repeat.

```bash
python benchmarks/benchmark_scaling.py \
  --build \
  --rust-threads 1 \
  --nprocesses 1 \
  --repeats 3
```

The scaling report includes median wall time, peak RSS, the HTSeq/TallySeq runtime ratio at every target, and a linear fit of runtime against millions of alignment records.

## Performance

An older benchmark used SRR5724993 aligned to GRCh38 with a GTF containing 1,065,949 genes and 58,663,336 alignment records:

| Tool | Time | Peak memory |
| --- | ---: | ---: |
| htseq-count | ~40 min | 1749 MB virtual, 150 MB resident |
| htseq-count-rust with `-o` | 3 min 15 s | 1947 MB virtual, 1351 MB resident |
| htseq-count-rust without `-o` | 1 min 54 s | 467 MB virtual, 135 MB resident |

These figures predate the current parity and paired-end work and should be rerun before using them as a current benchmark.
