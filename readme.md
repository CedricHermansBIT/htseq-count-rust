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

One unusual behavior is deliberately retained for HTSeq 2.1.2 parity: if paired mate 1 exists but lacks an `NH` tag, HTSeq 2.1.2 does not inspect an `NH` tag present only on mate 2. TallySeq mirrors that behavior.

## Dependencies and reproducible builds

`Cargo.lock` is committed and CI builds with `cargo build --release --locked`. Release packaging is validated natively on Linux x86_64/ARM64, macOS x86_64/ARM64, and Windows x86_64.

The implementation uses the lightweight pure-Rust `bam` crate on the normal SAM/BAM hot path. On Unix platforms, HTSlib is included for CRAM and compatibility conversion; the Windows package omits HTSlib so it can build natively. `rust-hdf5` is used for native H5AD/Loom output.

### Optional feature-tree export

The Graphviz feature-tree export is disabled by default and is not part of normal counting.

```bash
tallyseq --export-feature-tree feature_tree.dot reads.bam genes.gtf
```

`--export_feature_map` is retained as an alias. There is intentionally no `-f` short form because `-f` belongs to HTSeq's deprecated `--format` option.

## Real-data benchmark

A reproducible benchmark is available in `benchmarks/benchmark_real_data.py`. It downloads real paired-end Pasilla RNA-seq chromosome 4 BAM files and the matching Drosophila BDGP5.78 GTF from Zenodo record 61771 and verifies their published MD5 checksums.

Every benchmark scenario first requires exact count parity with HTSeq. A differing feature or special count invalidates that benchmark run.

The `full` profile covers single-end and paired-end data, all three overlap modes, all three strandedness modes, multimapper options, MAPQ thresholds, secondary/supplementary scoring, repeated feature types/IDs, and both paired input orderings.

```bash
python -m pip install "HTSeq==2.1.2"
python benchmarks/benchmark_real_data.py \
  --build \
  --profile full \
  --dataset GSM461177 \
  --dataset GSM461178 \
  --rust-threads 1 \
  --nprocesses 1 \
  --repeats 5
```

Downloads are cached under `benchmarks/data`; JSON and Markdown results are written under `benchmarks/results`.

The benchmark records the TallySeq thread/process settings, CPU model, logical CPU count, total RAM, platform/kernel, Rust/Cargo versions, Git commit, and filesystem provenance.

A second benchmark, `benchmarks/benchmark_scaling.py`, measures scaling with alignment count. It deterministically cycles the real Pasilla-derived single-end records to 100,000, 1,000,000, 5,000,000, and 10,000,000 alignments by default. Exact HTSeq/TallySeq count equality is required at every target and repeat.

```bash
python benchmarks/benchmark_scaling.py \
  --build \
  --dataset GSM461177 \
  --rust-threads 1 \
  --nprocesses 1 \
  --repeats 5 \
  --records 100000 \
  --records 1000000 \
  --records 5000000 \
  --records 10000000
```

The scaling report includes median wall time, peak RSS, the HTSeq/TallySeq runtime ratio at every target, and a linear fit of runtime against millions of alignment records.

## Performance

The release benchmark uses HTSeq 2.1.2 and TallySeq 0.1.0 with TallySeq `--threads 1` and both tools `-n 1`. Each real-data dataset/scenario combination is measured five times after exact normalized count equality is confirmed. The benchmark covers 31 scenarios on each of two Pasilla chromosome 4 datasets, for 62 dataset/scenario comparisons and 620 timed executions.

Across those comparisons, the median of the per-comparison median runtimes is 0.58 s for TallySeq and 14.59 s for HTSeq. The median pairwise runtime ratio is 26.1x, with a range of 18.7x to 29.1x. Median peak RSS is 28.1 MiB for TallySeq and 73.8 MiB for HTSeq.

The controlled read-count scaling benchmark keeps the real Pasilla alignment distribution fixed while cycling records to exact target sizes:

| Alignments | TallySeq median | HTSeq median | HTSeq/TallySeq | TallySeq RSS | HTSeq RSS |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100,000 | 0.53 s | 14.65 s | 27.6x | 25.8 MiB | 67.6 MiB |
| 1,000,000 | 1.29 s | 28.25 s | 21.9x | 25.8 MiB | 67.6 MiB |
| 5,000,000 | 5.27 s | 100.50 s | 19.1x | 25.9 MiB | 67.7 MiB |
| 10,000,000 | 10.02 s | 221.74 s | 22.1x | 25.8 MiB | 67.5 MiB |

Linear fits give 0.964 s per million alignments for TallySeq (R² = 0.9998) and 20.911 s per million for HTSeq (R² = 0.9931). All reported scaling runs passed exact normalized count equality.

The benchmark environment and toolchain are recorded in the generated JSON and manuscript provenance files, including CPU, logical CPU count, RAM, kernel/platform, Python, Rust/Cargo, Git commit, filesystem, thread/process settings, and execution order.
