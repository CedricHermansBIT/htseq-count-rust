# htseq-count-rust

A Rust implementation of `htseq-count` for counting aligned reads or read pairs against GTF/GFF features. The goal is count compatibility with HTSeq while keeping BAM/SAM processing fast and memory-efficient.

The implementation currently supports single-end and paired-end SAM/BAM input, all three HTSeq overlap modes, stranded counting, repeated feature types and ID attributes, multimapper handling, and secondary/supplementary alignment controls.

## Installation

Build from source with a recent stable Rust toolchain:

```bash
git clone https://github.com/CedricHermansBIT/htseq-count-rust
cd htseq-count-rust
cargo build --release --locked
```

The binary is created at `target/release/htseq_count_rust`.

## Usage

```bash
# Single-end
htseq_count_rust -s no alignments.bam genes.gtf

# Paired-end, query-name grouped input
htseq_count_rust -r name -s reverse alignments.name.bam genes.gtf

# Paired-end, coordinate-sorted input
htseq_count_rust -r pos -s reverse alignments.sorted.bam genes.gtf
```

Only one alignment file is processed per invocation.

### Main options

- `-m, --mode`: `union`, `intersection-strict`, or `intersection-nonempty`. Default: `union`.
- `-s, --stranded`: `yes`, `no`, or `reverse`. Default: `yes`.
- `-t, --type`: feature type to count. Default: `exon`. May be supplied multiple times.
- `-i, --idattr`: GTF/GFF attribute used as the feature ID. Default: `gene_id`. May be supplied multiple times; values are joined with `:`, matching HTSeq.
- `-a, --minaqual`: minimum MAPQ. Default: `10`.
- `--nonunique`: `none`, `all`, `fraction`, or `random`. Default: `none`.
- `--secondary-alignments`: `ignore` or `score`. Default: `ignore`.
- `--supplementary-alignments`: `ignore` or `score`. Default: `ignore`.
- `-r, --order`: paired-end input order, `name` or `pos`. Default: `name`. Ignored for single-end input.
- `--max-reads-in-buffer`: maximum number of unmatched mate keys retained for coordinate-sorted paired-end data. Default: `30000000`.
- `-n, --threads`: BAM decompression threads. Counting itself is currently single-threaded.
- `-c, --counts_output`: write counts to a file instead of stdout.
- `-o, --samout`: annotate single-end SAM output with `XF`. Paired-end `--samout` is not implemented yet.
- `-d, --delimiter`: output delimiter. Default: tab.
- `--extended-output`: also print the calculated number of uniquely mapped reads/fragments.

## Paired-end implementation

Two streaming strategies are used rather than loading all read names into memory:

- With `--order name`, only the current query-name group is retained. Alignments are paired using reciprocal mate status and coordinates, following HTSeq's name-sorted pairing behavior.
- With `--order pos`, unmatched records are kept in a hash map keyed by query name, read number, alignment position, mate position and template length. A record is removed as soon as its reciprocal mate is encountered. This follows the same general bounded-buffer strategy used by HTSeq for position-sorted input.

Secondary and supplementary records are removed before entering the pairing buffer when both are configured as `ignore`. Missing mates are still counted using the mate that is available, matching HTSeq.

## HTSeq compatibility testing

The `diagnostics` branch contains differential tests that run this binary and HTSeq 2.1.2 on the same generated input and compare every feature and special count.

Coverage includes:

- CIGAR `M`, `=`, `X`, insertions, deletions, skips and soft clipping
- `union`, `intersection-strict`, and `intersection-nonempty`
- stranded, unstranded and reverse-stranded counting
- repeated `-t` and `-i`
- MAPQ and `NH` handling
- secondary and supplementary alignments
- name-grouped and coordinate-sorted paired-end input
- missing and unmapped mates
- deterministic randomized differential cases
- official HTSeq test fixtures, including its position-sorted paired BAM

One unusual behavior is deliberately retained for compatibility with HTSeq 2.1.2: for paired reads, if mate 1 exists but has no `NH` tag, HTSeq's current implementation does not inspect an `NH` tag that is present only on mate 2. The Rust implementation mirrors that behavior so counts remain comparable.

Randomized cases where HTSeq itself raises an exception are reported as reference errors and are not treated as evidence of a Rust mismatch.

## Dependency status

The direct dependencies are kept at their current stable releases:

- `clap 4.6.7` with derive support
- `bam 0.1.4` (this is still the newest release of the pure-Rust `bam` crate)
- `rand 0.10.3`

`Cargo.lock` is committed and CI builds with `--locked`.

## Real-data benchmark

A reproducible benchmark is available in `benchmarks/benchmark_real_data.py`. It downloads real paired-end Pasilla RNA-seq chromosome 4 BAM files and the matching Drosophila BDGP5.78 GTF from Zenodo record 61771, verifies the published MD5 checksums, then benchmarks this implementation against HTSeq with identical counting options.

The benchmark treats count parity as a hard requirement. Feature IDs, special rows and numeric counts must be exactly equal. No tolerance is used. If even one count differs, the benchmark stops and reports the first differences instead of producing a successful performance report.

Run locally with:

```bash
python -m pip install HTSeq
python benchmarks/benchmark_real_data.py --build --repeats 3
```

To benchmark both included real BAMs:

```bash
python benchmarks/benchmark_real_data.py --build \\
  --dataset GSM461177 \\
  --dataset GSM461178 \\
  --repeats 3
```

Downloads are cached under `benchmarks/data`. Results are written as JSON and Markdown under `benchmarks/results`, including wall time, peak resident memory and the exact-parity result for every measured repeat. The script alternates which tool runs first between repeats to reduce systematic page-cache bias.

There is also a manual GitHub Actions workflow named `Real-data benchmark`, so the same benchmark can be launched from the Actions tab without preparing the dataset locally.

## Performance

An older benchmark used SRR5724993 aligned to GRCh38 with a GTF containing 1,065,949 genes and 58,663,336 alignment records:

| Tool | Time | Peak memory |
| --- | ---: | ---: |
| htseq-count | ~40 min | 1749 MB virtual, 150 MB resident |
| htseq-count-rust with `-o` | 3 min 15 s | 1947 MB virtual, 1351 MB resident |
| htseq-count-rust without `-o` | 1 min 54 s | 467 MB virtual, 135 MB resident |

These figures predate the current parity and paired-end work and should be rerun before using them as a current benchmark.
