# htseq-count-rust

This is a Rust implementation of the `htseq-count` tool from the HTSeq package. It is a tool for counting the number of reads mapping to each feature in a GFF file. It is designed to be used with the output of a read aligner such as STAR or HISAT2.

This implementation is designed to be faster than the original Python implementation, but it is not yet feature complete. It currently only supports single-end reads and only one idattr.

## Installation

To install the tool, you will need to have Rust installed. The best thing is to build the tool from source:

```bash
git clone https://github.com/CedricHermansBIT/htseq-count-rust
cd htseq-count-rust
cargo build --release
```

This will create a binary in `target/release/htseq-count-rust`. You can copy this binary to a location in your PATH, or you can run it directly from the `target/release` directory.

## Usage

The tool is used in the same way as the original `htseq-count` tool. You need to provide a GTF file and a SAM/BAM file. The SAM/BAM file can be piped in from a read aligner such as STAR or HISAT2. For example:

```bash
# Simple usage
htseq-count-rust alignments.bam genes.gtf

# With options
htseq-count-rust -s no -t exon -i gene_id -m intersection-strict -d , -a 10 -n 8 --nonunique all --secondary-alignments ignore --supplementary-alignments ignore -o alignments_out.sam alignments.bam genes.gtf
```

Note: currently, only one SAM/BAM file can be provided. If you want to run the tool on multiple files, you will need to run them separately.

### Options

The following options are currently supported:

- `-m` | `--mode`:     Mode. One of `union`, `intersection-strict`, or `intersection-nonempty`. Default is `union`. Note that intersection-nonempty does not fully give the same results as the original `htseq-count` tool.
- `-s` | `--stranded`: Strandedness. One of `yes`, `no`, or `reverse`. Default is `yes`.
- `-t` | `--type`:     Feature type. Default is `exon`.
- `-i` | `--idattr`:   Attribute to use as feature ID. Default is `gene_id`.
- `-a` | `--minaqual`: Read mapping quality. Reads with mapping quality less than this value will be ignored. Default is 10.
- `-n` | `--threads`:  Number of threads to use. Default is 4. Note that this is mainly for decompression of BAM files, and the actual counting is single-threaded. An additional thread is always used with the -o option for writing output to a SAM file.
- `-o` | `--samout`:   Write out all SAM alignment records into a SAM file. Each record will get an additional `XF` tag with the feature ID and the feature type.
- `-d` | `--delimiter`: Delimiter for feature IDs. Default is `\t`.
- `--nonunique`: One of `none`, `all`, `fraction` or `random`. How to handle non-unique features. Default is `none`.
- `--secondary-alignments`: One of `score` or `ignore`. Treat secondary alignments as distinct records. Default is to score them.
- `--supplementary-alignments`: One of `score` or `ignore`. Treat supplementary alignments as distinct records. Default is to score them.

## Performance

In the future, better performance comparisson will be made. For now, the following is a simple comparisson of the original `htseq-count` and this implementation.

### Test data
SRR5724993
The test data consists of a GTF file with 1,065,949 genes and a SAM file with 58,663,336 reads. The BAM file is derived from SRR5724993 and was aligned with HISAT2 to the human genome (GRCh38).

Note that memory usage highly depends on the write speed of the disk. Since the file is processed so fast, the writer can not keep up and the data to write is buffered in memory. This is why the memory usage is sometimes higher than the peak memory usage of the original `htseq-count` tool.

Tool | Time | Memory (peak)
--- | --- | ---
htseq-count | +-40m | 1749Mb Virtual, 150Mb Resident
htseq-count-rust | 3m15s | 1947Mb Virtual, 1351Mb Resident 
htseq-count-rust (without -o) | 1m54s | 467M Virtual, 135M Resident