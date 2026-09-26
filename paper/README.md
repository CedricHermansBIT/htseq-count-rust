# TallySeq preprint

This directory contains the TallySeq manuscript and the benchmark snapshots used in it.

## Benchmark data in the manuscript

The manuscript uses the repeated real-data benchmark and controlled read-count scaling benchmark.

The real-data benchmark covers 31 scenarios on each of two Pasilla chromosome 4 datasets, GSM461177 and GSM461178. Each dataset/scenario combination is measured five times for TallySeq and five times for HTSeq 2.1.2 after exact normalized count equality is confirmed. TallySeq uses `--threads 1`, and both tools use `-n 1`.

The scaling benchmark uses real single-end Pasilla alignments and deterministically cycles them to 100,000, 1,000,000, 5,000,000, and 10,000,000 records. It does not alter alignment coordinates, CIGAR strings, MAPQ values, tags, strands, or sequences. Each target and repeat must pass exact count equality before its timing is accepted.

The checked-in manuscript support files are:

- `benchmark_snapshot.tsv`: all 62 real-data dataset/scenario comparisons
- `benchmark_macros.tex`: aggregate real-data values used by the manuscript
- `benchmark_provenance.json`: real-data benchmark settings and system provenance
- `scaling_snapshot.tsv`: the four controlled read-count scaling targets
- `scaling_macros.tex`: scaling values used by the manuscript
- `scaling_provenance.json`: scaling settings, system provenance, and fitted coefficients

## Rebuilding benchmark values

Run the real-data benchmark:

```bash
python benchmarks/benchmark_real_data.py \
  --build \
  --profile full \
  --dataset GSM461177 \
  --dataset GSM461178 \
  --rust-threads 1 \
  --nprocesses 1 \
  --repeats 5
```

Run the controlled scaling benchmark:

```bash
python benchmarks/benchmark_scaling.py \
  --rust-bin target/release/tallyseq \
  --dataset GSM461177 \
  --rust-threads 1 \
  --nprocesses 1 \
  --repeats 5 \
  --records 100000 \
  --records 1000000 \
  --records 5000000 \
  --records 10000000
```

Regenerate all manuscript benchmark support files with:

```bash
python paper/update_benchmark.py \
  benchmarks/results/real_data_benchmark.json \
  --scaling-report benchmarks/results/scaling_benchmark.json
```

## Compile

With a standard TeX installation:

```bash
cd paper
latexmk -pdf main.tex
```

The manuscript uses standard packages plus `siunitx`, `natbib`, `authblk`, `booktabs`, `microtype`, and `hyperref`.
