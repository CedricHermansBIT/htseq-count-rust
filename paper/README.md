# TallySeq preprint

This directory contains the working preprint for TallySeq.

## Current status

The manuscript is an initial draft. The benchmark results currently embedded in the paper are preliminary CI measurements from the real-data benchmark on the Pasilla chromosome 4 data set. The full matrix covered 31 scenarios and required exact count equality with HTSeq 2.1.2.

Before submission, replace the preliminary CI benchmark with repeated measurements on dedicated hardware. The benchmark script already records wall time and peak RSS and rejects any scenario with unequal counts.

## Rebuilding benchmark values

After running:

```bash
python benchmarks/benchmark_real_data.py \
  --build \
  --profile full \
  --dataset GSM461177 \
  --dataset GSM461178 \
  --repeats 5
```

update the paper snapshot and LaTeX summary macros with:

```bash
python paper/update_benchmark.py benchmarks/results/real_data_benchmark.json
```

The script rewrites `benchmark_snapshot.tsv` and `benchmark_macros.tex`. Review the manuscript text and representative table after updating the benchmark.

## Compile

With a standard TeX installation:

```bash
cd paper
latexmk -pdf main.tex
```

The draft uses standard packages plus `siunitx`, `natbib`, `authblk`, `booktabs`, `microtype`, and `hyperref`.

## Preliminary benchmark provenance

The checked-in snapshot was extracted from GitHub Actions workflow run `36193504716` at commit `ad54b53a`. That run used one measured repetition for each of 31 scenarios. A later three-repeat paired-end union smoke run at commit `ac43dbdc` measured TallySeq at 0.46 s in all three repetitions and HTSeq at 18.77 s, 19.13 s, and 19.00 s, with all 15,687 count rows identical in every repeat.
