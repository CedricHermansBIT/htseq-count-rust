#!/usr/bin/env python3
"""Large-BAM scaling benchmark for TallySeq vs HTSeq.

This benchmark starts from the real Pasilla chromosome-4 data used by
benchmark_real_data.py, derives the same real single-end mate-1 alignments, and
then deterministically cycles those records to exact target record counts.

The generated BAMs therefore preserve real coordinates, CIGAR strings, MAPQ
values, tags, strand assignments, and sequence content while allowing controlled
scaling of alignment count. Exact HTSeq/TallySeq count equality is required for
every target size and repeat before performance results are reported.
"""

from __future__ import annotations

import argparse
import json
import statistics
from dataclasses import asdict
from pathlib import Path

from benchmark_real_data import (
    DOI,
    FILES,
    Measurement,
    counts,
    derive_single_end,
    download_verified,
    exact_differences,
    executable,
    get_rust_binary,
    human_bytes,
    md5sum,
    require_pysam,
    summarize,
    system_provenance,
    timed,
    version,
)

DEFAULT_RECORD_COUNTS = (100_000, 1_000_000, 5_000_000, 10_000_000)
COUNT_OPTIONS = [
    "-s",
    "no",
    "-m",
    "union",
    "-a",
    "10",
    "--nonunique",
    "none",
    "--secondary-alignments",
    "ignore",
    "--supplementary-alignments",
    "ignore",
    "-t",
    "exon",
    "-i",
    "gene_id",
]


def derive_scaled_single_end(
    source: Path,
    destination: Path,
    target_records: int,
) -> Path:
    """Cycle real single-end records until exactly target_records are written."""
    pysam = require_pysam()
    source_md5 = md5sum(source)
    metadata_path = destination.with_suffix(destination.suffix + ".scale.json")
    expected_metadata = {
        "generator_version": 1,
        "source_md5": source_md5,
        "target_records": target_records,
    }

    if destination.exists() and metadata_path.exists():
        try:
            metadata = json.loads(metadata_path.read_text())
        except (OSError, json.JSONDecodeError):
            metadata = None
        if metadata == expected_metadata:
            print(
                f"[cache] {destination.name} "
                f"({target_records:,} records, {human_bytes(destination.stat().st_size)})"
            )
            return destination

    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_name(destination.name + ".part")
    tmp.unlink(missing_ok=True)

    written = 0
    with pysam.AlignmentFile(str(source), "rb") as template:
        header = template.header.to_dict()

    with pysam.AlignmentFile(str(tmp), "wb", header=header) as writer:
        while written < target_records:
            before = written
            with pysam.AlignmentFile(str(source), "rb") as reader:
                for record in reader.fetch(until_eof=True):
                    writer.write(record)
                    written += 1
                    if written >= target_records:
                        break
            if written == before:
                raise RuntimeError(f"Source BAM {source} contains no records")

    tmp.replace(destination)
    metadata_path.write_text(json.dumps(expected_metadata, indent=2) + "\n")
    print(
        f"[derived] {destination.name}: {written:,} records "
        f"({human_bytes(destination.stat().st_size)})"
    )
    return destination


def linear_fit(points: list[tuple[int, float]]) -> dict:
    """OLS fit of wall seconds against millions of alignment records."""
    xs = [records / 1_000_000.0 for records, _ in points]
    ys = [seconds for _, seconds in points]
    x_mean = statistics.mean(xs)
    y_mean = statistics.mean(ys)
    denominator = sum((x - x_mean) ** 2 for x in xs)

    if denominator == 0:
        slope = 0.0
    else:
        slope = sum(
            (x - x_mean) * (y - y_mean) for x, y in zip(xs, ys)
        ) / denominator
    intercept = y_mean - slope * x_mean

    ss_total = sum((y - y_mean) ** 2 for y in ys)
    ss_residual = sum(
        (y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys)
    )
    r_squared = 1.0 if ss_total == 0 else 1.0 - ss_residual / ss_total

    return {
        "intercept_seconds": intercept,
        "seconds_per_million_records": slope,
        "r_squared": r_squared,
    }


def write_markdown(report: dict, path: Path) -> None:
    settings = report["benchmark_settings"]
    system = report["system"]
    lines = [
        "# TallySeq large-BAM scaling benchmark",
        "",
        f"Source: Zenodo DOI {DOI}",
        f"Dataset: {report['dataset']}",
        f"HTSeq version: {report['htseq_version']}",
        f"TallySeq version: {report['rust_version']}",
        f"TallySeq threads: {settings['tallyseq_threads']}",
        f"Processes (-n): {settings['nprocesses']}",
        f"Repeats: {settings['repeats']}",
        f"CPU: {system.get('cpu_model') or 'unknown'}",
        f"Logical CPUs: {system.get('cpu_count_logical')}",
        f"Platform: {system.get('platform')}",
        f"Rust: {system.get('rustc') or 'unknown'}",
        "",
        "Every reported run passed exact normalized count equality. "
        "No numeric tolerance was used.",
        "",
        "| Records | BAM size (MiB) | Tool | Runs | Median time (s) | Best (s) | Median RSS (MiB) | HTSeq/TallySeq speed ratio |",
        "| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |",
    ]

    for target in report["record_counts"]:
        summary = report["summaries"][str(target)]
        for tool in ("rust", "htseq"):
            row = summary[tool]
            ratio = summary["htseq_over_rust_time_ratio"]
            lines.append(
                f"| {target:,} | "
                f"{report['scaled_inputs'][str(target)]['bytes'] / 1024 / 1024:.1f} | "
                f"{tool} | {row['runs']} | "
                f"{row['median_wall_seconds']:.3f} | "
                f"{row['best_wall_seconds']:.3f} | "
                f"{row['median_max_rss_kib'] / 1024:.1f} | "
                f"{ratio:.2f}x |"
            )

    lines += ["", "## Runtime scaling fit", ""]
    for tool, label in (("rust", "TallySeq"), ("htseq", "HTSeq")):
        fit = report["linear_scaling"][tool]
        lines.append(
            f"- {label}: intercept {fit['intercept_seconds']:.3f} s; "
            f"{fit['seconds_per_million_records']:.3f} s per million records; "
            f"R^2={fit['r_squared']:.4f}."
        )

    lines += [
        "",
        "## Method",
        "",
        "The benchmark input is generated by repeatedly cycling the real derived "
        "single-end Pasilla mate-1 alignments until each exact target record "
        "count is reached. Alignment coordinates, CIGAR strings, MAPQ values, "
        "tags, strands, and sequences are not synthesized or modified by the "
        "scaling step.",
        "",
        "TallySeq and HTSeq are run with equivalent count options using union, "
        "unstranded exon counting by gene_id at MAPQ >= 10. Execution order is "
        "alternated by repeat to reduce systematic page-cache bias.",
        "",
    ]
    path.write_text("\n".join(lines))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".")
    ap.add_argument(
        "--dataset",
        choices=("GSM461177", "GSM461178"),
        default="GSM461177",
    )
    ap.add_argument("--cache-dir", default="benchmarks/data")
    ap.add_argument("--results-dir", default="benchmarks/results")
    ap.add_argument("--rust-bin")
    ap.add_argument("--htseq-bin", default="htseq-count")
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument(
        "--records",
        action="append",
        type=int,
        help=(
            "Target BAM record count; may be repeated. Defaults to "
            "100k, 1M, 5M, and 10M."
        ),
    )
    ap.add_argument(
        "--rust-threads",
        type=int,
        default=1,
        help="TallySeq BAM decoding threads. Default: 1.",
    )
    ap.add_argument(
        "--nprocesses",
        type=int,
        default=1,
        help="Alignment files processed concurrently by both tools. Default: 1.",
    )
    args = ap.parse_args()

    if args.repeats < 1:
        ap.error("--repeats must be >= 1")
    if args.rust_threads < 1:
        ap.error("--rust-threads must be >= 1")
    if args.nprocesses < 1:
        ap.error("--nprocesses must be >= 1")

    record_counts = sorted(set(args.records or DEFAULT_RECORD_COUNTS))
    if not record_counts or record_counts[0] < 1:
        ap.error("--records values must be positive")

    repo = Path(args.repo).expanduser().resolve()
    cache = Path(args.cache_dir)
    results = Path(args.results_dir)
    if not cache.is_absolute():
        cache = repo / cache
    if not results.is_absolute():
        results = repo / results
    cache.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    run_dir = results / "scaling_runs"
    run_dir.mkdir(parents=True, exist_ok=True)

    rust = get_rust_binary(repo, args.rust_bin, args.build)
    htseq = executable(args.htseq_bin)
    gtf = download_verified("gtf", cache)
    paired = download_verified(args.dataset, cache)
    real_single = derive_single_end(
        paired, cache / f"{args.dataset}.mate1.single.bam"
    )

    print(f"TallySeq: {rust} ({version(str(rust), repo)})")
    print(f"HTSeq:    {htseq} ({version(htseq, repo)})")
    print(f"Dataset:  {args.dataset}")
    print(f"Targets:  {', '.join(f'{x:,}' for x in record_counts)}")
    print(f"Repeats:  {args.repeats}")
    print(f"TallySeq threads: {args.rust_threads}")
    print(f"Processes (-n):   {args.nprocesses}")

    scaled_inputs: dict[int, Path] = {}
    for target in record_counts:
        scaled_inputs[target] = derive_scaled_single_end(
            real_single,
            cache / f"{args.dataset}.single.scaled.{target}.bam",
            target,
        )

    measurements: list[Measurement] = []

    for target in record_counts:
        bam = scaled_inputs[target]
        scenario = f"scale_{target}"
        rust_cmd = [
            str(rust),
            "--threads",
            str(args.rust_threads),
            "-n",
            str(args.nprocesses),
            *COUNT_OPTIONS,
            str(bam),
            str(gtf),
        ]
        htseq_cmd = [
            htseq,
            "-n",
            str(args.nprocesses),
            *COUNT_OPTIONS,
            str(bam),
            str(gtf),
        ]

        print(
            f"\n[target] {target:,} records "
            f"({human_bytes(bam.stat().st_size)})"
        )

        for repeat in range(1, args.repeats + 1):
            order = (
                [("rust", rust_cmd), ("htseq", htseq_cmd)]
                if repeat % 2
                else [("htseq", htseq_cmd), ("rust", rust_cmd)]
            )
            outputs: dict[str, Path] = {}

            for tool, command in order:
                stem = f"{args.dataset}.{scenario}.r{repeat}.{tool}"
                out = run_dir / f"{stem}.counts.tsv"
                err = run_dir / f"{stem}.stderr.txt"
                timing = run_dir / f"{stem}.time.tsv"
                measurement = timed(
                    tool,
                    args.dataset,
                    scenario,
                    repeat,
                    command,
                    out,
                    err,
                    timing,
                    repo,
                )
                measurements.append(measurement)
                outputs[tool] = out
                print(
                    f"  repeat {repeat} {tool:5s}: "
                    f"{measurement.wall_seconds:.3f}s, "
                    f"{measurement.max_rss_kib / 1024:.1f} MiB RSS"
                )
                if measurement.exit_code:
                    raise RuntimeError(
                        f"{tool} exited {measurement.exit_code}; see {err}"
                    )

            diff = exact_differences(outputs["rust"], outputs["htseq"])
            if diff:
                details = "\n".join(
                    f"  {feature}: rust={rv}, htseq={hv}"
                    for feature, rv, hv in diff[:25]
                )
                raise RuntimeError(
                    f"{target:,} records / repeat {repeat}: "
                    f"{len(diff)} exact count differences\n{details}"
                )

            print(
                f"  [exact] {len(counts(outputs['rust']))} count rows identical"
            )

    summaries: dict[str, dict] = {}
    rust_fit_points: list[tuple[int, float]] = []
    htseq_fit_points: list[tuple[int, float]] = []

    for target in record_counts:
        scenario = f"scale_{target}"
        rs = summarize(measurements, "rust", args.dataset, scenario)
        hs = summarize(measurements, "htseq", args.dataset, scenario)
        summaries[str(target)] = {
            "rust": rs,
            "htseq": hs,
            "htseq_over_rust_time_ratio": (
                hs["median_wall_seconds"] / rs["median_wall_seconds"]
            ),
        }
        rust_fit_points.append((target, rs["median_wall_seconds"]))
        htseq_fit_points.append((target, hs["median_wall_seconds"]))

    report = {
        "schema_version": 1,
        "benchmark": "large_bam_scaling",
        "zenodo_doi": DOI,
        "dataset": args.dataset,
        "exact_count_equality": True,
        "comparison": (
            "Exact feature-key and Decimal numeric-value equality; "
            "no tolerance. Row order is normalized."
        ),
        "record_counts": record_counts,
        "benchmark_settings": {
            "repeats": args.repeats,
            "tallyseq_threads": args.rust_threads,
            "nprocesses": args.nprocesses,
            "execution_order": (
                "alternated by repeat to reduce systematic page-cache bias"
            ),
            "count_options": COUNT_OPTIONS,
        },
        "rust_binary": str(rust),
        "rust_version": version(str(rust), repo),
        "htseq_binary": htseq,
        "htseq_version": version(htseq, repo),
        "system": system_provenance(repo),
        "source": {
            "paired_bam": str(paired),
            "paired_bam_md5": FILES[args.dataset][1],
            "derived_single_bam": str(real_single),
            "derived_single_bam_md5": md5sum(real_single),
            "gtf": str(gtf),
            "gtf_md5": FILES["gtf"][1],
        },
        "scaled_inputs": {
            str(target): {
                "path": str(path),
                "records": target,
                "bytes": path.stat().st_size,
                "generator": (
                    "real single-end records repeated cyclically without "
                    "modifying alignment fields"
                ),
            }
            for target, path in scaled_inputs.items()
        },
        "measurements": [asdict(row) for row in measurements],
        "summaries": summaries,
        "linear_scaling": {
            "rust": linear_fit(rust_fit_points),
            "htseq": linear_fit(htseq_fit_points),
        },
    }

    json_path = results / "scaling_benchmark.json"
    md_path = results / "scaling_benchmark.md"
    json_path.write_text(json.dumps(report, indent=2) + "\n")
    write_markdown(report, md_path)

    print(
        "\nPASS: every scaling target and repeat produced exactly identical counts."
    )
    print(f"JSON: {json_path}")
    print(f"Markdown: {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
