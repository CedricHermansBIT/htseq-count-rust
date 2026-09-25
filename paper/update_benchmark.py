#!/usr/bin/env python3
"""Update the preprint benchmark snapshot from a TallySeq benchmark JSON file."""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "report",
        nargs="?",
        default="benchmarks/results/real_data_benchmark.json",
    )
    ap.add_argument("--paper-dir", default="paper")
    args = ap.parse_args()

    report = json.loads(Path(args.report).read_text())
    paper = Path(args.paper_dir)
    paper.mkdir(parents=True, exist_ok=True)

    rows = []
    scenario_meta = {x["name"]: x for x in report["scenarios"]}
    for dataset in report["datasets"]:
        for scenario, summary in report["summaries"][dataset].items():
            rs = summary["rust"]
            hs = summary["htseq"]
            meta = scenario_meta[scenario]
            rows.append(
                {
                    "dataset": dataset,
                    "scenario": scenario,
                    "input": meta["input_kind"],
                    "mode": meta["mode"],
                    "stranded": meta["stranded"],
                    "tally_time": rs["median_wall_seconds"],
                    "tally_rss": rs["median_max_rss_kib"] / 1024.0,
                    "htseq_time": hs["median_wall_seconds"],
                    "htseq_rss": hs["median_max_rss_kib"] / 1024.0,
                    "speedup": summary["htseq_over_rust_time_ratio"],
                    "memory_ratio": (
                        hs["median_max_rss_kib"] / rs["median_max_rss_kib"]
                    ),
                }
            )

    header = [
        "dataset",
        "scenario",
        "input",
        "mode",
        "stranded",
        "TallySeq_seconds",
        "TallySeq_peak_RSS_MiB",
        "HTSeq_seconds",
        "HTSeq_peak_RSS_MiB",
        "speedup",
        "memory_ratio_HTSeq_over_TallySeq",
    ]
    with (paper / "benchmark_snapshot.tsv").open("w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write(
                "\t".join(
                    [
                        row["dataset"],
                        row["scenario"],
                        row["input"],
                        row["mode"],
                        row["stranded"],
                        f'{row["tally_time"]:.4f}',
                        f'{row["tally_rss"]:.2f}',
                        f'{row["htseq_time"]:.4f}',
                        f'{row["htseq_rss"]:.2f}',
                        f'{row["speedup"]:.4f}',
                        f'{row["memory_ratio"]:.4f}',
                    ]
                )
                + "\n"
            )

    speedups = [x["speedup"] for x in rows]
    mem_ratios = [x["memory_ratio"] for x in rows]
    tally_times = [x["tally_time"] for x in rows]
    htseq_times = [x["htseq_time"] for x in rows]
    tally_rss = [x["tally_rss"] for x in rows]
    htseq_rss = [x["htseq_rss"] for x in rows]

    macros = f"""\\newcommand{{\\BenchmarkScenarioCount}}{{{len(rows)}}}
\\newcommand{{\\BenchmarkExactScenarioCount}}{{{len(rows)}}}
\\newcommand{{\\MedianSpeedup}}{{{statistics.median(speedups):.1f}}}
\\newcommand{{\\MinimumSpeedup}}{{{min(speedups):.1f}}}
\\newcommand{{\\MaximumSpeedup}}{{{max(speedups):.1f}}}
\\newcommand{{\\MedianTallyTime}}{{{statistics.median(tally_times):.2f}}}
\\newcommand{{\\MedianHTSeqTime}}{{{statistics.median(htseq_times):.2f}}}
\\newcommand{{\\MedianTallyRSS}}{{{statistics.median(tally_rss):.1f}}}
\\newcommand{{\\MedianHTSeqRSS}}{{{statistics.median(htseq_rss):.1f}}}
\\newcommand{{\\MedianMemoryRatio}}{{{statistics.median(mem_ratios):.2f}}}
"""
    (paper / "benchmark_macros.tex").write_text(macros)
    print(f"Wrote {len(rows)} benchmark rows to {paper}")


if __name__ == "__main__":
    main()
