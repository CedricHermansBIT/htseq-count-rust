#!/usr/bin/env python3
"""Benchmark htseq-count-rust against HTSeq on real Pasilla RNA-seq data.

Downloads and verifies:
- Drosophila_melanogaster.BDGP5.78.gtf
- GSM461177_untreat_paired_chr4.bam and/or GSM461178_untreat_paired_chr4.bam

Source: Zenodo 61771, DOI 10.5281/zenodo.61771.

Every measured run is accepted only if the complete count table is exactly
equal between Rust and HTSeq. No numeric tolerance is used.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import statistics
import subprocess
import sys
import urllib.request
from dataclasses import asdict, dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path

RECORD = "61771"
DOI = "10.5281/zenodo.61771"
BASE = "https://zenodo.org/records/61771/files"

FILES = {
    "gtf": (
        "Drosophila_melanogaster.BDGP5.78.gtf",
        "e07b8dcbd4f4f7602feeba0d5ed1698e",
    ),
    "GSM461177": (
        "GSM461177_untreat_paired_chr4.bam",
        "2ebb239b2c5e6b69cdc58c71d9687f83",
    ),
    "GSM461178": (
        "GSM461178_untreat_paired_chr4.bam",
        "443034086a5b6d725b33236c7eb213af",
    ),
}

COMMON_OPTIONS = [
    "-r", "pos",
    "-s", "no",
    "-m", "union",
    "-t", "exon",
    "-i", "gene_id",
    "-a", "10",
    "--nonunique", "none",
    "--secondary-alignments", "ignore",
    "--supplementary-alignments", "ignore",
]


@dataclass
class Measurement:
    tool: str
    dataset: str
    repeat: int
    wall_seconds: float
    max_rss_kib: int
    exit_code: int


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        while True:
            block = fh.read(8 * 1024 * 1024)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def human_bytes(value: int) -> str:
    x = float(value)
    for unit in ("B", "KiB", "MiB", "GiB"):
        if x < 1024 or unit == "GiB":
            return f"{x:.1f} {unit}"
        x /= 1024
    return str(value)


def download_verified(key: str, cache: Path) -> Path:
    name, expected = FILES[key]
    cache.mkdir(parents=True, exist_ok=True)
    dest = cache / name

    if dest.exists():
        actual = md5sum(dest)
        if actual == expected:
            print(f"[cache] {name} ({human_bytes(dest.stat().st_size)})")
            return dest
        print(f"[cache] bad MD5 for {name}; downloading again")
        dest.unlink()

    partial = dest.with_name(dest.name + ".part")
    partial.unlink(missing_ok=True)
    url = f"{BASE}/{name}?download=1"
    print(f"[download] {name}")

    req = urllib.request.Request(
        url, headers={"User-Agent": "htseq-count-rust-benchmark/1.0"}
    )
    try:
        with urllib.request.urlopen(req) as response, partial.open("wb") as out:
            copied = 0
            next_report = 32 * 1024 * 1024
            while True:
                block = response.read(8 * 1024 * 1024)
                if not block:
                    break
                out.write(block)
                copied += len(block)
                if copied >= next_report:
                    print(f"  {human_bytes(copied)}")
                    next_report += 32 * 1024 * 1024
    except Exception:
        partial.unlink(missing_ok=True)
        raise

    actual = md5sum(partial)
    if actual != expected:
        partial.unlink(missing_ok=True)
        raise RuntimeError(
            f"MD5 mismatch for {name}: expected {expected}, got {actual}"
        )
    partial.replace(dest)
    print(f"[verified] {name}: MD5 {expected}")
    return dest


def executable(value: str) -> str:
    if os.path.sep in value:
        path = Path(value).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        return str(path)
    found = shutil.which(value)
    if not found:
        raise FileNotFoundError(f"{value!r} not found on PATH")
    return found


def get_rust_binary(repo: Path, requested: str | None, build: bool) -> Path:
    candidate = (
        Path(requested).expanduser()
        if requested
        else Path("target/release/htseq_count_rust")
    )
    if not candidate.is_absolute():
        candidate = repo / candidate

    if build or not candidate.exists():
        cargo = executable("cargo")
        print("[build] cargo build --release --locked")
        subprocess.run(
            [cargo, "build", "--release", "--locked"], cwd=repo, check=True
        )

    if not candidate.exists():
        raise FileNotFoundError(candidate)
    return candidate.resolve()


def version(command: str, cwd: Path) -> str:
    p = subprocess.run(
        [command, "--version"], cwd=cwd, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    lines = (p.stdout + "\n" + p.stderr).strip().splitlines()
    return lines[0] if lines else f"unknown (exit {p.returncode})"


def counts(path: Path) -> dict[str, Decimal]:
    result = {}
    with path.open() as fh:
        for line_no, raw in enumerate(fh, 1):
            raw = raw.rstrip("\r\n")
            if not raw:
                continue
            fields = raw.split("\t")
            if len(fields) < 2:
                raise RuntimeError(f"{path}:{line_no}: malformed count row")
            try:
                value = Decimal(fields[-1])
            except InvalidOperation as exc:
                raise RuntimeError(
                    f"{path}:{line_no}: invalid count {fields[-1]!r}"
                ) from exc
            if fields[0] in result:
                raise RuntimeError(f"{path}: duplicate feature {fields[0]!r}")
            result[fields[0]] = value
    return result


def exact_differences(a_path: Path, b_path: Path):
    a = counts(a_path)
    b = counts(b_path)
    result = []
    for key in sorted(set(a) | set(b)):
        if a.get(key) != b.get(key):
            result.append((key, a.get(key), b.get(key)))
    return result


def parse_time(path: Path):
    elapsed = rss = None
    with path.open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) != 2:
                continue
            if row[0] == "elapsed_seconds":
                elapsed = float(row[1])
            elif row[0] == "max_rss_kib":
                rss = int(row[1])
    if elapsed is None or rss is None:
        raise RuntimeError(f"Could not parse {path}")
    return elapsed, rss


def timed(tool, dataset, repeat, command, out, err, timing, cwd):
    if not Path("/usr/bin/time").exists():
        raise RuntimeError("GNU /usr/bin/time is required")
    wrapped = [
        "/usr/bin/time", "-f",
        "elapsed_seconds\t%e\nmax_rss_kib\t%M",
        "-o", str(timing), *command,
    ]
    with out.open("w") as stdout, err.open("w") as stderr:
        p = subprocess.run(
            wrapped, cwd=cwd, stdout=stdout, stderr=stderr, check=False
        )
    elapsed, rss = parse_time(timing)
    return Measurement(tool, dataset, repeat, elapsed, rss, p.returncode)


def stderr_tail(path: Path, n=20):
    for line in path.read_text(errors="replace").splitlines()[-n:]:
        print("    " + line)


def summary(rows, tool, dataset):
    selected = [x for x in rows if x.tool == tool and x.dataset == dataset]
    times = [x.wall_seconds for x in selected]
    rss = [x.max_rss_kib for x in selected]
    return {
        "runs": len(selected),
        "median_wall_seconds": statistics.median(times),
        "best_wall_seconds": min(times),
        "median_max_rss_kib": statistics.median(rss),
        "max_rss_kib": max(rss),
    }


def write_markdown(report, path: Path):
    lines = [
        "# Real-data benchmark",
        "",
        f"Source: Zenodo DOI {DOI}",
        f"HTSeq version: {report['htseq_version']}",
        f"Rust version: {report['rust_version']}",
        "",
        "Every reported run passed exact count equality.",
        "",
        "| Dataset | Tool | Runs | Median time (s) | Best time (s) | Median RSS (MiB) |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for dataset in report["datasets"]:
        for tool in ("rust", "htseq"):
            s = report["summaries"][dataset][tool]
            lines.append(
                f"| {dataset} | {tool} | {s['runs']} | "
                f"{s['median_wall_seconds']:.3f} | {s['best_wall_seconds']:.3f} | "
                f"{s['median_max_rss_kib'] / 1024:.1f} |"
            )
        ratio = report["summaries"][dataset]["htseq_over_rust_time_ratio"]
        lines.append(
            f"| {dataset} | HTSeq/Rust median ratio |  | {ratio:.2f}x |  |  |"
        )
    lines += [
        "",
        "Count options:",
        "",
        "    " + " ".join(COMMON_OPTIONS),
        "",
        "Downloads are cached under benchmarks/data and MD5-verified.",
        "",
    ]
    path.write_text("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".")
    ap.add_argument(
        "--dataset", action="append",
        choices=("GSM461177", "GSM461178"),
        help="May be repeated. Default: GSM461177."
    )
    ap.add_argument("--cache-dir", default="benchmarks/data")
    ap.add_argument("--results-dir", default="benchmarks/results")
    ap.add_argument("--rust-bin")
    ap.add_argument("--htseq-bin", default="htseq-count")
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument(
        "--warmup", action="store_true",
        help="Run an unmeasured exact-parity pass first."
    )
    args = ap.parse_args()

    if args.repeats < 1:
        ap.error("--repeats must be >= 1")

    repo = Path(args.repo).expanduser().resolve()
    cache = Path(args.cache_dir)
    results = Path(args.results_dir)
    if not cache.is_absolute():
        cache = repo / cache
    if not results.is_absolute():
        results = repo / results
    results.mkdir(parents=True, exist_ok=True)
    run_dir = results / "runs"
    run_dir.mkdir(parents=True, exist_ok=True)

    datasets = args.dataset or ["GSM461177"]
    rust = get_rust_binary(repo, args.rust_bin, args.build)
    htseq = executable(args.htseq_bin)

    print(f"Rust:  {rust}")
    print(f"HTSeq: {htseq} ({version(htseq, repo)})")
    print(f"Data:  Zenodo DOI {DOI}")

    gtf = download_verified("gtf", cache)
    bams = {key: download_verified(key, cache) for key in datasets}

    def tool_commands(dataset):
        tail = [str(bams[dataset]), str(gtf)]
        return (
            [str(rust), *COMMON_OPTIONS, *tail],
            [htseq, *COMMON_OPTIONS, *tail],
        )

    if args.warmup:
        for dataset in datasets:
            print(f"\n[warmup] {dataset}")
            rust_cmd, htseq_cmd = tool_commands(dataset)
            rust_out = run_dir / f"{dataset}.warmup.rust.tsv"
            htseq_out = run_dir / f"{dataset}.warmup.htseq.tsv"
            with rust_out.open("w") as out:
                rp = subprocess.run(
                    rust_cmd, cwd=repo, stdout=out,
                    stderr=subprocess.DEVNULL
                )
            with htseq_out.open("w") as out:
                hp = subprocess.run(
                    htseq_cmd, cwd=repo, stdout=out,
                    stderr=subprocess.DEVNULL
                )
            if rp.returncode or hp.returncode:
                raise RuntimeError(
                    f"Warmup failed: rust={rp.returncode}, htseq={hp.returncode}"
                )
            diff = exact_differences(rust_out, htseq_out)
            if diff:
                for row in diff[:50]:
                    print("  ", row)
                raise RuntimeError(
                    f"{dataset}: {len(diff)} exact count differences in warmup"
                )
            print(f"  [exact] {len(counts(rust_out))} rows identical")
            rust_out.unlink()
            htseq_out.unlink()

    measurements = []

    for dataset in datasets:
        print(f"\n[benchmark] {dataset}")
        rust_cmd, htseq_cmd = tool_commands(dataset)

        for repeat in range(1, args.repeats + 1):
            # Alternate order so page-cache effects do not always help one tool.
            order = (
                [("rust", rust_cmd), ("htseq", htseq_cmd)]
                if repeat % 2
                else [("htseq", htseq_cmd), ("rust", rust_cmd)]
            )
            outputs = {}

            for tool, command in order:
                stem = f"{dataset}.r{repeat}.{tool}"
                out = run_dir / f"{stem}.counts.tsv"
                err = run_dir / f"{stem}.stderr.txt"
                timing = run_dir / f"{stem}.time.tsv"
                m = timed(
                    tool, dataset, repeat, command, out, err, timing, repo
                )
                measurements.append(m)
                outputs[tool] = out
                print(
                    f"  repeat {repeat} {tool:5s}: "
                    f"{m.wall_seconds:.3f}s, {m.max_rss_kib / 1024:.1f} MiB RSS"
                )
                if m.exit_code:
                    print(f"[FAIL] {tool} exited {m.exit_code}")
                    stderr_tail(err)
                    return 2

            diff = exact_differences(outputs["rust"], outputs["htseq"])
            if diff:
                print(
                    f"[FAIL] {dataset} repeat {repeat}: "
                    f"{len(diff)} exact count differences"
                )
                for feature, rv, hv in diff[:50]:
                    print(f"  {feature}: rust={rv}, htseq={hv}")
                if len(diff) > 50:
                    print(f"  ... {len(diff) - 50} more")
                return 3
            print(f"  [exact] {len(counts(outputs['rust']))} rows identical")

    report = {
        "schema_version": 1,
        "zenodo_record": RECORD,
        "zenodo_doi": DOI,
        "datasets": datasets,
        "count_options": COMMON_OPTIONS,
        "exact_count_equality": True,
        "rust_binary": str(rust),
        "rust_version": version(str(rust), repo),
        "htseq_binary": htseq,
        "htseq_version": version(htseq, repo),
        "system": {
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "cpu_count": os.cpu_count(),
        },
        "files": {
            "gtf": {
                "path": str(gtf),
                "bytes": gtf.stat().st_size,
                "md5": FILES["gtf"][1],
            },
            **{
                d: {
                    "path": str(bams[d]),
                    "bytes": bams[d].stat().st_size,
                    "md5": FILES[d][1],
                }
                for d in datasets
            },
        },
        "measurements": [asdict(x) for x in measurements],
        "summaries": {},
    }

    for dataset in datasets:
        rs = summary(measurements, "rust", dataset)
        hs = summary(measurements, "htseq", dataset)
        report["summaries"][dataset] = {
            "rust": rs,
            "htseq": hs,
            "htseq_over_rust_time_ratio":
                hs["median_wall_seconds"] / rs["median_wall_seconds"],
        }

    json_path = results / "real_data_benchmark.json"
    md_path = results / "real_data_benchmark.md"
    json_path.write_text(json.dumps(report, indent=2) + "\n")
    write_markdown(report, md_path)

    print("\nPASS: all measured count tables were exactly identical.")
    print(f"JSON: {json_path}")
    print(f"Markdown: {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
