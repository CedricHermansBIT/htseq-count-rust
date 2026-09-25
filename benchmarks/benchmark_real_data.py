#!/usr/bin/env python3
"""Real-data differential benchmark for TallySeq vs HTSeq.

The benchmark downloads real Pasilla RNA-seq data and its matching Ensembl GTF
from Zenodo record 61771, verifies published MD5 checksums, and runs a matrix of
single-end and paired-end counting scenarios.

Correctness is a hard gate. For every scenario and every repeat, feature IDs,
special rows, and numeric count values must be exactly equal between Rust and
HTSeq. No floating-point tolerance is used. Performance results are only
reported for scenarios that pass exact count parity.

The source BAM is paired-end. A reproducible single-end BAM is derived from
mate 1 of the same real alignments by clearing pairing metadata; coordinates,
CIGAR strings, MAPQ values, tags, strands, and read sequences remain real.
A query-name sorted paired BAM is also derived automatically for --order name.
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
from typing import Iterable

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

MODES = ("union", "intersection-strict", "intersection-nonempty")
STRANDEDNESS = ("no", "yes", "reverse")


@dataclass(frozen=True)
class Scenario:
    name: str
    input_kind: str
    mode: str = "union"
    stranded: str = "no"
    order: str | None = None
    feature_types: tuple[str, ...] = ("exon",)
    idattrs: tuple[str, ...] = ("gene_id",)
    minaqual: int = 10
    nonunique: str = "none"
    secondary: str = "ignore"
    supplementary: str = "ignore"
    category: str = "core"

    def options(self) -> list[str]:
        out: list[str] = []
        if self.input_kind == "paired":
            out += ["-r", self.order or "pos"]
        out += [
            "-s", self.stranded,
            "-m", self.mode,
            "-a", str(self.minaqual),
            "--nonunique", self.nonunique,
            "--secondary-alignments", self.secondary,
            "--supplementary-alignments", self.supplementary,
        ]
        for feature_type in self.feature_types:
            out += ["-t", feature_type]
        for attr in self.idattrs:
            out += ["-i", attr]
        return out


@dataclass
class Measurement:
    tool: str
    dataset: str
    scenario: str
    repeat: int
    wall_seconds: float
    max_rss_kib: int
    exit_code: int
    gtf_parse_seconds: float | None = None
    index_build_seconds: float | None = None
    annotation_total_seconds: float | None = None
    counting_seconds: float | None = None


def scenario_matrix(profile: str) -> list[Scenario]:
    scenarios: list[Scenario] = []

    # Core semantic matrix:
    # single + paired x all overlap modes x all strandedness modes = 18 cases.
    for input_kind in ("single", "paired"):
        for mode in MODES:
            for stranded in STRANDEDNESS:
                scenarios.append(
                    Scenario(
                        name=f"{input_kind}_{mode}_{stranded}",
                        input_kind=input_kind,
                        mode=mode,
                        stranded=stranded,
                        order="pos" if input_kind == "paired" else None,
                        category="core",
                    )
                )

    if profile == "core":
        return scenarios

    # Target options are intentionally not a full Cartesian product. Each one
    # isolates a behavior that can change HTSeq counts or paired-end semantics.
    targeted = [
        Scenario(
            "single_union_no_nonunique_all",
            "single",
            nonunique="all",
            category="options",
        ),
        Scenario(
            "single_union_no_nonunique_fraction",
            "single",
            nonunique="fraction",
            category="options",
        ),
        Scenario(
            "paired_union_no_nonunique_all",
            "paired",
            order="pos",
            nonunique="all",
            category="options",
        ),
        Scenario(
            "paired_union_no_nonunique_fraction",
            "paired",
            order="pos",
            nonunique="fraction",
            category="options",
        ),
        Scenario(
            "single_union_no_minaqual_0",
            "single",
            minaqual=0,
            category="options",
        ),
        Scenario(
            "paired_union_no_minaqual_0",
            "paired",
            order="pos",
            minaqual=0,
            category="options",
        ),
        Scenario(
            "single_union_no_secondary_supplementary_score",
            "single",
            secondary="score",
            supplementary="score",
            category="options",
        ),
        Scenario(
            "paired_union_no_secondary_supplementary_score",
            "paired",
            order="pos",
            secondary="score",
            supplementary="score",
            category="options",
        ),
        Scenario(
            "single_union_no_multi_type_exon_cds",
            "single",
            feature_types=("exon", "CDS"),
            category="options",
        ),
        Scenario(
            "paired_union_no_multi_type_exon_cds",
            "paired",
            order="pos",
            feature_types=("exon", "CDS"),
            category="options",
        ),
        Scenario(
            "single_union_no_multi_id_gene_transcript",
            "single",
            idattrs=("gene_id", "transcript_id"),
            category="options",
        ),
        Scenario(
            "paired_union_no_multi_id_gene_transcript",
            "paired",
            order="pos",
            idattrs=("gene_id", "transcript_id"),
            category="options",
        ),
        Scenario(
            "paired_union_no_name_sorted",
            "paired",
            order="name",
            category="options",
        ),
    ]
    scenarios.extend(targeted)
    return scenarios


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
        url, headers={"User-Agent": "tallyseq-benchmark/2.0"}
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
        else Path("target/release/tallyseq")
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
        [command, "--version"],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = (p.stdout + "\n" + p.stderr).strip().splitlines()
    return lines[0] if lines else f"unknown (exit {p.returncode})"


def require_pysam():
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError(
            "pysam is required to derive the benchmark inputs. "
            "Installing HTSeq with pip installs pysam as a dependency."
        ) from exc
    return pysam


def derive_single_end(source: Path, destination: Path) -> Path:
    """Create a real single-end BAM from mate 1 of the paired source BAM."""
    pysam = require_pysam()
    fingerprint = destination.with_suffix(destination.suffix + ".source.md5")
    source_md5 = md5sum(source)

    if (
        destination.exists()
        and fingerprint.exists()
        and fingerprint.read_text().strip() == source_md5
    ):
        print(f"[cache] {destination.name} ({human_bytes(destination.stat().st_size)})")
        return destination

    tmp = destination.with_name(destination.name + ".part")
    tmp.unlink(missing_ok=True)

    with pysam.AlignmentFile(str(source), "rb") as reader:
        with pysam.AlignmentFile(str(tmp), "wb", header=reader.header) as writer:
            for record in reader.fetch(until_eof=True):
                if record.is_paired and not record.is_read1:
                    continue

                # Preserve the real alignment, sequence, tags, CIGAR, MAPQ and
                # strand. Remove only pair-specific SAM metadata so both tools
                # treat this as single-end input.
                record.flag &= ~(
                    0x1   # PAIRED
                    | 0x2 # PROPER_PAIR
                    | 0x8 # MATE_UNMAPPED
                    | 0x20 # MATE_REVERSE
                    | 0x40 # READ1
                    | 0x80 # READ2
                )
                record.next_reference_id = -1
                record.next_reference_start = -1
                record.template_length = 0
                writer.write(record)

    tmp.replace(destination)
    fingerprint.write_text(source_md5 + "\n")
    print(f"[derived] {destination.name} from real mate-1 alignments")
    return destination


def derive_name_sorted(source: Path, destination: Path) -> Path:
    pysam = require_pysam()
    fingerprint = destination.with_suffix(destination.suffix + ".source.md5")
    source_md5 = md5sum(source)

    if (
        destination.exists()
        and fingerprint.exists()
        and fingerprint.read_text().strip() == source_md5
    ):
        print(f"[cache] {destination.name} ({human_bytes(destination.stat().st_size)})")
        return destination

    tmp = destination.with_name(destination.name + ".part")
    tmp.unlink(missing_ok=True)
    pysam.sort("-n", "-o", str(tmp), str(source))
    tmp.replace(destination)
    fingerprint.write_text(source_md5 + "\n")
    print(f"[derived] {destination.name} query-name sorted")
    return destination


def inspect_bam(path: Path) -> dict:
    pysam = require_pysam()
    stats = {
        "records": 0,
        "paired_records": 0,
        "secondary_records": 0,
        "supplementary_records": 0,
        "unmapped_records": 0,
        "nh_multimapped_records": 0,
    }
    with pysam.AlignmentFile(str(path), "rb") as reader:
        for record in reader.fetch(until_eof=True):
            stats["records"] += 1
            stats["paired_records"] += int(record.is_paired)
            stats["secondary_records"] += int(record.is_secondary)
            stats["supplementary_records"] += int(record.is_supplementary)
            stats["unmapped_records"] += int(record.is_unmapped)
            if record.has_tag("NH") and record.get_tag("NH") > 1:
                stats["nh_multimapped_records"] += 1
    return stats


def counts(path: Path) -> dict[str, Decimal]:
    result: dict[str, Decimal] = {}
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
    """Compare normalized tables with no numeric tolerance."""
    a = counts(a_path)
    b = counts(b_path)
    result = []
    for key in sorted(set(a) | set(b)):
        if a.get(key) != b.get(key):
            result.append((key, a.get(key), b.get(key)))
    return result


def parse_phase_timings(path: Path) -> dict[str, float]:
    timings: dict[str, float] = {}
    prefixes = {
        "__timing_gtf_parse_seconds": "gtf_parse_seconds",
        "__timing_index_build_seconds": "index_build_seconds",
        "__timing_annotation_total_seconds": "annotation_total_seconds",
        "__timing_counting_seconds": "counting_seconds",
    }
    for line in path.read_text(errors="replace").splitlines():
        for prefix, key in prefixes.items():
            if line.startswith(prefix + "\t"):
                try:
                    timings[key] = float(line.split("\t", 1)[1])
                except ValueError:
                    pass
    return timings


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


def timed(
    tool: str,
    dataset: str,
    scenario: str,
    repeat: int,
    command: list[str],
    out: Path,
    err: Path,
    timing: Path,
    cwd: Path,
) -> Measurement:
    if not Path("/usr/bin/time").exists():
        raise RuntimeError("GNU /usr/bin/time is required")
    wrapped = [
        "/usr/bin/time",
        "-f",
        "elapsed_seconds\t%e\nmax_rss_kib\t%M",
        "-o",
        str(timing),
        *command,
    ]
    env = os.environ.copy()
    if tool == "rust":
        env["HTSEQ_COUNT_RUST_TIMINGS"] = "1"

    with out.open("w") as stdout, err.open("w") as stderr:
        p = subprocess.run(
            wrapped,
            cwd=cwd,
            stdout=stdout,
            stderr=stderr,
            check=False,
            env=env,
        )

    elapsed, rss = parse_time(timing)
    phases = parse_phase_timings(err) if tool == "rust" else {}
    return Measurement(
        tool=tool,
        dataset=dataset,
        scenario=scenario,
        repeat=repeat,
        wall_seconds=elapsed,
        max_rss_kib=rss,
        exit_code=p.returncode,
        **phases,
    )


def stderr_tail(path: Path, n: int = 20):
    for line in path.read_text(errors="replace").splitlines()[-n:]:
        print("    " + line)


def summarize(
    rows: Iterable[Measurement],
    tool: str,
    dataset: str,
    scenario: str,
):
    selected = [
        x
        for x in rows
        if x.tool == tool
        and x.dataset == dataset
        and x.scenario == scenario
    ]
    times = [x.wall_seconds for x in selected]
    rss = [x.max_rss_kib for x in selected]
    return {
        "runs": len(selected),
        "median_wall_seconds": statistics.median(times),
        "best_wall_seconds": min(times),
        "median_max_rss_kib": statistics.median(rss),
        "max_rss_kib": max(rss),
    }


def write_markdown(report: dict, path: Path):
    scenarios_by_name = {
        row["name"]: row for row in report["scenarios"]
    }
    lines = [
        "# Real-data benchmark",
        "",
        f"Source: Zenodo DOI {DOI}",
        f"HTSeq version: {report['htseq_version']}",
        f"Rust version: {report['rust_version']}",
        "",
        "Every reported run passed exact normalized count equality. "
        "No numeric tolerance was used.",
        "",
        "| Dataset | Scenario | Input | Mode | Strand | Tool | Runs | Median time (s) | Best (s) | Median RSS (MiB) |",
        "| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: |",
    ]

    for dataset in report["datasets"]:
        for scenario_name, scenario_summary in report["summaries"][dataset].items():
            scenario = scenarios_by_name[scenario_name]
            for tool in ("rust", "htseq"):
                s = scenario_summary[tool]
                lines.append(
                    f"| {dataset} | {scenario_name} | "
                    f"{scenario['input_kind']} | {scenario['mode']} | "
                    f"{scenario['stranded']} | {tool} | {s['runs']} | "
                    f"{s['median_wall_seconds']:.3f} | "
                    f"{s['best_wall_seconds']:.3f} | "
                    f"{s['median_max_rss_kib'] / 1024:.1f} |"
                )

    lines += [
        "",
        "## Coverage",
        "",
        "- Core matrix: single-end and paired-end x union, intersection-strict, and intersection-nonempty x no/yes/reverse strandedness.",
        "- Targeted options: nonunique all/fraction, MAPQ 0, secondary/supplementary score, repeated feature types, repeated ID attributes, and paired position/name ordering.",
        "- nonunique=random is intentionally excluded from exact parity benchmarking because HTSeq and Rust use independent random-number generators, so exact feature assignment is not deterministic across implementations.",
        "- Paired --samout is not benchmarked because paired SAM annotation is not implemented in TallySeq yet.",
        "",
        "The single-end BAM is reproducibly derived from read 1 of the real paired-end BAM by removing only pairing metadata.",
        "",
    ]
    path.write_text("\n".join(lines))


def input_for_scenario(
    scenario: Scenario,
    paired_pos: Path,
    paired_name: Path,
    single: Path,
) -> Path:
    if scenario.input_kind == "single":
        return single
    if scenario.order == "name":
        return paired_name
    return paired_pos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".")
    ap.add_argument(
        "--dataset",
        action="append",
        choices=("GSM461177", "GSM461178"),
        help="May be repeated. Default: GSM461177.",
    )
    ap.add_argument("--cache-dir", default="benchmarks/data")
    ap.add_argument("--results-dir", default="benchmarks/results")
    ap.add_argument("--rust-bin")
    ap.add_argument("--htseq-bin", default="htseq-count")
    ap.add_argument("--build", action="store_true")
    ap.add_argument(
        "--profile",
        choices=("core", "full"),
        default="full",
        help="core = 18 SE/PE mode/strand cases; full adds targeted option cases.",
    )
    ap.add_argument("--repeats", type=int, default=1)
    ap.add_argument(
        "--scenario",
        action="append",
        help="Run only the named scenario. May be repeated.",
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
    cache.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    run_dir = results / "runs"
    run_dir.mkdir(parents=True, exist_ok=True)

    datasets = args.dataset or ["GSM461177"]
    rust = get_rust_binary(repo, args.rust_bin, args.build)
    htseq = executable(args.htseq_bin)

    scenarios = scenario_matrix(args.profile)
    if args.scenario:
        wanted = set(args.scenario)
        scenarios = [s for s in scenarios if s.name in wanted]
        missing = wanted - {s.name for s in scenarios}
        if missing:
            ap.error("unknown scenario(s): " + ", ".join(sorted(missing)))

    print(f"Rust:     {rust}")
    print(f"HTSeq:    {htseq} ({version(htseq, repo)})")
    print(f"Data:     Zenodo DOI {DOI}")
    print(f"Profile:  {args.profile}")
    print(f"Scenarios:{len(scenarios)}")

    gtf = download_verified("gtf", cache)
    paired_pos = {
        key: download_verified(key, cache) for key in datasets
    }

    inputs: dict[str, dict[str, Path]] = {}
    input_stats: dict[str, dict[str, dict]] = {}

    for dataset in datasets:
        source = paired_pos[dataset]
        single = derive_single_end(
            source, cache / f"{dataset}.mate1.single.bam"
        )
        paired_name = derive_name_sorted(
            source, cache / f"{dataset}.queryname.bam"
        )
        inputs[dataset] = {
            "single": single,
            "paired_pos": source,
            "paired_name": paired_name,
        }
        input_stats[dataset] = {
            "single": inspect_bam(single),
            "paired_pos": inspect_bam(source),
            "paired_name": inspect_bam(paired_name),
        }

    measurements: list[Measurement] = []

    for dataset in datasets:
        print(f"\n[dataset] {dataset}")
        for scenario in scenarios:
            bam = input_for_scenario(
                scenario,
                inputs[dataset]["paired_pos"],
                inputs[dataset]["paired_name"],
                inputs[dataset]["single"],
            )
            options = scenario.options()
            rust_cmd = [str(rust), *options, str(bam), str(gtf)]
            htseq_cmd = [htseq, *options, str(bam), str(gtf)]

            print(
                f"\n[scenario] {scenario.name} "
                f"({scenario.input_kind}, {scenario.mode}, "
                f"strand={scenario.stranded})"
            )

            for repeat in range(1, args.repeats + 1):
                # Alternate execution order so page cache does not
                # systematically favor the second executable.
                order = (
                    [("rust", rust_cmd), ("htseq", htseq_cmd)]
                    if repeat % 2
                    else [("htseq", htseq_cmd), ("rust", rust_cmd)]
                )
                outputs: dict[str, Path] = {}

                for tool, command in order:
                    stem = (
                        f"{dataset}.{scenario.name}.r{repeat}.{tool}"
                    )
                    out = run_dir / f"{stem}.counts.tsv"
                    err = run_dir / f"{stem}.stderr.txt"
                    timing = run_dir / f"{stem}.time.tsv"
                    m = timed(
                        tool,
                        dataset,
                        scenario.name,
                        repeat,
                        command,
                        out,
                        err,
                        timing,
                        repo,
                    )
                    measurements.append(m)
                    outputs[tool] = out
                    print(
                        f"  repeat {repeat} {tool:5s}: "
                        f"{m.wall_seconds:.3f}s, "
                        f"{m.max_rss_kib / 1024:.1f} MiB RSS"
                    )
                    if tool == "rust" and m.annotation_total_seconds is not None:
                        print(
                            "    phases: "
                            f"gtf_parse={m.gtf_parse_seconds:.3f}s, "
                            f"index={m.index_build_seconds:.3f}s, "
                            f"count={m.counting_seconds:.3f}s"
                        )
                    if m.exit_code:
                        print(f"[FAIL] {tool} exited {m.exit_code}")
                        stderr_tail(err)
                        return 2

                diff = exact_differences(
                    outputs["rust"], outputs["htseq"]
                )
                if diff:
                    print(
                        f"[FAIL] {dataset} / {scenario.name} / "
                        f"repeat {repeat}: {len(diff)} exact count "
                        "differences"
                    )
                    for feature, rv, hv in diff[:50]:
                        print(
                            f"  {feature}: rust={rv}, htseq={hv}"
                        )
                    if len(diff) > 50:
                        print(f"  ... {len(diff) - 50} more")
                    return 3

                print(
                    f"  [exact] "
                    f"{len(counts(outputs['rust']))} count rows identical"
                )

    report = {
        "schema_version": 2,
        "zenodo_record": RECORD,
        "zenodo_doi": DOI,
        "datasets": datasets,
        "profile": args.profile,
        "exact_count_equality": True,
        "comparison": (
            "Exact feature-key and Decimal numeric-value equality; "
            "no tolerance. Row order is normalized."
        ),
        "excluded_from_exact_matrix": [
            "nonunique=random: independent RNGs make exact ambiguous-feature "
            "assignment nondeterministic across implementations",
            "paired samout: not yet implemented in TallySeq",
        ],
        "rust_binary": str(rust),
        "rust_version": version(str(rust), repo),
        "htseq_binary": htseq,
        "htseq_version": version(htseq, repo),
        "system": {
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "cpu_count": os.cpu_count(),
        },
        "scenarios": [asdict(s) for s in scenarios],
        "input_stats": input_stats,
        "files": {
            "gtf": {
                "path": str(gtf),
                "bytes": gtf.stat().st_size,
                "md5": FILES["gtf"][1],
            },
            **{
                dataset: {
                    "source_paired_bam": str(inputs[dataset]["paired_pos"]),
                    "single_derived_bam": str(inputs[dataset]["single"]),
                    "name_sorted_paired_bam": str(
                        inputs[dataset]["paired_name"]
                    ),
                    "source_bytes": inputs[dataset][
                        "paired_pos"
                    ].stat().st_size,
                    "source_md5": FILES[dataset][1],
                }
                for dataset in datasets
            },
        },
        "measurements": [asdict(x) for x in measurements],
        "summaries": {},
    }

    for dataset in datasets:
        report["summaries"][dataset] = {}
        for scenario in scenarios:
            rs = summarize(
                measurements, "rust", dataset, scenario.name
            )
            hs = summarize(
                measurements, "htseq", dataset, scenario.name
            )
            report["summaries"][dataset][scenario.name] = {
                "rust": rs,
                "htseq": hs,
                "htseq_over_rust_time_ratio": (
                    hs["median_wall_seconds"]
                    / rs["median_wall_seconds"]
                ),
            }

    json_path = results / "real_data_benchmark.json"
    md_path = results / "real_data_benchmark.md"
    json_path.write_text(json.dumps(report, indent=2) + "\n")
    write_markdown(report, md_path)

    print(
        "\nPASS: every scenario and repeat produced exactly identical counts."
    )
    print(f"JSON: {json_path}")
    print(f"Markdown: {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
