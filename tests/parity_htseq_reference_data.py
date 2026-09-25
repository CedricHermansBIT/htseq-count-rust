#!/usr/bin/env python3
import argparse
import gzip
import math
import shutil
import subprocess
import tempfile
from pathlib import Path

SPECIAL = {
    "__no_feature",
    "__ambiguous",
    "__too_low_aQual",
    "__not_aligned",
    "__alignment_not_unique",
}

def run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=False)

def parse_counts(text):
    counts = {}
    for line in text.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue
        try:
            counts[fields[0]] = float(fields[-1])
        except ValueError:
            pass
    return counts

def compare(a, b):
    keys = sorted(set(a) | set(b))
    return [(key, a.get(key), b.get(key)) for key in keys
            if a.get(key) is None or b.get(key) is None
            or not math.isclose(a[key], b[key], rel_tol=1e-9, abs_tol=1e-9)]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-dir", required=True)
    ap.add_argument("--rust-bin", default="target/release/tallyseq")
    ap.add_argument("--htseq-bin", default="htseq-count")
    args = ap.parse_args()

    root = Path.cwd()
    ref = Path(args.reference_dir).resolve() / "example_data"
    rust = (root / args.rust_bin).resolve()
    htseq = shutil.which(args.htseq_bin) or args.htseq_bin

    fixtures = []
    fixtures.append((
        "upstream_bamfile_no_qualities_default",
        ref / "bamfile_no_qualities.sam",
        ref / "bamfile_no_qualities.gtf",
        [],
    ))

    fixtures.append((
        "upstream_paired_position_sorted",
        ref / "SRR001432_head_sorted.bam",
        ref / "bamfile_no_qualities.gtf",
        ["-r", "pos", "-s", "no"],
    ))

    yeast_gz = ref / "Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz"
    with tempfile.TemporaryDirectory(prefix="htseq-ref-") as td:
        yeast_gtf = Path(td) / "yeast.gtf"
        with gzip.open(yeast_gz, "rt") as src, yeast_gtf.open("w") as dst:
            shutil.copyfileobj(src, dst)

        yeast_sam = ref / "yeast_RNASeq_excerpt_withNH.sam"
        common_score = [
            "-m", "intersection-nonempty",
            "--secondary-alignments", "score",
            "--supplementary-alignments", "score",
        ]
        fixtures.extend([
            (
                "upstream_yeast_nonunique_none",
                yeast_sam,
                yeast_gtf,
                common_score + ["--nonunique", "none"],
            ),
            (
                "upstream_yeast_nonunique_all",
                yeast_sam,
                yeast_gtf,
                common_score + ["--nonunique", "all"],
            ),
            (
                "upstream_yeast_nonunique_fraction",
                yeast_sam,
                yeast_gtf,
                common_score + ["--nonunique", "fraction"],
            ),
            (
                "upstream_yeast_ignore_secondary",
                yeast_sam,
                yeast_gtf,
                [
                    "-m", "intersection-nonempty",
                    "--nonunique", "none",
                    "--secondary-alignments", "ignore",
                    "--supplementary-alignments", "score",
                ],
            ),
        ])

        failures = 0
        for name, sam, gtf, opts in fixtures:
            rust_cmd = [str(rust)] + opts + [str(sam), str(gtf)]
            htseq_cmd = [htseq] + opts + [str(sam), str(gtf)]
            rr = run(rust_cmd, root)
            hr = run(htseq_cmd, root)

            if rr.returncode or hr.returncode:
                failures += 1
                print(f"[DIFF] {name}: exit rust={rr.returncode}, htseq={hr.returncode}")
                if rr.stderr.strip():
                    print("  rust stderr tail:", " | ".join(rr.stderr.splitlines()[-5:]))
                if hr.stderr.strip():
                    print("  htseq stderr tail:", " | ".join(hr.stderr.splitlines()[-5:]))
                continue

            diff = compare(parse_counts(rr.stdout), parse_counts(hr.stdout))
            if diff:
                failures += 1
                print(f"[DIFF] {name}: {len(diff)} differing rows")
                for row in diff[:25]:
                    print(" ", row)
                if len(diff) > 25:
                    print(f"  ... {len(diff) - 25} more")
            else:
                print(f"[ OK ] {name}")

        print(f"Official HTSeq fixture differences: {failures}/{len(fixtures)}")
        raise SystemExit(1 if failures else 0)

if __name__ == "__main__":
    main()
