#!/usr/bin/env python3
"""
Small deterministic edge cases that have historically exposed differences
between interval implementations and HTSeq's GenomicArrayOfSets semantics.
"""
from parity_against_htseq import exon, sam, opts

EDGE_CASES = [
    (
        "one_base_feature_exact",
        [exon(100, 100, "geneA")],
        [sam("r1", pos=100, cigar="1M")],
        opts(),
    ),
    (
        "one_base_before_feature",
        [exon(101, 101, "geneA")],
        [sam("r1", pos=100, cigar="1M")],
        opts(),
    ),
    (
        "one_base_after_feature",
        [exon(100, 100, "geneA")],
        [sam("r1", pos=101, cigar="1M")],
        opts(),
    ),
    (
        "adjacent_same_id",
        [exon(100, 104, "geneA"), exon(105, 109, "geneA")],
        [sam("r1", pos=103, cigar="5M")],
        opts(mode="intersection-strict"),
    ),
    (
        "overlapping_same_id",
        [exon(100, 106, "geneA"), exon(104, 110, "geneA")],
        [sam("r1", pos=102, cigar="7M")],
        opts(mode="intersection-strict"),
    ),
    (
        "adjacent_different_ids_union",
        [exon(100, 104, "geneA"), exon(105, 109, "geneB")],
        [sam("r1", pos=104, cigar="2M")],
        opts(mode="union"),
    ),
    (
        "adjacent_different_ids_strict",
        [exon(100, 104, "geneA"), exon(105, 109, "geneB")],
        [sam("r1", pos=104, cigar="2M")],
        opts(mode="intersection-strict"),
    ),
    (
        "contained_overlap_strict",
        [exon(100, 120, "geneA"), exon(105, 110, "geneB")],
        [sam("r1", pos=106, cigar="3M")],
        opts(mode="intersection-strict"),
    ),
    (
        "contained_overlap_nonempty",
        [exon(100, 120, "geneA"), exon(105, 110, "geneB")],
        [sam("r1", pos=103, cigar="10M")],
        opts(mode="intersection-nonempty"),
    ),
    (
        "deletion_crosses_gene_boundary_union",
        [exon(100, 104, "geneA"), exon(108, 112, "geneB")],
        [sam("r1", pos=100, cigar="5M3D5M")],
        opts(mode="union"),
    ),
    (
        "splice_gap_is_not_counted",
        [exon(100, 104, "geneA"), exon(108, 112, "geneB")],
        [sam("r1", pos=100, cigar="5M3N5M")],
        opts(mode="intersection-strict"),
    ),
    (
        "hardclip_padding_do_not_shift",
        [exon(100, 109, "geneA")],
        [sam("r1", pos=100, cigar="2H5M1P5M")],
        opts(),
    ),
    (
        "reverse_unstranded",
        [exon(100, 109, "geneA", strand="+")],
        [sam("r1", flag=16, pos=100, cigar="10M")],
        opts(stranded="no"),
    ),
    (
        "opposite_strand_ignored",
        [exon(100, 109, "geneA", strand="+")],
        [sam("r1", flag=16, pos=100, cigar="10M")],
        opts(stranded="yes"),
    ),
]


if __name__ == "__main__":
    import argparse
    import math
    import shutil
    import subprocess
    import tempfile
    from pathlib import Path

    def run(cmd, cwd):
        return subprocess.run(
            cmd, cwd=cwd, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

    def parse_counts(text):
        out = {}
        for line in text.splitlines():
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            try:
                out[fields[0]] = float(fields[-1])
            except ValueError:
                pass
        return out

    def delta(a, b):
        keys = sorted(set(a) | set(b))
        return [
            (key, a.get(key), b.get(key))
            for key in keys
            if a.get(key) is None
            or b.get(key) is None
            or not math.isclose(a[key], b[key], rel_tol=1e-9, abs_tol=1e-9)
        ]

    ap = argparse.ArgumentParser()
    ap.add_argument("--rust-bin", default="target/release/tallyseq")
    ap.add_argument("--htseq-bin", default="htseq-count")
    args = ap.parse_args()

    root = Path.cwd()
    rust = str((root / args.rust_bin).resolve())
    htseq = shutil.which(args.htseq_bin) or args.htseq_bin
    failures = 0

    with tempfile.TemporaryDirectory(prefix="htseq-rust-edge-") as td:
        td = Path(td)
        for name, gtfs, sams, op in EDGE_CASES:
            case = td / name
            case.mkdir()
            gtf = case / "features.gtf"
            sam_path = case / "reads.sam"
            gtf.write_text("\n".join(gtfs) + "\n")
            sam_path.write_text(
                "@HD\tVN:1.6\tSO:queryname\n"
                "@SQ\tSN:chr1\tLN:1000\n"
                + "\n".join(sams) + "\n"
            )

            rr = run([rust] + op + [str(sam_path), str(gtf)], root)
            hr = run([htseq] + op + [str(sam_path), str(gtf)], root)

            if rr.returncode or hr.returncode:
                failures += 1
                print(f"[DIFF] {name}: exit rust={rr.returncode}, htseq={hr.returncode}")
                if rr.stderr.strip():
                    print("  rust:", rr.stderr.strip().replace("\n", " | "))
                if hr.stderr.strip():
                    print("  htseq:", hr.stderr.strip().replace("\n", " | "))
                continue

            d = delta(parse_counts(rr.stdout), parse_counts(hr.stdout))
            if d:
                failures += 1
                print(f"[DIFF] {name}")
                for row in d:
                    print(" ", row)
            else:
                print(f"[ OK ] {name}")

    print(f"Edge-case differences: {failures}/{len(EDGE_CASES)}")
    raise SystemExit(1 if failures else 0)
