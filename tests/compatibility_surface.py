#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import tempfile
from decimal import Decimal
from pathlib import Path

import numpy as np
import pysam


def run(cmd, cwd, *, input_text=None, env=None):
    return subprocess.run(
        [str(x) for x in cmd],
        cwd=cwd,
        text=True,
        input=input_text,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
        check=False,
    )


def require_ok(label, result):
    if result.returncode != 0:
        raise RuntimeError(
            f"{label} exited {result.returncode}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def parse_table(text):
    rows = {}
    header = None
    for raw in text.splitlines():
        fields = raw.split("\t")
        if fields and fields[0] == "":
            header = fields
            continue
        if not fields or len(fields) < 2:
            continue
        rows[fields[0]] = fields[1:]
    return header, rows


def numeric_tail(fields, n):
    return [Decimal(x) for x in fields[-n:]]


def compare_tabular(label, rust_text, htseq_text, n_samples, metadata_cols=0):
    rust_header, rust_rows = parse_table(rust_text)
    ht_header, ht_rows = parse_table(htseq_text)

    if set(rust_rows) != set(ht_rows):
        raise AssertionError(
            f"{label}: row IDs differ\n"
            f"rust-only={sorted(set(rust_rows)-set(ht_rows))}\n"
            f"htseq-only={sorted(set(ht_rows)-set(rust_rows))}"
        )

    for key in sorted(rust_rows):
        rr = rust_rows[key]
        hr = ht_rows[key]
        if rr[:metadata_cols] != hr[:metadata_cols]:
            raise AssertionError(
                f"{label}: metadata differs for {key}: {rr[:metadata_cols]} != {hr[:metadata_cols]}"
            )
        if numeric_tail(rr, n_samples) != numeric_tail(hr, n_samples):
            raise AssertionError(
                f"{label}: counts differ for {key}: "
                f"{numeric_tail(rr,n_samples)} != {numeric_tail(hr,n_samples)}"
            )
    return rust_header, ht_header


def write_fixture(root: Path):
    gtf = root / "features.gtf"
    gtf.write_text(
        'chr1\tcompat\texon\t100\t149\t.\t+\t.\t'
        'gene_id "geneA"; gene_name "Alpha"; transcript_id "txA"; note "semi;colon";\n'
        'chr1\tcompat\texon\t200\t249\t.\t-\t.\t'
        'gene_id "geneB"; gene_name "Beta"; transcript_id "txB"; note "plain";\n'
    )

    sam1 = root / "reads1.sam"
    sam1.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n"
        "@SQ\tSN:chr1\tLN:1000\n"
        "a1\t0\tchr1\t105\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tNH:i:1\n"
        "b1\t16\tchr1\t205\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tNH:i:1\n"
        "amb1\t0\tchr1\t145\t60\t60M\t*\t0\t0\t"
        + "A"*60 + "\t" + "I"*60 + "\tNH:i:1\n"
    )

    sam2 = root / "reads2.sam"
    sam2.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n"
        "@SQ\tSN:chr1\tLN:1000\n"
        "a2\t0\tchr1\t110\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tNH:i:1\n"
        "a3\t0\tchr1\t120\t0\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tNH:i:1\n"
    )

    paired = root / "paired.sam"
    paired.write_text(
        "@HD\tVN:1.6\tSO:queryname\n"
        "@SQ\tSN:chr1\tLN:1000\n"
        "p1\t99\tchr1\t105\t60\t10M\t=\t125\t30\tAAAAAAAAAA\tIIIIIIIIII\tNH:i:1\n"
        "p1\t147\tchr1\t125\t60\t10M\t=\t105\t-30\tAAAAAAAAAA\tIIIIIIIIII\tNH:i:1\n"
        "p2\t99\tchr1\t205\t60\t10M\t=\t225\t30\tAAAAAAAAAA\tIIIIIIIIII\tNH:i:1\n"
        "p2\t147\tchr1\t225\t60\t10M\t=\t205\t-30\tAAAAAAAAAA\tIIIIIIIIII\tNH:i:1\n"
    )

    return gtf, sam1, sam2, paired


def test_cli_compat(repo, rust, htseq, gtf, sam1):
    version = run([rust, "--version"], repo)
    require_ok("Rust --version", version)
    if not version.stdout.strip():
        raise AssertionError("--version produced no output")
    print("[ OK ] version")

    # HTSeq accepts -f for old scripts but ignores it in modern releases.
    # Deliberately lie about the format and ensure auto-detection still wins.
    rr = run([rust, "-f", "bam", "-s", "no", sam1, gtf], repo)
    hr = run([htseq, "-f", "bam", "-s", "no", sam1, gtf], repo)
    # HTSeq itself ignores -f and therefore succeeds on the SAM input.
    require_ok("Rust deprecated -f", rr)
    require_ok("HTSeq deprecated -f", hr)
    compare_tabular("deprecated -f", rr.stdout, hr.stdout, 1)
    print("[ OK ] deprecated format flag")

    quiet = run([rust, "-q", "-s", "no", sam1, gtf], repo)
    require_ok("Rust quiet", quiet)
    noisy_fragments = ("GFF lines processed", "records processed", "Creating IntervalTree")
    if any(fragment in quiet.stderr for fragment in noisy_fragments):
        raise AssertionError(f"--quiet still emitted progress:\n{quiet.stderr}")
    print("[ OK ] quiet")


def test_tabular(repo, rust, htseq, root, gtf, sam1, sam2):
    cases = [
        (
            "multi-file",
            ["-s","no","-n","2",str(sam1),str(sam2),str(gtf)],
            2, 0,
        ),
        (
            "with-header",
            ["-s","no","--with-header",str(sam1),str(sam2),str(gtf)],
            2, 0,
        ),
        (
            "metadata",
            [
                "-s","no",
                "--additional-attr","gene_name",
                "--add-chromosome-info",
                str(sam1),str(gtf),
            ],
            1, 2,
        ),
        (
            "feature-query",
            [
                "-s","no",
                "--feature-query",'gene_name == "Alpha"',
                str(sam1),str(gtf),
            ],
            1, 0,
        ),
        (
            "multi-id",
            [
                "-s","no",
                "-i","gene_id","-i","transcript_id",
                str(sam1),str(gtf),
            ],
            1, 0,
        ),
    ]

    for label, opts, n_samples, metadata_cols in cases:
        rr = run([rust, *opts], repo)
        hr = run([htseq, *opts], repo)
        require_ok(f"Rust {label}", rr)
        require_ok(f"HTSeq {label}", hr)
        rh, hh = compare_tabular(label, rr.stdout, hr.stdout, n_samples, metadata_cols)
        if label == "with-header" and rh != hh:
            raise AssertionError(f"{label}: headers differ: {rh!r} != {hh!r}")
        print(f"[ OK ] {label}")

    gz = root / "features.gtf.gz"
    with gzip.open(gz, "wt") as out:
        out.write(gtf.read_text())
    rr = run([rust, "-s","no",str(sam1),str(gz)], repo)
    hr = run([htseq, "-s","no",str(sam1),str(gz)], repo)
    require_ok("Rust gzip GTF", rr)
    require_ok("HTSeq gzip GTF", hr)
    compare_tabular("gzip GTF", rr.stdout, hr.stdout, 1)
    print("[ OK ] gzip GTF")

    sam_text = sam1.read_text()
    rr = run([rust,"-s","no","-",str(gtf)], repo, input_text=sam_text)
    hr = run([htseq,"-s","no","-",str(gtf)], repo, input_text=sam_text)
    require_ok("Rust stdin", rr)
    require_ok("HTSeq stdin", hr)
    compare_tabular("stdin", rr.stdout, hr.stdout, 1)
    print("[ OK ] stdin")

    odd = root / "alignment.data"
    shutil.copyfile(sam1, odd)
    rr = run([rust,"-s","no",str(odd),str(gtf)], repo)
    hr = run([htseq,"-s","no",str(odd),str(gtf)], repo)
    require_ok("Rust auto format", rr)
    require_ok("HTSeq auto format", hr)
    compare_tabular("auto format", rr.stdout, hr.stdout, 1)
    print("[ OK ] auto format")

    rust_append = root / "rust_append.tsv"
    ht_append = root / "ht_append.tsv"
    rust_append.write_text("PREEXISTING\n")
    ht_append.write_text("PREEXISTING\n")
    rr = run([
        rust, "-s", "no", "--append-output",
        "-c", rust_append, sam1, gtf
    ], repo)
    hr = run([
        htseq, "-s", "no", "--append-output",
        "-c", ht_append, sam1, gtf
    ], repo)
    require_ok("Rust append output", rr)
    require_ok("HTSeq append output", hr)
    if not rust_append.read_text().startswith("PREEXISTING\n"):
        raise AssertionError("Rust --append-output truncated existing content")
    if not ht_append.read_text().startswith("PREEXISTING\n"):
        raise AssertionError("HTSeq --append-output test fixture failed")
    compare_tabular(
        "append output",
        "\n".join(rust_append.read_text().splitlines()[1:]),
        "\n".join(ht_append.read_text().splitlines()[1:]),
        1,
    )
    print("[ OK ] append output")


def xf_rows(path, mode="r"):
    rows = []
    with pysam.AlignmentFile(str(path), mode) as bam:
        for rec in bam.fetch(until_eof=True):
            xf = rec.get_tag("XF") if rec.has_tag("XF") else None
            rows.append((rec.query_name, rec.flag, rec.reference_start, xf))
    return rows


def test_samout(repo, rust, htseq, root, gtf, sam1, sam2, paired):
    # Single-end and multiple samout destinations.
    rust_a = root / "rust_a.sam"
    rust_b = root / "rust_b.sam"
    ht_a = root / "ht_a.sam"
    ht_b = root / "ht_b.sam"

    rr = run([
        rust,"-s","no",
        "-o",rust_a,"-o",rust_b,
        sam1,sam2,gtf
    ], repo)
    hr = run([
        htseq,"-s","no",
        "-o",ht_a,"-o",ht_b,
        sam1,sam2,gtf
    ], repo)
    require_ok("Rust multiple samout", rr)
    require_ok("HTSeq multiple samout", hr)
    if xf_rows(rust_a) != xf_rows(ht_a) or xf_rows(rust_b) != xf_rows(ht_b):
        raise AssertionError("multiple --samout XF assignments differ")
    print("[ OK ] multiple samout")

    # Paired name-sorted SAM output.
    rust_pair = root / "rust_pair.sam"
    ht_pair = root / "ht_pair.sam"
    rr = run([rust,"-s","no","-r","name","-o",rust_pair,paired,gtf], repo)
    hr = run([htseq,"-s","no","-r","name","-o",ht_pair,paired,gtf], repo)
    require_ok("Rust paired samout", rr)
    require_ok("HTSeq paired samout", hr)
    if xf_rows(rust_pair) != xf_rows(ht_pair):
        raise AssertionError(
            f"paired SAM XF differs\nrust={xf_rows(rust_pair)}\nhtseq={xf_rows(ht_pair)}"
        )
    print("[ OK ] paired samout SAM")

    # Paired BAM samout format.
    rust_pair_bam = root / "rust_pair.bam"
    ht_pair_bam = root / "ht_pair.bam"
    rr = run([
        rust,"-s","no","-r","name","-p","BAM",
        "-o",rust_pair_bam,paired,gtf
    ], repo)
    hr = run([
        htseq,"-s","no","-r","name","-p","BAM",
        "-o",ht_pair_bam,paired,gtf
    ], repo)
    require_ok("Rust paired BAM samout", rr)
    require_ok("HTSeq paired BAM samout", hr)
    if xf_rows(rust_pair_bam, "rb") != xf_rows(ht_pair_bam, "rb"):
        raise AssertionError("paired BAM XF assignments differ")
    print("[ OK ] paired samout BAM")

    # Coordinate-sort the same pairs and repeat the paired SAM check.
    paired_bam = root / "paired.bam"
    with pysam.AlignmentFile(str(paired), "r") as src:
        with pysam.AlignmentFile(str(paired_bam), "wb", template=src) as dst:
            for rec in src.fetch(until_eof=True):
                dst.write(rec)
    pos_bam = root / "paired.pos.bam"
    pysam.sort("-o", str(pos_bam), str(paired_bam))

    rust_pos = root / "rust_pos.sam"
    ht_pos = root / "ht_pos.sam"
    rr = run([rust,"-s","no","-r","pos","-o",rust_pos,pos_bam,gtf], repo)
    hr = run([htseq,"-s","no","-r","pos","-o",ht_pos,pos_bam,gtf], repo)
    require_ok("Rust pos samout", rr)
    require_ok("HTSeq pos samout", hr)
    if xf_rows(rust_pos) != xf_rows(ht_pos):
        raise AssertionError("position-sorted paired XF assignments differ")
    print("[ OK ] position-sorted paired samout")


def test_cram(repo, rust, htseq, root, gtf, sam1):
    fasta = root / "ref.fa"
    fasta.write_text(">chr1\n" + "A"*1000 + "\n")
    pysam.faidx(str(fasta))

    cram = root / "reads.cram"
    with pysam.AlignmentFile(str(sam1), "r") as src:
        header = src.header.to_dict()
        header["SQ"][0]["UR"] = str(fasta.resolve())
        with pysam.AlignmentFile(
            str(cram),
            "wc",
            header=header,
            reference_filename=str(fasta),
        ) as dst:
            for rec in src.fetch(until_eof=True):
                dst.write(rec)

    rr = run([rust,"-s","no",cram,gtf], repo)
    hr = run([htseq,"-s","no",cram,gtf], repo)
    require_ok("Rust CRAM", rr)
    require_ok("HTSeq CRAM", hr)
    compare_tabular("CRAM", rr.stdout, hr.stdout, 1)
    print("[ OK ] CRAM")


def test_formats(repo, rust, htseq, root, gtf, sam1, sam2):
    import scipy.io
    import scipy.sparse
    import anndata
    import loompy

    # MTX sparse, direct comparison against HTSeq.
    rust_mtx = root / "rust.mtx"
    ht_mtx = root / "ht.mtx"
    common = [
        "-s","no","--counts-output-sparse",
        sam1,sam2,gtf,
    ]
    rr = run([rust,*common,"-c",rust_mtx], repo)
    hr = run([htseq,*common,"-c",ht_mtx], repo)
    require_ok("Rust MTX", rr)
    require_ok("HTSeq MTX", hr)
    rmat = scipy.io.mmread(rust_mtx).toarray() if scipy.sparse.issparse(scipy.io.mmread(rust_mtx)) else np.asarray(scipy.io.mmread(rust_mtx))
    htmp = scipy.io.mmread(ht_mtx)
    hmat = htmp.toarray() if scipy.sparse.issparse(htmp) else np.asarray(htmp)
    if rmat.shape != hmat.shape or not np.array_equal(rmat, hmat):
        raise AssertionError(f"MTX differs: rust={rmat} htseq={hmat}")
    print("[ OK ] MTX sparse")

    # H5AD sparse.
    rust_h5 = root / "rust.h5ad"
    ht_h5 = root / "ht.h5ad"
    rr = run([rust,*common,"-c",rust_h5], repo)
    hr = run([htseq,*common,"-c",ht_h5], repo)
    require_ok("Rust H5AD", rr)
    require_ok("HTSeq H5AD", hr)
    ra = anndata.read_h5ad(rust_h5)
    ha = anndata.read_h5ad(ht_h5)
    rx = ra.X.toarray() if scipy.sparse.issparse(ra.X) else np.asarray(ra.X)
    hx = ha.X.toarray() if scipy.sparse.issparse(ha.X) else np.asarray(ha.X)
    if not np.array_equal(rx, hx):
        raise AssertionError("H5AD X differs")
    if list(ra.obs_names) != list(ha.obs_names) or list(ra.var_names) != list(ha.var_names):
        raise AssertionError("H5AD axes differ")
    print("[ OK ] H5AD sparse")

    # Loom. HTSeq accepts the same count matrix; sparse flag is not needed.
    rust_loom = root / "rust.loom"
    ht_loom = root / "ht.loom"
    common_dense = ["-s","no",sam1,sam2,gtf]
    rr = run([rust,*common_dense,"-c",rust_loom], repo)
    hr = run([htseq,*common_dense,"-c",ht_loom], repo)
    require_ok("Rust Loom", rr)
    require_ok("HTSeq Loom", hr)
    with loompy.connect(str(rust_loom), mode="r") as rds, loompy.connect(str(ht_loom), mode="r") as hds:
        if rds.shape != hds.shape or not np.array_equal(rds[:, :], hds[:, :]):
            raise AssertionError("Loom matrix differs")
    print("[ OK ] Loom")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rust-bin", default="target/release/htseq_count_rust")
    ap.add_argument("--htseq-bin", default="htseq-count")
    ap.add_argument("--skip-formats", action="store_true")
    args = ap.parse_args()

    repo = Path.cwd()
    rust = str((repo / args.rust_bin).resolve())
    htseq = shutil.which(args.htseq_bin) or args.htseq_bin

    with tempfile.TemporaryDirectory(prefix="htseq-rust-compat-") as td:
        root = Path(td)
        gtf, sam1, sam2, paired = write_fixture(root)
        test_cli_compat(repo, rust, htseq, gtf, sam1)
        test_tabular(repo, rust, htseq, root, gtf, sam1, sam2)
        test_cram(repo, rust, htseq, root, gtf, sam1)
        test_samout(repo, rust, htseq, root, gtf, sam1, sam2, paired)
        if not args.skip_formats:
            test_formats(repo, rust, htseq, root, gtf, sam1, sam2)

    print("Compatibility surface: PASS")


if __name__ == "__main__":
    main()
