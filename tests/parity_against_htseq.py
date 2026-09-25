#!/usr/bin/env python3
import argparse
import math
import re
import random
import shutil
import subprocess
import tempfile
from pathlib import Path

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

def qlen(cigar):
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in "MIS=X")

def sam(name, pos=100, cigar="10M", flag=0, mapq=60, tags=()):
    n = qlen(cigar)
    seq = "A" * n if n else "*"
    qual = "I" * n if n else "*"
    fields = [name, str(flag), "chr1", str(pos), str(mapq), cigar, "*", "0", "0", seq, qual]
    return "\t".join(fields + list(tags))

def paired_sam(name, flag, pos, mate_pos, cigar="10M", mapq=60, tags=(), tlen=0):
    n = qlen(cigar)
    seq = "A" * n if n else "*"
    qual = "I" * n if n else "*"
    fields = [
        name, str(flag), "chr1", str(pos), str(mapq), cigar, "=",
        str(mate_pos), str(tlen), seq, qual,
    ]
    return "\t".join(fields + list(tags))

def paired_unmapped_sam(name, flag, mate_pos, mapq=255, tags=()):
    fields = [
        name, str(flag), "*", "0", str(mapq), "*", "chr1",
        str(mate_pos), "0", "*", "*",
    ]
    return "\t".join(fields + list(tags))

def exon(start, end, gene, strand="+", gene_name=None, extra="", chrom="chr1", feature_type="exon"):
    gene_name = gene_name or f"{gene}_name"
    attrs = f'gene_id "{gene}"; gene_name "{gene_name}";'
    if extra:
        attrs += " " + extra
    return f"{chrom}\tparity\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"

def opts(mode="union", nonunique="none", stranded="no", idattrs=("gene_id",),
         secondary="ignore", supplementary="ignore", minaqual=10,
         feature_types=("exon",), order=None):
    out = ["-s", stranded]
    for feature_type in feature_types:
        out += ["-t", feature_type]
    out += ["-m", mode, "-a", str(minaqual), "--nonunique", nonunique]
    for attr in idattrs:
        out += ["-i", attr]
    if secondary is not None:
        out += ["--secondary-alignments", secondary]
    if supplementary is not None:
        out += ["--supplementary-alignments", supplementary]
    if order is not None:
        out += ["-r", order]
    return out

CASES = [
    ("basic_union",
     [exon(100,109,"geneA")],
     [sam("r1")],
     opts()),

    ("spliced_union",
     [exon(100,104,"geneA"), exon(115,119,"geneA")],
     [sam("r1", cigar="5M10N5M")],
     opts()),

    ("multiple_feature_types",
     [exon(100,109,"geneA"), exon(120,129,"geneB", feature_type="pseudogene")],
     [sam("r1", pos=100), sam("r2", pos=120)],
     opts(feature_types=("exon", "pseudogene"))),

    ("eq_x_cigar",
     [exon(100,109,"geneA")],
     [sam("r1", cigar="5=1X4=")],
     opts()),

    ("leading_insertion",
     [exon(100,107,"geneA"), exon(108,108,"geneB")],
     [sam("r1", cigar="2I8M")],
     opts()),

    ("leading_deletion",
     [exon(102,109,"geneA"), exon(110,110,"geneB")],
     [sam("r1", cigar="2D8M")],
     opts()),

    ("internal_insertion",
     [exon(100,107,"geneA"), exon(108,109,"geneB")],
     [sam("r1", cigar="4M2I4M")],
     opts()),

    ("internal_deletion",
     [exon(100,103,"geneA"), exon(106,109,"geneA")],
     [sam("r1", cigar="4M2D4M")],
     opts()),

    ("soft_clipped",
     [exon(100,107,"geneA")],
     [sam("r1", cigar="3S8M2S")],
     opts()),

    ("intersection_nonempty_changing_sets",
     [exon(100,106,"geneA"), exon(104,110,"geneB")],
     [sam("r1", cigar="11M")],
     opts(mode="intersection-nonempty")),

    ("intersection_nonempty_ignores_empty_gap",
     [exon(100,102,"geneA"), exon(108,110,"geneA")],
     [sam("r1", cigar="11M")],
     opts(mode="intersection-nonempty")),

    ("intersection_strict_empty_gap",
     [exon(100,102,"geneA"), exon(108,110,"geneA")],
     [sam("r1", cigar="11M")],
     opts(mode="intersection-strict")),

    ("nh_nonunique_all",
     [exon(100,109,"geneA")],
     [sam("r1", tags=("NH:i:2",))],
     opts(nonunique="all")),

    ("nh_nonunique_fraction",
     [exon(100,109,"geneA")],
     [sam("r1", tags=("NH:i:2",))],
     opts(nonunique="fraction")),

    ("fraction_three_features",
     [exon(100,109,"geneA"), exon(100,109,"geneB"), exon(100,109,"geneC")],
     [sam("r1")],
     opts(nonunique="fraction")),

    ("multiple_idattr",
     [exon(100,109,"geneA", extra='exon_number "7";')],
     [sam("r1")],
     opts(idattrs=("gene_id","exon_number"))),

    ("low_mapq_plus_nh",
     [exon(100,109,"geneA")],
     [sam("r1", mapq=0, tags=("NH:i:2",))],
     opts(nonunique="none", minaqual=10)),

    ("unmapped_flag_with_reference_fields",
     [exon(100,109,"geneA")],
     [sam("r1", flag=4, pos=100, cigar="10M")],
     opts()),

    ("reverse_stranded",
     [exon(100,109,"geneA", strand="-")],
     [sam("r1", flag=16, pos=100, cigar="10M")],
     opts(stranded="yes")),

    ("same_id_opposite_strands",
     [exon(100,109,"geneA", strand="+"), exon(100,109,"geneA", strand="-")],
     [sam("r1", flag=0, pos=100, cigar="10M")],
     opts(stranded="yes")),

    ("paired_end_same_gene",
     [exon(100,300,"geneA")],
     [paired_sam("pair1", 99, 120, 220), paired_sam("pair1", 147, 220, 120)],
     opts(order="name")),

    ("paired_end_stranded_yes",
     [exon(100,300,"geneA", strand="+")],
     [paired_sam("pair1", 99, 120, 220), paired_sam("pair1", 147, 220, 120)],
     opts(stranded="yes", order="name")),

    ("paired_end_stranded_reverse",
     [exon(100,300,"geneA", strand="-")],
     [paired_sam("pair1", 99, 120, 220), paired_sam("pair1", 147, 220, 120)],
     opts(stranded="reverse", order="name")),

    ("paired_end_ambiguous_union",
     [exon(110,150,"geneA"), exon(210,250,"geneB")],
     [paired_sam("pair1", 99, 120, 220), paired_sam("pair1", 147, 220, 120)],
     opts(order="name")),

    ("paired_end_missing_mate",
     [exon(100,180,"geneA")],
     [paired_sam("pair1", 99, 120, 220)],
     opts(order="name")),

    ("paired_end_one_unmapped",
     [exon(100,180,"geneA")],
     [
         paired_sam("pair1", 73, 120, 0),
         paired_unmapped_sam("pair1", 133, 120),
     ],
     opts(order="name")),

    ("paired_end_nh_none",
     [exon(100,300,"geneA")],
     [
         paired_sam("pair1", 99, 120, 220, tags=("NH:i:2",)),
         paired_sam("pair1", 147, 220, 120),
     ],
     opts(nonunique="none", order="name")),

    ("paired_end_nh_all",
     [exon(100,300,"geneA")],
     [
         paired_sam("pair1", 99, 120, 220, tags=("NH:i:2",)),
         paired_sam("pair1", 147, 220, 120),
     ],
     opts(nonunique="all", order="name")),

    ("paired_end_position_sorted_interleaved",
     [exon(90,330,"geneA")],
     [
         paired_sam("pairA", 99, 100, 300),
         paired_sam("pairB", 99, 150, 250),
         paired_sam("pairB", 147, 250, 150),
         paired_sam("pairA", 147, 300, 100),
     ],
     opts(order="pos")),

    ("paired_end_position_sorted_missing_mate",
     [exon(90,180,"geneA")],
     [paired_sam("pairA", 99, 100, 300)],
     opts(order="pos")),

    ("paired_end_spliced_cigar",
     [exon(100,104,"geneA"), exon(115,119,"geneA"), exon(200,209,"geneA")],
     [
         paired_sam("pair1", 99, 100, 200, cigar="5M10N5M"),
         paired_sam("pair1", 147, 200, 100),
     ],
     opts(order="name")),

    ("paired_end_low_mapq",
     [exon(100,300,"geneA")],
     [
         paired_sam("pair1", 99, 120, 220, mapq=5),
         paired_sam("pair1", 147, 220, 120),
     ],
     opts(order="name")),

    ("paired_end_secondary_ignore",
     [exon(100,300,"geneA")],
     [
         paired_sam("pair1", 99, 120, 220),
         paired_sam("pair1", 147, 220, 120),
         paired_sam("pair1", 355, 120, 220),
         paired_sam("pair1", 403, 220, 120),
     ],
     opts(secondary="ignore", supplementary="ignore", order="name")),

    ("paired_end_secondary_score",
     [exon(100,300,"geneA")],
     [
         paired_sam("pair1", 99, 120, 220),
         paired_sam("pair1", 147, 220, 120),
         paired_sam("pair1", 355, 120, 220),
         paired_sam("pair1", 403, 220, 120),
     ],
     opts(secondary="score", supplementary="ignore", order="name")),

    ("paired_end_supplementary_score",
     [exon(100,300,"geneA")],
     [
         paired_sam("pair1", 99, 120, 220),
         paired_sam("pair1", 147, 220, 120),
         paired_sam("pair1", 2147, 120, 220),
         paired_sam("pair1", 2195, 220, 120),
     ],
     opts(secondary="ignore", supplementary="score", order="name")),

    ("paired_end_intersection_nonempty",
     [exon(100,250,"geneA"), exon(210,230,"geneB")],
     [paired_sam("pair1", 99, 120, 220), paired_sam("pair1", 147, 220, 120)],
     opts(mode="intersection-nonempty", order="name")),

    ("paired_end_intersection_strict",
     [exon(100,250,"geneA"), exon(215,225,"geneB")],
     [paired_sam("pair1", 99, 120, 220), paired_sam("pair1", 147, 220, 120)],
     opts(mode="intersection-strict", order="name")),

    ("default_idattr",
     [exon(100,109,"ENSG_TEST", gene_name="GENE_TEST")],
     [sam("r1")],
     ["-s","no","-t","exon","-m","union","-a","10","--nonunique","none",
      "--secondary-alignments","ignore","--supplementary-alignments","ignore"]),

    ("default_secondary",
     [exon(100,109,"geneA")],
     [sam("r1", flag=256)],
     opts(secondary=None, supplementary=None)),

    ("default_supplementary",
     [exon(100,109,"geneA")],
     [sam("r1", flag=2048)],
     opts(secondary=None, supplementary=None)),

    ("annotation_contig_not_in_sam_header",
     [exon(100,109,"geneA"), exon(100,109,"geneB", chrom="chr2")],
     [sam("r1")],
     opts()),
]

def run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

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
    return [(k, a.get(k), b.get(k)) for k in keys
            if a.get(k) is None or b.get(k) is None
            or not math.isclose(a[k], b[k], rel_tol=1e-9, abs_tol=1e-9)]

def random_cigar(rng):
    # All are valid single-end CIGARs. They deliberately mix operations that
    # consume reference, query, both, or neither.
    choices = [
        "1M", "5M", "10M", "3M2N4M", "4M1D5M", "2I8M",
        "3S7M", "5=1X4=", "2M1I3M2D4M", "1S3M4N2M1I2M",
    ]
    return rng.choice(choices)

def fuzz_cases(seed, count):
    rng = random.Random(seed)
    generated = []
    modes = ["union", "intersection-strict", "intersection-nonempty"]
    stranded_values = ["no", "yes", "reverse"]
    nonunique_values = ["none", "all", "fraction"]

    for idx in range(count):
        gtfs = []
        gene_count = rng.randint(1, 4)
        for gi in range(gene_count):
            gene = f"gene{gi}"
            strand = rng.choice(["+", "-"])
            exon_count = rng.randint(1, 3)
            for _ in range(exon_count):
                start = rng.randint(85, 135)
                length = rng.randint(1, 20)
                gtfs.append(exon(start, start + length - 1, gene, strand=strand))

        cigar = random_cigar(rng)
        pos = rng.randint(90, 130)
        flag = 16 if rng.choice([False, True]) else 0
        mapq = rng.choice([0, 5, 10, 20, 60])
        tags = ()
        if rng.random() < 0.25:
            tags = ("NH:i:2",)

        mode = rng.choice(modes)
        stranded = rng.choice(stranded_values)
        nonunique = rng.choice(nonunique_values)
        generated.append((
            f"fuzz_{idx:04d}",
            gtfs,
            [sam("r1", pos=pos, cigar=cigar, flag=flag, mapq=mapq, tags=tags)],
            opts(
                mode=mode,
                nonunique=nonunique,
                stranded=stranded,
                minaqual=10,
            ),
        ))
    return generated

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rust-bin", default="target/release/htseq_count_rust")
    ap.add_argument("--htseq-bin", default="htseq-count")
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--case", action="append", default=[])
    ap.add_argument("--fuzz", type=int, default=0,
                    help="Add N deterministic randomized single-end parity cases.")
    ap.add_argument("--fuzz-seed", type=int, default=1337)
    args = ap.parse_args()

    root = Path.cwd()
    rust = root / args.rust_bin
    if args.build or not rust.exists():
        subprocess.run(["cargo", "build", "--release"], cwd=root, check=True)
    rust = str(rust.resolve())
    htseq = shutil.which(args.htseq_bin) or args.htseq_bin

    hv = run([htseq, "--version"], root)
    print("HTSeq:", (hv.stdout + hv.stderr).strip())
    print("Rust binary:", rust)
    print()

    cases = list(CASES)
    if args.fuzz:
        cases.extend(fuzz_cases(args.fuzz_seed, args.fuzz))
    if args.case:
        wanted = set(args.case)
        cases = [case for case in cases if case[0] in wanted]

    different = 0
    with tempfile.TemporaryDirectory(prefix="htseq-rust-parity-") as td:
        root_tmp = Path(td)
        for name, gtfs, sams, op in cases:
            case = root_tmp / name
            case.mkdir()
            gtf = case / "features.gtf"
            sf = case / "reads.sam"
            gtf.write_text("\n".join(gtfs) + "\n")
            sf.write_text("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n"
                          + "\n".join(sams) + "\n")

            rr = run([rust] + op + [str(sf), str(gtf)], root)
            hr = run([htseq] + op + [str(sf), str(gtf)], root)

            if rr.returncode or hr.returncode:
                different += 1
                print(f"[DIFF] {name}: exit rust={rr.returncode}, htseq={hr.returncode}")
                if name.startswith("fuzz_"):
                    print("  opts:", " ".join(op))
                    print("  GTF:")
                    for line in gtfs:
                        print("   ", line)
                    print("  SAM:")
                    for line in sams:
                        print("   ", line)
                if rr.stderr.strip():
                    print("  rust stderr:", rr.stderr.strip().replace("\n", " | "))
                if hr.stderr.strip():
                    print("  htseq stderr:", hr.stderr.strip().replace("\n", " | "))
                continue

            d = delta(parse_counts(rr.stdout), parse_counts(hr.stdout))
            if d:
                different += 1
                print(f"[DIFF] {name}")
                if name.startswith("fuzz_"):
                    print("  opts:", " ".join(op))
                    print("  GTF:")
                    for line in gtfs:
                        print("   ", line)
                    print("  SAM:")
                    for line in sams:
                        print("   ", line)
                for row in d:
                    print(" ", row)
            else:
                print(f"[ OK ] {name}")

    print()
    print(f"Differing cases: {different}/{len(cases)}")
    raise SystemExit(1 if different else 0)

if __name__ == "__main__":
    main()
