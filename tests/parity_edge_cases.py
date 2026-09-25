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
