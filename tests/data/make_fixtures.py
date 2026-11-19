#!/usr/bin/env python3
"""Generate the VCF/gVCF test fixtures and the expected-call matrix.

This script is the source of truth for tests/data. Every expected value below
was derived by hand from the specification (see README.md "判讀邏輯"); the Rust
tests only read the generated files.

Requirements: python3, bgzip (htslib) and bcftools on PATH.
Usage: python3 tests/data/make_fixtures.py
"""

import gzip
import os
import subprocess
import sys

OUT = os.path.dirname(os.path.abspath(__file__))

# hg38 chr19:44908584-44908923 (1-based, inclusive), UCSC hg38. The same 340 bp
# window is identical in hg19 (offset +503257) and hs1 (offset +2824696).
WIN_START = 44908584
SEQ = (
    "GGAGGAACAACTGACCCCGGTGGCGGAGGAGACGCGGGCACGGCTGTCCAAGGAGCTGCAGGCGGCGCAGGCCCGGCTGG"
    "GCGCGGACATGGAGGACGTGTGCGGCCGCCTGGTGCAGTACCGCGGCGAGGTGCAGGCCATGCTCGGCCAGAGCACCGAG"
    "GAGCTGCGGGTGCGCCTCGCCTCCCACCTGCGCAAGCTGCGTAAGCGGCTCCTCCGCGATGCCGATGACCTGCAGAAGCG"
    "CCTGGCAGTGTACCAGGCCGGGGCCCGCGAGGGCGCCGAGCGCGGCCTCAGCGCCATCCGCGAGCGCCTGGGGCCCCTGG"
    "TGGAACAGGGCCGCGTGCGG"
)
assert len(SEQ) == 340

P1, P2 = 44908684, 44908822  # hg38 rs429358, rs7412

BUILDS = {
    "hg19": {"len": 59128983, "off": 503257},
    "hg38": {"len": 58617616, "off": 0},
    "hs1": {"len": 61707364, "off": 2824696},
}


def b(pos38):
    """Reference base at an hg38 position inside the window."""
    return SEQ[pos38 - WIN_START]


def seq(start38, end38):
    return SEQ[start38 - WIN_START:end38 - WIN_START + 1]


assert b(P1) == "T" and b(P2) == "C"
assert seq(44908680, 44908684) == "CGTGT"
assert seq(44908822, 44908823) == "CG"
assert seq(44908683, 44908685) == "GTG"
assert seq(44908821, 44908822) == "GC"

GT_FMT = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
PS_FMT = '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">'
DP_FMT = '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">'
AD_FMT = '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">'
GQ_FMT = '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">'
MIN_DP_FMT = '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">'
PL_FMT = '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">'
END_INFO = '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">'
PASS_FLT = '##FILTER=<ID=PASS,Description="All filters passed">'
LOWQUAL_FLT = '##FILTER=<ID=LowQual,Description="Low quality">'


def contig(name="chr19", build="hg38", length=True):
    if length:
        return f"##contig=<ID={name},length={BUILDS[build]['len']}>"
    return f"##contig=<ID={name}>"


def vcf(meta, samples, records, eol="\n"):
    lines = ["##fileformat=VCFv4.2"] + list(meta)
    cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if samples:
        cols += ["FORMAT"] + list(samples)
    lines.append("\t".join(cols))
    for r in records:
        lines.append("\t".join(str(x) for x in r))
    return eol.join(lines) + eol


def write(name, text):
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        fh.write(text)


EXPECTED = []  # (file, ref, sample, APOE, rs429358, rs7412, d1, d2, d3, d4, phased, note_contains)
ERRORS = []  # (file, ref, error_contains)

NA4 = ("NA", "NA", "NA", "NA")
DOSES = {
    "ε1/ε1": (2, 0, 0, 0), "ε1/ε2": (1, 1, 0, 0), "ε1/ε3": (1, 0, 1, 0), "ε1/ε4": (1, 0, 0, 1),
    "ε2/ε2": (0, 2, 0, 0), "ε2/ε3": (0, 1, 1, 0), "ε2/ε4": (0, 1, 0, 1),
    "ε3/ε3": (0, 0, 2, 0), "ε3/ε4": (0, 0, 1, 1), "ε4/ε4": (0, 0, 0, 2),
}
AMBIG = "ε2/ε4 or ε1/ε3"


def expect(fname, refs, sample, apoe, s1, s2, phased="no", notes=()):
    doses = NA4 if apoe == "NA" else DOSES[apoe]
    for ref in refs:
        EXPECTED.append((fname, ref, sample, apoe, s1, s2, *map(str, doses), phased, "|".join(notes)))


def gt_records(chrom, p1, p2, gts1, gts2, fmt="GT", id1="rs429358", id2="rs7412",
               ref1="T", alt1="C", ref2="C", alt2="T", flt1="PASS", flt2="PASS"):
    return [
        (chrom, p1, id1, ref1, alt1, 100, flt1, ".", fmt, *gts1),
        (chrom, p2, id2, ref2, alt2, 100, flt2, ".", fmt, *gts2),
    ]


# Six genotypes of table.md, in table order: (sample, rs429358 GT, rs7412 GT, APOE, display1, display2)
TABLE = [
    ("e2e2", "0/0", "1/1", "ε2/ε2", "T/T", "T/T"),
    ("e3e3", "0/0", "0/0", "ε3/ε3", "T/T", "C/C"),
    ("e2e3", "0/0", "0/1", "ε2/ε3", "T/T", "T/C"),
    ("e2e4", "0/1", "0/1", "ε2/ε4", "T/C", "T/C"),
    ("e3e4", "0/1", "0/0", "ε3/ε4", "T/C", "C/C"),
    ("e4e4", "1/1", "0/0", "ε4/ε4", "C/C", "C/C"),
]


def table_fixture(fname, build, chrom, contig_line):
    off = BUILDS[build]["off"]
    samples = [t[0] for t in TABLE]
    recs = gt_records(chrom, P1 + off, P2 + off, [t[1] for t in TABLE], [t[2] for t in TABLE])
    write(fname, vcf([PASS_FLT, GT_FMT, contig_line], samples, recs))
    for s, _, _, apoe, d1, d2 in TABLE:
        notes = (AMBIG,) if apoe == "ε2/ε4" else ()
        expect(fname, [build, "auto"], s, apoe, d1, d2, notes=notes)


# ---------------------------------------------------------------- VCF, unphased, hg38
F = "vcf_hg38_unphased.vcf"
rows = [
    # sample, rs429358 GT, rs7412 GT, APOE, display1, display2, notes
    ("e2e2", "0/0", "1/1", "ε2/ε2", "T/T", "T/T", ()),
    ("e2e3", "0/0", "0/1", "ε2/ε3", "T/T", "T/C", ()),
    ("e2e4", "0/1", "0/1", "ε2/ε4", "T/C", "T/C", (AMBIG,)),
    ("e3e3", "0/0", "0/0", "ε3/ε3", "T/T", "C/C", ()),
    ("e3e4", "0/1", "0/0", "ε3/ε4", "T/C", "C/C", ()),
    ("e4e4", "1/1", "0/0", "ε4/ε4", "C/C", "C/C", ()),
    ("e1e1", "1/1", "1/1", "ε1/ε1", "C/C", "T/T", ()),
    ("e1e2", "0/1", "1/1", "ε1/ε2", "T/C", "T/T", ()),
    ("e1e4", "1/1", "0/1", "ε1/ε4", "C/C", "T/C", ()),
    ("e3e4_rev", "1/0", "0/0", "ε3/ε4", "T/C", "C/C", ()),
    ("e2e3_rev", "0/0", "1/0", "ε2/ε3", "T/T", "T/C", ()),
    ("nocall_429", "./.", "0/0", "NA", "./.", "C/C", ("rs429358: no call (GT ./.)",)),
    ("nocall_7412", "0/1", "./.", "NA", "T/C", "./.", ("rs7412: no call (GT ./.)",)),
    ("partial", "0/.", "0/0", "NA", "T/.", "C/C", ("rs429358: no call (GT 0/.)",)),
    ("nocall_both", "./.", "./.", "NA", "./.", "./.", ("rs429358: no call", "rs7412: no call")),
]
write(F, vcf([PASS_FLT, GT_FMT, contig()], [r[0] for r in rows],
             gt_records("chr19", P1, P2, [r[1] for r in rows], [r[2] for r in rows])))
for s, _, _, apoe, d1, d2, notes in rows:
    expect(F, ["hg38", "auto"], s, apoe, d1, d2, notes=notes)

# ---------------------------------------------------------------- VCF, phased GT, hg38
F = "vcf_hg38_phased.vcf"
rows = [
    # sample, rs429358 GT:PS, rs7412 GT:PS, APOE, d1, d2, phased, notes
    ("p_e2e4", "0|1:.", "1|0:.", "ε2/ε4", "T/C", "T/C", "yes", ()),
    ("p_e1e3", "1|0:.", "1|0:.", "ε1/ε3", "T/C", "T/C", "yes", ()),
    ("p_e1e3_b", "0|1:.", "0|1:.", "ε1/ε3", "T/C", "T/C", "yes", ()),
    ("p_e2e4_b", "1|0:.", "0|1:.", "ε2/ε4", "T/C", "T/C", "yes", ()),
    ("p_e3e4", "0|1:.", "0|0:.", "ε3/ε4", "T/C", "C/C", "yes", ()),
    ("p_e2e2", "0|0:.", "1|1:.", "ε2/ε2", "T/T", "T/T", "yes", ()),
    ("p_ps_same", "1|0:100", "1|0:100", "ε1/ε3", "T/C", "T/C", "yes", ()),
    ("p_ps_diff", "1|0:100", "1|0:200", "ε2/ε4", "T/C", "T/C", "no", ("phase sets differ", AMBIG)),
    ("p_ps_one", "1|0:100", "1|0:.", "ε2/ε4", "T/C", "T/C", "no", ("phase sets differ", AMBIG)),
    ("p_mixed", "1|0:.", "1/0:.", "ε2/ε4", "T/C", "T/C", "no", (AMBIG,)),
    ("p_hom_diffps", "1|1:100", "0|1:200", "ε1/ε4", "C/C", "T/C", "no", ()),
    ("p_nocall", ".|.:.", "0|0:.", "NA", "./.", "C/C", "no", ("rs429358: no call (GT .|.)",)),
]
write(F, vcf([PASS_FLT, GT_FMT, PS_FMT, contig()], [r[0] for r in rows],
             gt_records("chr19", P1, P2, [r[1] for r in rows], [r[2] for r in rows], fmt="GT:PS")))
for s, _, _, apoe, d1, d2, ph, notes in rows:
    expect(F, ["hg38"], s, apoe, d1, d2, phased=ph, notes=notes)

# ---------------------------------------------------------------- builds: hg19 (b37 names), hg38, hs1
table_fixture("vcf_hg19_b37.vcf", "hg19", "19", contig("19", "hg19"))
table_fixture("vcf_hg38_table.vcf", "hg38", "chr19", contig("chr19", "hg38"))
table_fixture("vcf_hs1.vcf", "hs1", "chr19", contig("chr19", "hs1"))
for f, ref in [("vcf_hg19_b37.vcf", "hg38"), ("vcf_hg19_b37.vcf", "hs1"), ("vcf_hs1.vcf", "hg19"),
               ("vcf_hs1.vcf", "hg38"), ("vcf_hg38_table.vcf", "hg19"), ("vcf_hg38_table.vcf", "hs1")]:
    ERRORS.append((f, ref, "mismatch"))

# Same table as a plain (single-member) gzip file.
with open(os.path.join(OUT, "vcf_hg38_table.vcf"), "rb") as src:
    data = src.read()
with open(os.path.join(OUT, "vcf_hg38_table_gzip.vcf.gz"), "wb") as dst:
    dst.write(gzip.compress(data, mtime=0))
for s, _, _, apoe, d1, d2 in TABLE:
    expect("vcf_hg38_table_gzip.vcf.gz", ["hg38", "auto"], s, apoe, d1, d2,
           notes=(AMBIG,) if apoe == "ε2/ε4" else ())

# RefSeq accession contig name (hg38) with length.
F = "vcf_hg38_refseq.vcf"
write(F, vcf([PASS_FLT, GT_FMT, contig("NC_000019.10", "hg38")], ["acc_e3e4", "acc_e2e3"],
             gt_records("NC_000019.10", P1, P2, ["0/1", "0/0"], ["0/0", "0/1"])))
expect(F, ["hg38", "auto"], "acc_e3e4", "ε3/ε4", "T/C", "C/C")
expect(F, ["hg38", "auto"], "acc_e2e3", "ε2/ε3", "T/T", "T/C")
ERRORS.append((F, "hg19", "mismatch"))

# hg19 accession contig without length: build detected from the accession alone.
F = "vcf_hg19_accession_nolength.vcf"
off = BUILDS["hg19"]["off"]
write(F, vcf([PASS_FLT, GT_FMT, contig("NC_000019.9", "hg19", length=False)], ["acc_e4e4"],
             gt_records("NC_000019.9", P1 + off, P2 + off, ["1/1"], ["0/0"])))
expect(F, ["hg19", "auto"], "acc_e4e4", "ε4/ε4", "C/C", "C/C")
ERRORS.append((F, "hg38", "mismatch"))

# ---------------------------------------------------------------- gVCF
GATK_META = [
    '##ALT=<ID=NON_REF,Description="Represents any possible alternative allele not already represented at this location by REF and ALT">',
    PASS_FLT, AD_FMT, DP_FMT, GQ_FMT, GT_FMT, MIN_DP_FMT, PL_FMT,
    "##GVCFBlock0-20=minGQ=0(inclusive),maxGQ=20(exclusive)",
    "##GVCFBlock20-100=minGQ=20(inclusive),maxGQ=100(exclusive)",
    END_INFO, contig(),
]

F = "gvcf_gatk_hg38.g.vcf"
write(F, vcf(GATK_META, ["NA12878"], [
    ("chr19", 44908600, ".", b(44908600), "<NON_REF>", ".", ".", "END=44908683", "GT:DP:GQ:MIN_DP:PL", "0/0:35:99:30:0,90,1350"),
    ("chr19", P1, "rs429358", "T", "C,<NON_REF>", 812.64, ".", "DP=33;ExcessHet=3.0103", "GT:AD:DP:GQ:PL",
     "0/1:15,18,0:33:99:841,0,567,886,621,1507"),
    ("chr19", 44908685, ".", b(44908685), "<NON_REF>", ".", ".", "END=44908900", "GT:DP:GQ:MIN_DP:PL", "0/0:40:99:32:0,99,1485"),
]))
expect(F, ["hg38", "auto"], "NA12878", "ε3/ε4", "T/C", "C/C")

F = "gvcf_bcftools_hg38.g.vcf"
write(F, vcf([
    '##ALT=<ID=*,Description="Represents allele(s) other than observed.">',
    END_INFO, '##INFO=<ID=MinDP,Number=1,Type=Integer,Description="Minimum per-sample depth in this gVCF block">',
    GT_FMT, PL_FMT, contig(),
], ["S1"], [
    ("chr19", 44908650, ".", b(44908650), "<*>", 0, ".", "END=44908850;MinDP=20", "GT:PL", "0/0:0,60,255"),
]))
expect(F, ["hg38", "auto"], "S1", "ε3/ε3", "T/T", "C/C")

F = "gvcf_deepvariant_hg38.g.vcf"
write(F, vcf([
    '##ALT=<ID=*,Description="Represents any possible alternative allele at this location">',
    PASS_FLT, END_INFO, GT_FMT, GQ_FMT, DP_FMT, MIN_DP_FMT, AD_FMT,
    '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions">', PL_FMT, contig(),
], ["DV1"], [
    ("chr19", 44908600, ".", b(44908600), "<*>", 0, ".", "END=44908683", "GT:GQ:MIN_DP:PL", "0/0:50:28:0,84,839"),
    ("chr19", P1, ".", "T", "C,<*>", 45.2, "PASS", ".", "GT:GQ:DP:AD:VAF:PL", "1/1:44:30:0,30,0:1,0:45,46,0,990,990,990"),
    ("chr19", 44908685, ".", b(44908685), "<*>", 0, ".", "END=44908900", "GT:GQ:MIN_DP:PL", "0/0:50:31:0,90,899"),
]))
expect(F, ["hg38", "auto"], "DV1", "ε4/ε4", "C/C", "C/C")

F = "gvcf_gap_hg38.g.vcf"
write(F, vcf(GATK_META, ["G1"], [
    ("chr19", 44908600, ".", b(44908600), "<NON_REF>", ".", ".", "END=44908800", "GT:DP:GQ:MIN_DP:PL", "0/0:30:90:25:0,75,1125"),
    ("chr19", 44908850, ".", b(44908850), "<NON_REF>", ".", ".", "END=44909000", "GT:DP:GQ:MIN_DP:PL", "0/0:31:90:26:0,78,1170"),
]))
expect(F, ["hg38", "auto"], "G1", "NA", "T/T", "./.", notes=("rs7412: not covered by gVCF",))

F = "gvcf_nocall_block_hg38.g.vcf"
write(F, vcf(GATK_META, ["G2"], [
    ("chr19", 44908600, ".", b(44908600), "<NON_REF>", ".", ".", "END=44908800", "GT:DP:GQ:MIN_DP:PL", "0/0:30:90:25:0,75,1125"),
    ("chr19", 44908801, ".", b(44908801), "<NON_REF>", ".", ".", "END=44908900", "GT:DP:GQ:MIN_DP:PL", "./.:5:0:3:0,0,0"),
]))
expect(F, ["hg38", "auto"], "G2", "NA", "T/T", "./.", notes=("rs7412: no call (GT ./.)",))

F = "gvcf_zero_depth_hg38.g.vcf"
write(F, vcf(GATK_META, ["G3"], [
    ("chr19", 44908600, ".", b(44908600), "<NON_REF>", ".", ".", "END=44908700", "GT:DP:GQ:MIN_DP:PL", "0/0:12:0:0:0,0,0"),
    ("chr19", 44908701, ".", b(44908701), "<NON_REF>", ".", ".", "END=44908900", "GT:DP:GQ:PL", "0/0:0:0:0,0,0"),
]))
expect(F, ["hg38", "auto"], "G3", "NA", "./.", "./.",
       notes=("rs429358: zero depth in gVCF reference block", "rs7412: zero depth in gVCF reference block"))

F = "gvcf_strelka_hg38.g.vcf"
write(F, vcf([
    PASS_FLT, END_INFO,
    '##INFO=<ID=BLOCKAVG_min30p3a,Number=0,Type=Flag,Description="Non-variant multi-site block">',
    GT_FMT, '##FORMAT=<ID=GQX,Number=1,Type=Integer,Description="Empirically calibrated GQ">', DP_FMT,
    '##FORMAT=<ID=DPF,Number=1,Type=Integer,Description="Basecalls filtered from input prior to site genotyping">',
    MIN_DP_FMT, contig(),
], ["ST1"], [
    ("chr19", 44908600, ".", b(44908600), ".", ".", "PASS", "END=44908900;BLOCKAVG_min30p3a", "GT:GQX:DP:DPF:MIN_DP", "0/0:99:35:0:30"),
]))
expect(F, ["hg38", "auto"], "ST1", "ε3/ε3", "T/T", "C/C")

# ---------------------------------------------------------------- wildtype
F = "vcf_wildtype_explicit_hg38.vcf"
write(F, vcf([PASS_FLT, GT_FMT, contig()], ["W1", "W2"],
             gt_records("chr19", P1, P2, ["0/0", "0/0"], ["0/0", "0/0"])))
expect(F, ["hg38", "auto"], "W1", "ε3/ε3", "T/T", "C/C")
expect(F, ["hg38", "auto"], "W2", "ε3/ε3", "T/T", "C/C")

F = "vcf_wildtype_absent_hg38.vcf"
alt = {"A": "G", "C": "T", "G": "A", "T": "C"}
write(F, vcf([PASS_FLT, GT_FMT, contig()], ["A1", "A2"], [
    ("chr1", 1000000, ".", "A", "G", 50, "PASS", ".", "GT", "0/1", "0/0"),
    ("chr19", 44908600, ".", b(44908600), alt[b(44908600)], 50, "PASS", ".", "GT", "0/1", "1/1"),
    ("chr19", 44908900, ".", b(44908900), alt[b(44908900)], 50, "PASS", ".", "GT", "0/0", "0/1"),
]))
for s in ["A1", "A2"]:
    expect(F, ["hg38", "auto"], s, "ε3/ε3", "T/T", "C/C", notes=(
        "rs429358: not in VCF; assumed reference (T/T)", "rs7412: not in VCF; assumed reference (C/C)"))

F = "vcf_no_chr19_hg38.vcf"
write(F, vcf([PASS_FLT, GT_FMT, contig()], ["N1"], [
    ("chr1", 1000000, ".", "A", "G", 50, "PASS", ".", "GT", "0/1"),
    ("chr2", 2000000, ".", "C", "T", 50, "PASS", ".", "GT", "1/1"),
]))
expect(F, ["hg38", "auto"], "N1", "NA", "./.", "./.", notes=(
    "rs429358: VCF has no chr19 records", "rs7412: VCF has no chr19 records"))

# No ##contig header: the build cannot be detected.
F = "vcf_no_contig_hg38.vcf"
write(F, vcf([PASS_FLT, GT_FMT], ["U1"], [
    ("chr19", 44908600, ".", b(44908600), alt[b(44908600)], 50, "PASS", ".", "GT", "0/1"),
    ("chr19", P1, "rs429358", "T", "C", 50, "PASS", ".", "GT", "0/1"),
]))
expect(F, ["hg38"], "U1", "ε3/ε4", "T/C", "C/C", notes=("rs7412: not in VCF; assumed reference (C/C)",))
expect(F, ["auto"], "U1", "NA", "T/C", "./.", notes=(
    "rs7412: not in VCF and reference build not confirmed", "hg38 assumed"))

F = "vcf_no_contig_both_hg38.vcf"
write(F, vcf([PASS_FLT, GT_FMT], ["U2"], gt_records("chr19", P1, P2, ["0/1"], ["0/0"])))
expect(F, ["hg38"], "U2", "ε3/ε4", "T/C", "C/C")
expect(F, ["auto"], "U2", "ε3/ε4", "T/C", "C/C", notes=("hg38 assumed",))

# ---------------------------------------------------------------- input robustness
F = "vcf_crlf_hg38.vcf"
write(F, vcf([PASS_FLT, GT_FMT, contig()], ["C1", "C2"],
             gt_records("chr19", P1, P2, ["0/1", "0/0"], ["0/0", "1/1"]), eol="\r\n"))
expect(F, ["hg38", "auto"], "C1", "ε3/ε4", "T/C", "C/C")
expect(F, ["hg38", "auto"], "C2", "ε2/ε2", "T/T", "T/T")

F = "vcf_filter_hg38.vcf"
write(F, vcf([PASS_FLT, LOWQUAL_FLT, GT_FMT, contig()], ["F1"],
             gt_records("chr19", P1, P2, ["0/1"], ["0/0"], flt1="LowQual", flt2=".")))
expect(F, ["hg38"], "F1", "ε3/ε4", "T/C", "C/C", notes=("rs429358: FILTER=LowQual",))

F = "vcf_ploidy_hg38.vcf"
write(F, vcf([PASS_FLT, GT_FMT, contig()], ["H1", "H2"],
             gt_records("chr19", P1, P2, ["1", "0/0/1"], ["0/0", "0/0"])))
expect(F, ["hg38"], "H1", "NA", "C", "C/C", notes=("rs429358: non-diploid genotype (1)",))
expect(F, ["hg38"], "H2", "NA", "T/T/C", "C/C", notes=("rs429358: non-diploid genotype (0/0/1)",))

F = "vcf_ref_mismatch_hg38.vcf"
write(F, vcf([PASS_FLT, GT_FMT, contig()], ["R1"],
             gt_records("chr19", P1, P2, ["0/1"], ["0/0"], ref1="G", alt1="A")))
expect(F, ["hg38"], "R1", "NA", "./.", "C/C", notes=("rs429358: REF mismatch (expected T, found G)",))

F = "vcf_format_variants_hg38.vcf"
write(F, vcf([PASS_FLT, GT_FMT, DP_FMT, AD_FMT, contig()], ["V1", "V2"], [
    ("chr19", P1, "rs429358", "T", "C", 60, "PASS", ".", "DP:GT", "20:0/1", "18:0/0"),
    ("chr19", P2, "rs7412", "c", "t", 60, "PASS", ".", "GT:AD:DP", "0/0", "0/1:10,9:19"),
]))
expect(F, ["hg38"], "V1", "ε3/ε4", "T/C", "C/C")
expect(F, ["hg38"], "V2", "ε2/ε3", "T/T", "T/C")

F = "vcf_no_gt_hg38.vcf"
write(F, vcf([PASS_FLT, GT_FMT, DP_FMT, contig()], ["NG1"], [
    ("chr19", P1, "rs429358", "T", "C", 60, "PASS", ".", "DP", "25"),
    ("chr19", P2, "rs7412", "C", "T", 60, "PASS", ".", "GT", "0/0"),
]))
expect(F, ["hg38"], "NG1", "NA", "./.", "C/C", notes=("rs429358: no GT field",))

# ---------------------------------------------------------------- indels, MNPs, multi-allelic
DEL = "allele deleted by overlapping deletion"

F = "indel_hg38.vcf"
samples = ["i_del_het", "i_del_hom", "i_ins_het", "i_ins_snp", "i_del_after_7412", "i_del_and_snp",
           "i_snp_7412_del_after"]
gts = {  # R1 del(spans rs429358), R2 SNP rs429358, R3 ins after rs429358, R4 SNP rs7412, R5 del after rs7412
    "i_del_het": ["0/1", "0/0", "0/0", "0/0", "0/0"],
    "i_del_hom": ["1/1", "0/0", "0/0", "0/0", "0/0"],
    "i_ins_het": ["0/0", "0/0", "0/1", "0/0", "0/0"],
    "i_ins_snp": ["0/0", "0/1", "0/1", "0/0", "0/0"],
    "i_del_after_7412": ["0/0", "0/0", "0/0", "0/0", "0/1"],
    "i_del_and_snp": ["0/1", "0/1", "0/0", "0/0", "0/0"],
    "i_snp_7412_del_after": ["0/0", "0/0", "0/0", "0/1", "0/1"],
}
recs = [
    ("chr19", 44908680, ".", seq(44908680, 44908684), "C", 60, "PASS", ".", "GT"),
    ("chr19", P1, "rs429358", "T", "C", 60, "PASS", ".", "GT"),
    ("chr19", P1, ".", "T", "TGC", 60, "PASS", ".", "GT"),
    ("chr19", P2, "rs7412", "C", "T", 60, "PASS", ".", "GT"),
    ("chr19", P2, ".", seq(44908822, 44908823), "C", 60, "PASS", ".", "GT"),
]
write(F, vcf([PASS_FLT, GT_FMT, contig()], samples,
             [r + tuple(gts[s][i] for s in samples) for i, r in enumerate(recs)]))
expect(F, ["hg38"], "i_del_het", "NA", "T/del", "C/C", notes=("rs429358: " + DEL,))
expect(F, ["hg38"], "i_del_hom", "NA", "del/del", "C/C", notes=("rs429358: " + DEL,))
expect(F, ["hg38"], "i_ins_het", "ε3/ε3", "T/T", "C/C")
expect(F, ["hg38"], "i_ins_snp", "ε3/ε4", "T/C", "C/C")
expect(F, ["hg38"], "i_del_after_7412", "ε3/ε3", "T/T", "C/C")
expect(F, ["hg38"], "i_del_and_snp", "NA", "C/del", "C/C", notes=("rs429358: " + DEL,))
expect(F, ["hg38"], "i_snp_7412_del_after", "ε2/ε3", "T/T", "T/C")

F = "indel_star_hg38.vcf"
samples = ["s_star_het", "s_star_c", "s_snp_only", "s_ref"]
gts = {
    "s_star_het": ["0/1", "0/2", "0/0"],
    "s_star_c": ["0/1", "1/2", "0/0"],
    "s_snp_only": ["0/0", "0/1", "0/1"],
    "s_ref": ["0/0", "0/0", "0/0"],
}
recs = [
    ("chr19", 44908680, ".", seq(44908680, 44908684), "C", 60, "PASS", ".", "GT"),
    ("chr19", P1, "rs429358", "T", "C,*", 60, "PASS", ".", "GT"),
    ("chr19", P2, "rs7412", "C", "T", 60, "PASS", ".", "GT"),
]
write(F, vcf([PASS_FLT, GT_FMT, contig()], samples,
             [r + tuple(gts[s][i] for s in samples) for i, r in enumerate(recs)]))
expect(F, ["hg38"], "s_star_het", "NA", "T/del", "C/C", notes=("rs429358: " + DEL,))
expect(F, ["hg38"], "s_star_c", "NA", "C/del", "C/C", notes=("rs429358: " + DEL,))
expect(F, ["hg38"], "s_snp_only", "ε2/ε4", "T/C", "T/C", notes=(AMBIG,))
expect(F, ["hg38"], "s_ref", "ε3/ε3", "T/T", "C/C")

F = "indel_mnp_complex_hg38.vcf"
samples = ["m_mnp", "m_mnp7412", "m_complex"]
gts = {
    "m_mnp": ["0/1", "0/0", "0/0"],
    "m_mnp7412": ["0/0", "0/0", "1/1"],
    "m_complex": ["0/0", "0/1", "0/0"],
}
recs = [
    ("chr19", 44908683, ".", seq(44908683, 44908684), "GC", 60, "PASS", ".", "GT"),
    ("chr19", 44908683, ".", seq(44908683, 44908685), "A", 60, "PASS", ".", "GT"),
    ("chr19", 44908821, ".", seq(44908821, 44908822), "AT", 60, "PASS", ".", "GT"),
]
write(F, vcf([PASS_FLT, GT_FMT, contig()], samples,
             [r + tuple(gts[s][i] for s in samples) for i, r in enumerate(recs)]))
expect(F, ["hg38"], "m_mnp", "ε3/ε4", "T/C", "C/C")
expect(F, ["hg38"], "m_mnp7412", "ε2/ε2", "T/T", "T/T")
expect(F, ["hg38"], "m_complex", "NA", "T/?", "C/C", notes=("rs429358: unresolvable allele",))

F = "multiallelic_hg38.vcf"
write(F, vcf([PASS_FLT, GT_FMT, contig()], ["ma_1_2", "ma_0_1", "ma_0_2", "ma_2_2"],
             gt_records("chr19", P1, P2, ["1/2", "0/1", "0/2", "2/2"], ["0/0"] * 4, alt1="C,G")))
expect(F, ["hg38"], "ma_1_2", "NA", "C/G", "C/C", notes=("rs429358: unexpected allele G",))
expect(F, ["hg38"], "ma_0_1", "ε3/ε4", "T/C", "C/C")
expect(F, ["hg38"], "ma_0_2", "NA", "T/G", "C/C", notes=("rs429358: unexpected allele G",))
expect(F, ["hg38"], "ma_2_2", "NA", "G/G", "C/C", notes=("rs429358: unexpected allele G",))

F = "split_multiallelic_hg38.vcf"
samples = ["sp_c_het", "sp_c_g", "sp_conflict"]
gts = {
    "sp_c_het": ["0/1", "0/0", "0/0"],
    "sp_c_g": ["1/0", "0/1", "0/0"],
    "sp_conflict": ["1/1", "0/1", "0/0"],
}
recs = [
    ("chr19", P1, ".", "T", "C", 60, "PASS", ".", "GT"),
    ("chr19", P1, ".", "T", "G", 60, "PASS", ".", "GT"),
    ("chr19", P2, ".", "C", "T", 60, "PASS", ".", "GT"),
]
write(F, vcf([PASS_FLT, GT_FMT, contig()], samples,
             [r + tuple(gts[s][i] for s in samples) for i, r in enumerate(recs)]))
expect(F, ["hg38"], "sp_c_het", "ε3/ε4", "T/C", "C/C")
expect(F, ["hg38"], "sp_c_g", "NA", "C/G", "C/C", notes=("rs429358: unexpected allele G",))
expect(F, ["hg38"], "sp_conflict", "NA", "?/?", "C/C", notes=("rs429358: conflicting overlapping records",))

F = "sv_del_hg38.vcf"
write(F, vcf([
    '##ALT=<ID=DEL,Description="Deletion">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    END_INFO, PASS_FLT, GT_FMT, contig(),
], ["SV1"], [
    ("chr19", 44908600, "sv1", b(44908600), "<DEL>", 50, "PASS", "SVTYPE=DEL;END=44908900;SVLEN=-300", "GT", "0/1"),
]))
expect(F, ["hg38"], "SV1", "NA", "T/del", "C/del", notes=("rs429358: " + DEL, "rs7412: " + DEL))

# ---------------------------------------------------------------- bgzip, multiple BGZF blocks
F = "vcf_hg38_multiblock.vcf.gz"
plain = os.path.join(OUT, "vcf_hg38_multiblock.vcf")
lines = []
for c in range(1, 19):
    for i in range(250):
        lines.append(("chr%d" % c, 100000 + i * 997, ".", "A", "G", 50, "PASS", "DP=30", "GT", "0/1", "0/0", "1/1"))
for i in range(200):
    lines.append(("chr19", 1000000 + i * 1009, ".", "A", "G", 50, "PASS", "DP=30", "GT", "0/1", "0/0", "1/1"))
lines += gt_records("chr19", P1, P2, ["0/1", "1/1", "0/0"], ["0/0", "0/0", "0/1"])
for i in range(200):
    lines.append(("chr19", 50000000 + i * 1009, ".", "A", "G", 50, "PASS", "DP=30", "GT", "0/1", "0/0", "1/1"))
for i in range(300):
    lines.append(("chr20", 100000 + i * 997, ".", "A", "G", 50, "PASS", "DP=30", "GT", "0/1", "0/0", "1/1"))
write("vcf_hg38_multiblock.vcf", vcf([PASS_FLT, GT_FMT, contig()], ["B1", "B2", "B3"], lines))
assert os.path.getsize(plain) > 3 * 65536
subprocess.run(["bgzip", "-f", plain], check=True)
with open(os.path.join(OUT, F), "rb") as fh:
    blocks = fh.read().count(b"\x1f\x8b\x08\x04")
assert blocks >= 4, blocks
expect(F, ["hg38", "auto"], "B1", "ε3/ε4", "T/C", "C/C")
expect(F, ["hg38", "auto"], "B2", "ε4/ε4", "C/C", "C/C")
expect(F, ["hg38", "auto"], "B3", "ε2/ε3", "T/T", "T/C")

# ---------------------------------------------------------------- encodings and container quirks
def table_expect(fname, refs):
    for s, _, _, apoe, d1, d2 in TABLE:
        expect(fname, refs, s, apoe, d1, d2, notes=(AMBIG,) if apoe == "ε2/ε4" else ())


# UTF-8 BOM before ##fileformat (files saved by Windows editors).
F = "vcf_bom_hg38.vcf"
write(F, "﻿" + open(os.path.join(OUT, "vcf_hg38_table.vcf"), encoding="utf-8").read())
table_expect(F, ["hg38", "auto"])

# gzip followed by zero padding (ignored silently) or by non-gzip garbage
# (decoded data kept; a warning is emitted).
with open(os.path.join(OUT, "vcf_hg38_table_gzip.vcf.gz"), "rb") as fh:
    gz = fh.read()
with open(os.path.join(OUT, "vcf_hg38_gzip_zero_padding.vcf.gz"), "wb") as fh:
    fh.write(gz + b"\0" * 512)
table_expect("vcf_hg38_gzip_zero_padding.vcf.gz", ["hg38", "auto"])
with open(os.path.join(OUT, "vcf_hg38_gzip_trailing_garbage.vcf.gz"), "wb") as fh:
    fh.write(gz + b"trailing garbage\n")
table_expect("vcf_hg38_gzip_trailing_garbage.vcf.gz", ["hg38", "auto"])
# A gzip member cut short is still a fatal I/O error.
with open(os.path.join(OUT, "err_truncated_gzip.vcf.gz"), "wb") as fh:
    fh.write(gz[: len(gz) // 2])

# VCF 4.4 explicit phasing prefix on the first allele ("|0|1", "/0/1", "|1").
F = "vcf_hg38_gt_prefix_v44.vcf"
text = vcf([PASS_FLT, GT_FMT, contig()], ["g44_phased", "g44_unphased", "g44_haploid"],
           gt_records("chr19", P1, P2, ["|0|1", "/0/1", "|1"], ["|1|0", "/0/1", "0/0"]))
write(F, text.replace("##fileformat=VCFv4.2", "##fileformat=VCFv4.4", 1))
expect(F, ["hg38"], "g44_phased", "ε2/ε4", "T/C", "T/C", phased="yes")
expect(F, ["hg38"], "g44_unphased", "ε2/ε4", "T/C", "T/C", notes=(AMBIG,))
expect(F, ["hg38"], "g44_haploid", "NA", "C", "C/C", notes=("rs429358: non-diploid genotype",))

# ---------------------------------------------------------------- fatal errors
write("err_sites_only.vcf", vcf([PASS_FLT, contig()], [], [
    ("chr19", P1, "rs429358", "T", "C", 50, "PASS", "."),
]))
ERRORS.append(("err_sites_only.vcf", "hg38", "no sample columns"))
write("err_not_vcf.txt", "hello\nthis is not a VCF file\n")
ERRORS.append(("err_not_vcf.txt", "hg38", "#CHROM"))
write("err_empty.vcf", "")
ERRORS.append(("err_empty.vcf", "hg38", "#CHROM"))
write("err_no_chrom_header.vcf", "##fileformat=VCFv4.2\n" + "\t".join(
    ["chr19", str(P1), "rs429358", "T", "C", "50", "PASS", ".", "GT", "0/1"]) + "\n")
ERRORS.append(("err_no_chrom_header.vcf", "hg38", "#CHROM"))
write("err_space_delimited.vcf", "##fileformat=VCFv4.2\n#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1\n"
      "chr19 44908684 rs429358 T C 50 PASS . GT 0/1\n")
ERRORS.append(("err_space_delimited.vcf", "hg38", "tab"))
write("err_truncated_record.vcf", vcf([PASS_FLT, GT_FMT, contig()], ["S1"], [
    ("chr19", P1, "rs429358", "T", "C"),
]))
ERRORS.append(("err_truncated_record.vcf", "hg38", "line|columns"))
write("err_sample_count.vcf", vcf([PASS_FLT, GT_FMT, contig()], ["S1", "S2"], [
    ("chr19", P1, "rs429358", "T", "C", 50, "PASS", ".", "GT", "0/1"),
]))
ERRORS.append(("err_sample_count.vcf", "hg38", "line|columns"))
subprocess.run(["bcftools", "view", "--no-version", "-Ob", "-o", os.path.join(OUT, "err_bcf.bcf"),
                os.path.join(OUT, "vcf_hg38_table.vcf")], check=True)
ERRORS.append(("err_bcf.bcf", "hg38", "BCF"))
ERRORS.append(("err_truncated_gzip.vcf.gz", "hg38", "I/O error"))
ERRORS.append(("does_not_exist.vcf", "hg38", "cannot open"))

# ---------------------------------------------------------------- matrix files
with open(os.path.join(OUT, "expected_calls.tsv"), "w", newline="\n") as fh:
    fh.write("#file\tref\tsample\tAPOE\trs429358\trs7412\te1_dose\te2_dose\te3_dose\te4_dose\tphased\tnote_contains\n")
    for row in EXPECTED:
        fh.write("\t".join(row) + "\n")
with open(os.path.join(OUT, "expected_errors.tsv"), "w", newline="\n") as fh:
    fh.write("#file\tref\terror_contains\n")
    for row in ERRORS:
        fh.write("\t".join(row) + "\n")

COLS = ["sample", "APOE", "rs429358 (chr19:{p1})", "rs7412 (chr19:{p2})",
        "e1_dose", "e2_dose", "e3_dose", "e4_dose", "phased", "note"]
with open(os.path.join(OUT, "expected_headers.tsv"), "w", newline="\n") as fh:
    fh.write("#build\tcolumns...\n")
    for build in ["hg19", "hg38", "hs1"]:
        off = BUILDS[build]["off"]
        cols = [c.format(p1=P1 + off, p2=P2 + off) for c in COLS]
        fh.write(build + "\t" + "\t".join(cols) + "\n")

print(f"{len(EXPECTED)} expected calls, {len(ERRORS)} expected errors", file=sys.stderr)
