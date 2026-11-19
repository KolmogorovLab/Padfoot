"""Tiers, gene-SV pairs, whole-gene labels and the gene table (by_gene)."""

import os
from collections import defaultdict

import pytest

from padfoot import annot
from toy_genome import exon


def DEL(sid, start, end, chrom="chr1", svtype="DEL", **kw):
    return dict(id=sid, chrom=chrom, pos=start, end=end, svlen=end - start, svtype=svtype, **kw)


def BND(sid, pos, mate, strands="+-", chrom="chr1", **kw):
    return dict(id=sid, chrom=chrom, pos=pos, mate=mate, svtype="BND", strands=strands, **kw)


def pairs(sv):
    return {p.gene.name: p for p in sv.pairs}


E2 = exon("GENEA", 2)
D1 = exon("GENED", 1)       # entirely 5' UTR


@pytest.mark.parametrize("record,gene,cls,tier", [
    (DEL("S", E2[0] + 10, E2[0] + 20), "GENEA", "frameshift", "HIGH"),
    (DEL("S", D1[0] + 10, D1[0] + 22), "GENED", "exon_altering", "MODERATE"),     # UTR-only exon
    (DEL("S", 103_000, 106_000), "GENEA", "between_exons", "HIGH"),
    (DEL("S", 105_500, 106_000), "GENEA", "intronic", "MODERATE"),
    (DEL("S", 132_000, 133_000), "GENEA", "promoter/UTR", "MODERATE"),
    (DEL("S", 105_500, 160_000), "GENEA", "Gene-NonCoding", "HIGH"),               # breakend in gene body
    (DEL("S", 105_500, 160_000, svtype="DUP", strands="-+"), "GENEA", "Gene-NonCoding", "MODERATE"),  # partial DUP
    (DEL("S", 132_000, 160_000), "GENEA", "Gene-NonCoding", "MODERATE"),           # padding only
    (DEL("S", 105_500, 205_000), "GENEB", "between_genes", "HIGH"),
    (BND("S", 105_500, ("chr2", 105_000)), "GENEG", "possible_fusion", "HIGH"),
    (BND("S", 105_500, ("chr2", 95_000)), "GENEG", "between_genes", "MODERATE"),  # padding only, no fusion
    (dict(id="S", chrom="chr1", pos=E2[0] + 50, svtype="sBND"), "GENEA", "exon_altering", "HIGH"),
    (dict(id="S", chrom="chr1", pos=106_000, svlen=40, svtype="INS", alt="A" * 40), "GENEA", "intronic", "MODERATE"),
])
def test_pair_class_and_tier(classify, record, gene, cls, tier):
    sv = next(iter(classify([record]).values()))
    p = pairs(sv)[gene]
    assert (p.cls, p.tier) == (cls, tier)


def test_tier_does_not_depend_on_cancer_genes(classify):
    sv = classify([DEL("S", 105_500, 106_000)])["S"]
    annot.cancer_annot_svs([sv], cancer_genes={"GENEA": "TSG"}, fusion=defaultdict(list))
    assert sv.impact == "MODERATE"
    assert sv.cancer[0] == "TSG"


def test_sv_tier_is_most_severe_pair(classify):
    # bp1 in GENEA intron 3 (HIGH for GENEA), bp2 in GENEB's left padding only (MODERATE)
    sv = classify([DEL("S", 105_500, 195_000)])["S"]
    assert pairs(sv)["GENEB"].tier == "MODERATE"
    assert sv.impact == "HIGH"


def test_every_overlapping_gene_gets_a_pair(classify):
    # 316,000 and 318,000 are in GENED's body and in GENEC's right padding
    sv = classify([DEL("S", 316_000, 318_000)])["S"]
    assert sv.genes[2] == "intronic"
    assert (pairs(sv)["GENED"].cls, pairs(sv)["GENEC"].cls) == ("intronic", "promoter/UTR")


# --- fusions -----------------------------------------------------------------------

@pytest.mark.parametrize("record,partners", [
    (BND("F", 105_500, ("chr2", 105_000), strands="+-"), "GENEA::GENEG"),
    (BND("F", 105_000, ("chr1", 105_500), strands="+-", chrom="chr2"), "GENEG::GENEA"),
])
def test_fusion_partners_are_5p_3p(classify, record, partners):
    sv = next(iter(classify([record]).values()))
    assert sv.genes[2] == "possible_fusion"
    assert sv.fusion == partners


def test_no_fusion_with_incompatible_orientation(classify):
    sv = classify([BND("F", 105_500, ("chr2", 105_000), strands="++")])["F_1"]
    assert sv.genes[2] == "between_genes"
    assert sv.fusion == ""


# --- whole genes ---------------------------------------------------------------------

def test_whole_gene_label_without_tier(classify):
    # GENEC and GENED lie entirely inside; bp1 is in GENEC's left padding
    sv = classify([DEL("W", 290_000, 336_000)])["W"]
    ps = pairs(sv)
    assert set(ps) == {"GENEC", "GENED"}
    assert all(p.labels == ["whole_gene_deletion"] and p.tier == "" for p in ps.values())
    assert sv.impact == ""


def test_whole_gene_duplication_and_inversion(classify):
    dup = classify([DEL("W", 290_000, 336_000, svtype="DUP", strands="-+")])["W"]
    inv = classify([DEL("W", 290_000, 336_000, svtype="INV", strands="++")])["W"]
    assert pairs(dup)["GENED"].labels == ["whole_gene_duplication"]
    assert pairs(inv)["GENED"].labels == ["whole_gene_inversion"]


def test_whole_gene_size_cap(toy_annotation, parse_svs):
    svs = parse_svs([DEL("W", 290_000, 336_000)])
    annot.annot_SVS(toy_annotation, svs, defaultdict(list), whole_gene_max=10_000)
    ps = pairs(svs[0])
    assert set(ps) == {"GENEC"}                 # only the padding hit remains
    assert (ps["GENEC"].labels, ps["GENEC"].tier) == ([], "MODERATE")


def test_alt_isoform_exon_label(classify):
    sv = classify([DEL("S", 106_550, 106_600)])["S"]
    assert pairs(sv)["GENEA"].labels == ["alt_isoform_exon"]


# --- gene table ----------------------------------------------------------------------

def run(toy_annotation, parse_svs, records):
    svs = parse_svs(records)
    by_gene = defaultdict(list)
    annot.annot_SVS(toy_annotation, svs, by_gene)
    return svs, by_gene


def test_gene_impact_is_most_severe_per_haplotype(toy_annotation, parse_svs):
    # v0.1 kept the first SV seen
    _, by_gene = run(toy_annotation, parse_svs, [
        DEL("I", 105_500, 106_000, hp="1|1"),            # intronic, MODERATE
        DEL("X", E2[0] + 10, E2[0] + 20, hp="1|1"),      # coding exon, HIGH
        DEL("U", 132_000, 133_000, hp="0|0"),            # promoter/UTR, MODERATE, unphased
    ])
    g = next(r for r in by_gene.values() if r.Gene_symbol == "GENEA")
    assert (g.impact(1), g.impact(2), g.impact(0)) == ("HIGH", "", "MODERATE")


def test_pair_haplotypes_follow_breakends(toy_annotation, parse_svs):
    svs, _ = run(toy_annotation, parse_svs, [BND("F", 105_500, ("chr2", 105_000), hp="1|2")])
    ps = pairs(svs[0])
    assert ps["GENEA"].haplotypes == {1}
    assert ps["GENEG"].haplotypes == {2}


def test_outputs(toy_annotation, parse_svs, tmp_path):
    svs, by_gene = run(toy_annotation, parse_svs, [
        BND("F", 105_500, ("chr2", 105_000), hp="1|1"), DEL("W", 290_000, 336_000, hp="2|2")])
    annot.output_svs(svs, str(tmp_path))
    annot.output_genes(by_gene, str(tmp_path))
    annot.output_gene_sv_pairs(svs, str(tmp_path))

    rows = [l.rstrip("\n").split("\t") for l in open(tmp_path / "gene_sv_pairs.tsv")]
    head, rows = rows[0], [dict(zip(rows[0], r)) for r in rows[1:]]
    assert head == annot.GENE_SV_PAIR_COLUMNS
    a = next(r for r in rows if r["gene"] == "GENEA")
    assert (a["class"], a["tier"], a["bp1_location"], a["fusion"]) == ("possible_fusion", "HIGH", "intron3", "GENEA::GENEG")
    d = next(r for r in rows if r["gene"] == "GENED")
    assert (d["class"], d["tier"], d["labels"], d["haplotype"]) == ("", "", "whole_gene_deletion", "2")

    svs_out = [l.rstrip("\n").split("\t") for l in open(tmp_path / "annotated_svs.tsv")]
    assert "Score" not in svs_out[0]
    assert all(len(r) == len(svs_out[0]) for r in svs_out)
    genes_out = [l.rstrip("\n").split("\t") for l in open(tmp_path / "by_gene.tsv")]
    assert "Score" not in genes_out[0]
    assert all(len(r) == len(genes_out[0]) for r in genes_out)
    gd = dict(zip(genes_out[0], next(r for r in genes_out if r[0] == "GENED")))
    assert gd["HP2 SV"] == "-:[W,chr1:290000-336000,whole_gene_deletion]"
