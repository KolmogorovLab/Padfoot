"""Output tables: annotated_svs.tsv and by_gene.tsv columns."""

from collections import defaultdict

import pytest

from padfoot import annot
from padfoot.annotation import Gene as AnnotGene


def DEL(i, s, e, **kw):
    return dict(id=i, chrom="chr1", pos=s, end=e, svtype="DEL", strands="+-", **kw)


def run(toy_annotation, parse_svs, tmp_path, records, fusion=None):
    svs = parse_svs(records)
    by_gene = defaultdict(list)
    annot.annot_SVS(toy_annotation, svs, by_gene)
    if fusion is not None:
        annot.cancer_annot_svs(svs, {}, fusion)
    annot.output_svs(svs, str(tmp_path))
    annot.output_genes(by_gene, str(tmp_path))
    return read(tmp_path / "annotated_svs.tsv", "SV ID"), read(tmp_path / "by_gene.tsv", "Gene Name")


def read(path, key):
    lines = [l.rstrip("\n").split("\t") for l in open(path)]
    assert all(len(r) == len(lines[0]) for r in lines)
    return lines[0], {r[lines[0].index(key)]: dict(zip(lines[0], r)) for r in lines[1:]}


RECORDS = [
    dict(id="F", chrom="chr1", pos=105_500, mate=("chr2", 105_000), svtype="BND", strands="+-",
         hp="1|1", cluster="severus_7", detailed_type="Reciprocal_tra"),
    DEL("W", 290_000, 336_000, hp="2|2"),             # GENEC and GENED lie inside
    DEL("N", 420_500, 421_000, hp="1|1"),             # nested GENEF; GENEE overlaps both ends too
    DEL("A", 106_600, 108_000, hp="1|1"),             # bp1 in an exon of the CCDS isoform only
]


def test_sv_columns(toy_annotation, parse_svs, tmp_path):
    (head, svs), _ = run(toy_annotation, parse_svs, tmp_path, RECORDS, fusion={"GENEA": ["GENEG"]})
    assert head == annot.SV_COLUMNS
    f = svs["F_1"]
    assert (f["Type"], f["fusion_partners"], f["known_fusion"]) == ("oncogenic_fusion", "GENEA::GENEG", "yes")
    assert (f["cluster_id"], f["detailed_type"]) == ("severus_7", "Reciprocal_tra")
    w = svs["W"]
    assert w["whole_gene_labels"] == "whole_gene_deletion:GENEC,GENED"
    assert (w["known_fusion"], w["fusion_partners"], w["cluster_id"]) == ("no", "", "")
    n = svs["N"]
    assert (n["Gene Name1"], n["Other_genes1"], n["Other_genes2"]) == ("GENEF", "GENEE", "GENEE")
    assert svs["A"]["alt_isoform_exon"] == "GENEA"


@pytest.mark.parametrize("hp,label", [("1|1", "HP1"), ("2|2", "HP2"), ("0|1", "HP1"), ("2|0", "HP2"),
                                      ("0|0", "HP0"), (None, "HP0"), ("1|2", "HP1|HP2")])
def test_haplotype_column(toy_annotation, parse_svs, tmp_path, hp, label):
    (_, svs), _ = run(toy_annotation, parse_svs, tmp_path, [DEL("H", 420_500, 421_000, hp=hp)])
    assert svs["H"]["haplotype"] == label


def test_gene_columns(toy_annotation, parse_svs, tmp_path):
    _, (head, genes) = run(toy_annotation, parse_svs, tmp_path,
                           RECORDS + [DEL("I", 105_500, 105_600, hp="1|1")])
    assert head == annot.GENE_COLUMNS
    a = genes["GENEA"]
    assert (a["n_SVs"], a["in_cluster"]) == ("3", "yes")          # F, A, I
    assert a["gene_id"].startswith("ENSG")
    assert (genes["GENEC"]["n_SVs"], genes["GENEC"]["in_cluster"]) == ("1", "no")


def test_same_name_at_two_loci_gets_two_rows():
    by_gene = defaultdict(list)
    for gid, chrom in (("ENSG1", "chrX"), ("ENSG1_PAR_Y", "chrY")):
        g = AnnotGene(gid, "SHOX", chrom, "+")
        g.start, g.end = 100, 200
        annot.add_gene(by_gene, g)
    assert sorted((r.ref_id, r.gene_id) for r in by_gene.values()) == [("chrX", "ENSG1"), ("chrY", "ENSG1_PAR_Y")]
