"""Compact annotation: build from GFF3, transcript selection and breakend lookup."""

import gzip

import pytest

from padfoot.annotation import FORMAT_HEADER, is_compact_annotation
from toy_genome import exon


def gene(annotation, name):
    return next(g for g in annotation.genes if g.name == name)


def hit(annotation, chrom, pos, name):
    return next(h for h in annotation.lookup(chrom, pos) if h.gene.name == name)


# --- build -------------------------------------------------------------------

def test_header_records_source(toy_table):
    with gzip.open(toy_table, "rt") as f:
        first = f.readline()
    assert first.startswith(FORMAT_HEADER)
    assert "source=toy.gff3" in first and "md5=" in first and "built=" in first
    assert is_compact_annotation(toy_table)


def test_only_protein_coding_genes(toy_annotation):
    names = {g.name for g in toy_annotation.genes}
    assert names == {"GENEA", "GENEB", "GENEC", "GENED", "GENEE", "GENEF", "GENEG"}


def test_transcript_selection(toy_annotation):
    # MANE + CCDS kept, untagged transcript dropped
    a = gene(toy_annotation, "GENEA")
    assert len(a.transcripts) == 2
    assert "MANE_Select" in a.naming.tags


def test_untagged_gene_falls_back_to_longest(toy_annotation):
    g = gene(toy_annotation, "GENEG")
    assert len(g.transcripts) == 1
    assert g.naming.end == 109_999


def test_gene_body_and_padding(toy_annotation):
    a = gene(toy_annotation, "GENEA")
    assert (a.start, a.end) == (100_000, 128_999)
    assert toy_annotation.lookup("chr1", 90_000)        # start - 10 kb
    assert not toy_annotation.lookup("chr1", 89_999)
    assert toy_annotation.lookup("chr1", 138_999)       # end + 10 kb
    assert not toy_annotation.lookup("chr1", 139_000)


def test_chrom_without_prefix(toy_annotation):
    assert [h.gene.name for h in toy_annotation.lookup("2", 105_000)] == ["GENEG"]


# --- breakend labels -----------------------------------------------------------

@pytest.mark.parametrize("pos,label", [
    (exon("GENEA", 1)[0], "exon1"),
    (exon("GENEA", 13)[1], "exon13"),
    (105_000, "intron3"),
    (108_000, "intron4"),          # v0.1: promoter/utr (len(exons) slip)
    (127_000, "intron12"),         # v0.1: intron2 (last digit only)
    (92_000, "promoter/utr"),
    (132_000, "promoter/utr"),     # v0.1: intron3 (right padding)
])
def test_plus_strand_labels(toy_annotation, pos, label):
    assert hit(toy_annotation, "chr1", pos, "GENEA").label == label


def test_minus_strand_numbers_follow_transcript(toy_annotation):
    # GENEB (-): exon1 is right-most, intron1 lies between exon1 and exon2
    assert hit(toy_annotation, "chr1", exon("GENEB", 1)[0] + 10, "GENEB").label == "exon1"
    assert hit(toy_annotation, "chr1", 215_000, "GENEB").label == "intron1"
    assert hit(toy_annotation, "chr1", 205_000, "GENEB").label == "intron2"


def test_flank_sides_are_strand_aware(toy_annotation):
    assert hit(toy_annotation, "chr1", 92_000, "GENEA").side == "upstream"
    assert hit(toy_annotation, "chr1", 132_000, "GENEA").side == "downstream"
    assert hit(toy_annotation, "chr1", 225_000, "GENEB").side == "upstream"
    assert hit(toy_annotation, "chr1", 195_000, "GENEB").side == "downstream"


def test_in_body(toy_annotation):
    assert hit(toy_annotation, "chr1", 105_000, "GENEA").in_body
    assert not hit(toy_annotation, "chr1", 92_000, "GENEA").in_body


def test_coding_and_utr_only_exons(toy_annotation):
    assert hit(toy_annotation, "chr1", exon("GENEA", 2)[0], "GENEA").coding
    d1 = hit(toy_annotation, "chr1", exon("GENED", 1)[0] + 10, "GENED")
    assert d1.label == "exon1" and not d1.coding


def test_alt_isoform_exon(toy_annotation):
    h = hit(toy_annotation, "chr1", 106_600, "GENEA")       # exon of GENEA-202 (CCDS)
    assert h.label == "intron3" and h.alt_isoform_exon
    h = hit(toy_annotation, "chr1", 110_600, "GENEA")       # exon of untagged GENEA-203
    assert h.label == "intron5" and not h.alt_isoform_exon


# --- overlapping genes ---------------------------------------------------------

def test_overlapping_padded_genes_are_all_returned(toy_annotation):
    assert [h.gene.name for h in toy_annotation.lookup("chr1", 312_000)] == ["GENEC", "GENED"]
    assert [h.gene.name for h in toy_annotation.lookup("chr1", 316_000)] == ["GENEC", "GENED"]


def test_nested_gene(toy_annotation):
    hits = {h.gene.name: h for h in toy_annotation.lookup("chr1", 421_000)}
    assert set(hits) == {"GENEE", "GENEF"}
    assert hits["GENEE"].label == "intron1"
    assert hits["GENEF"].label == "intron1"
