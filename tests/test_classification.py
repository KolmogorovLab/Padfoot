"""
SV consequence classes on the toy genome (tests/toy_genome.py).

Each case states the class expected under the agreed v0.2 rules (same classes
as v0.1, bugs fixed). Cases that v0.1 gets wrong are marked xfail(strict=True)
with the bug as the reason: when a fix lands, the case starts passing, the
strict xfail turns into a failure, and the marker must be removed in the same
commit.
"""

import pytest

from toy_genome import exon

V01_BUG = lambda reason: pytest.mark.xfail(strict=True, raises=AssertionError, reason="v0.1 bug: " + reason)  # noqa: E731

E2 = exon("GENEA", 2)          # 102,400-102,599


def DEL(sid, chrom, start, end):
    return dict(id=sid, chrom=chrom, pos=start, end=end, svlen=end - start, svtype="DEL")


CASES = [
    # --- classes v0.1 already gets right ------------------------------------
    pytest.param(DEL("A", "chr1", 105_500, 106_000), "intronic", id="intronic"),
    pytest.param(DEL("A", "chr1", E2[0] + 50, E2[0] + 100), "frameshift", id="exon-frameshift"),
    pytest.param(DEL("A", "chr1", 92_000, 93_000), "promoter/UTR", id="left-flank"),
    pytest.param(DEL("A", "chr1", 105_500, 160_000), "Gene-NonCoding", id="gene-to-nogene"),
    pytest.param(DEL("A", "chr1", 105_500, 205_000), "between_genes", id="between-genes"),
    pytest.param(DEL("A", "chr1", 103_000, 106_000), "between_exons", id="intron2-to-intron3"),
    pytest.param(DEL("A", "chr1", 150_000, 160_000), "", id="no-gene"),
    pytest.param(dict(id="A", chrom="chr1", pos=105_500, mate=("chr2", 105_000), svtype="BND", strands="+-"),
                 "possible_fusion", id="fusion-BND"),
    pytest.param(dict(id="A", chrom="chr1", pos=E2[0] + 50, svtype="sBND"), "exon_altering", id="sBND-in-exon"),

    # --- v0.1 bugs -----------------------------------------------------------
    pytest.param(DEL("A", "chr1", E2[0] + 50, E2[0] + 101), "exon_altering", id="exon-inframe",
                 marks=V01_BUG("frame test uses svlen // 3 instead of svlen % 3")),
    pytest.param(dict(id="A", chrom="chr1", pos=E2[0] + 50, svlen=50, svtype="INS", alt="A" * 50),
                 "frameshift", id="INS-frameshift",
                 marks=V01_BUG("insertions get length 0; SVLEN should be used")),
    pytest.param(DEL("A", "chr1", 103_000, 127_000), "between_exons", id="intron2-to-intron12",
                 marks=V01_BUG("intron label keeps only the last digit (intron12 -> intron2)")),
    pytest.param(DEL("A", "chr1", 108_000, 108_500), "intronic", id="intron4",
                 marks=V01_BUG("len(exons) instead of len(exons[0]): intron 4 / exon 5 labelled promoter/utr")),
    pytest.param(DEL("A", "chr1", 132_000, 133_000), "promoter/UTR", id="right-flank",
                 marks=V01_BUG("right-hand padding labelled as intron")),
    pytest.param(DEL("A", "chr1", 316_000, 318_000), "intronic", id="overlapping-padded-genes",
                 marks=V01_BUG("bisect lookup misses genes whose padded intervals overlap")),
    pytest.param(dict(id="A", chrom="chr1", pos=105_500, mate=("chr2", 95_000), svtype="BND", strands="+-"),
                 "between_genes", id="fusion-needs-gene-body",
                 marks=V01_BUG("possible_fusion called when one end is only in the padding")),
]


@pytest.mark.parametrize("record,expected", CASES)
def test_class(classify, record, expected):
    sv = next(iter(classify([record]).values()))
    assert sv.genes[2] == expected


def test_breakend_labels(classify):
    sv = classify([DEL("A", "chr1", 105_500, 106_000)])["A"]
    assert sv.genes[0] == ("GENEA", "intron3", "+")


def test_minus_strand_exon_numbers(classify):
    # GENEB is on the minus strand: its right-most exon is exon 1
    s, _ = exon("GENEB", 1)
    sv = classify([DEL("B", "chr1", s + 10, s + 20)])["B"]
    assert sv.genes[0] == ("GENEB", "exon1", "-")
