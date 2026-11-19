"""Known-fusion lookup (annot.cancer_annot_svs)."""

from collections import defaultdict

import pytest

from padfoot import annot

FUSION_BND = dict(id="F", chrom="chr1", pos=105_500, mate=("chr2", 105_000), svtype="BND", strands="+-")


@pytest.mark.xfail(strict=True, raises=AssertionError,
                   reason="v0.1 bug: oncogenic_fusion never assigned (== instead of =, wrong key)")
@pytest.mark.parametrize("fusions", [{"GENEA": ["GENEG"]}, {"GENEG": ["GENEA"]}], ids=["listed-5p", "listed-3p"])
def test_known_fusion(classify, fusions):
    sv = classify([FUSION_BND])["F_1"]
    assert sv.genes[2] == "possible_fusion"
    annot.cancer_annot_svs([sv], cancer_genes={}, fusion=defaultdict(list, fusions))
    assert sv.genes[2] == "oncogenic_fusion"
