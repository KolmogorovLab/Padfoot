"""Cancer gene list: roles and the known-fusion lookup (labels only)."""

import os
from collections import defaultdict

import pytest

from padfoot import annot

FUSION_BND = dict(id="F", chrom="chr1", pos=105_500, mate=("chr2", 105_000), svtype="BND", strands="+-")
BUNDLED = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "beds", "hg38_cancer_genes.tsv")


def write_list(path, rows, header=("original_symbol", "Gene_symbol", "Role_in_cancer", "fusion", "Supp", "nSupp")):
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    return str(path)


def run(classify, records, cancer_genes, fusion):
    svs = classify(records)
    annot.cancer_annot_svs(list(svs.values()), cancer_genes, fusion)
    return svs


# --- known fusions -------------------------------------------------------------

@pytest.mark.parametrize("fusion", [{"GENEA": ["GENEG"]}, {"GENEG": ["GENEA"]}], ids=["listed-5p", "listed-3p"])
def test_known_fusion(classify, fusion):
    sv = run(classify, [FUSION_BND], {}, fusion)["F_1"]
    assert sv.genes[2] == "oncogenic_fusion"
    assert sv.fusion == "GENEA::GENEG"
    assert {p.gene.name: (p.cls, p.tier) for p in sv.pairs} == {
        "GENEA": ("oncogenic_fusion", "HIGH"), "GENEG": ("oncogenic_fusion", "HIGH")}
    assert sv.impact == "HIGH"


def test_unlisted_pair_stays_possible_fusion(classify):
    sv = run(classify, [FUSION_BND], {}, {"GENEA": ["GENEB"]})["F_1"]
    assert sv.genes[2] == "possible_fusion"


def test_known_pair_without_compatible_orientation_is_not_a_fusion(classify):
    sv = run(classify, [dict(FUSION_BND, strands="++")], {}, {"GENEA": ["GENEG"]})["F_1"]
    assert sv.genes[2] == "between_genes"


def test_roles_are_labels_only(classify):
    records = [FUSION_BND, dict(id="D", chrom="chr1", pos=105_000, end=106_000, svtype="DEL", strands="+-")]
    plain = {k: (sv.genes[2], sv.impact, [p.tier for p in sv.pairs]) for k, sv in classify(records).items()}
    svs = run(classify, records, {"GENEA": "TSG", "GENEG": "oncogene"}, {})
    assert {k: (sv.genes[2], sv.impact, [p.tier for p in sv.pairs]) for k, sv in svs.items()} == plain
    assert svs["F_1"].cancer == ["TSG", "oncogene"]
    assert svs["D"].cancer == ["TSG", "TSG"]


def test_gene_table_role(classify):
    svs = classify([FUSION_BND])
    by_gene = defaultdict(list)
    for sv in svs.values():
        for p in sv.pairs:
            annot.add_gene(by_gene, p.gene.name, p.gene.chrom, p.gene.start, p.gene.end)
    annot.cancer_annot_genes(by_gene, {"GENEG": "oncogene"})
    assert (by_gene["GENEG"].cancer, by_gene["GENEA"].cancer) == ("oncogene", "")


# --- list loading ----------------------------------------------------------------

def test_load_cancer_genes(tmp_path):
    path = write_list(tmp_path / "c.tsv", [
        ("ABL1", "ABL1", "oncogene", "BCR, ETV6", "COSMIC", "1"),
        ("BCR", "BCR", "", "ABL1", "COSMIC", "1"),
        ("X", "X", "TSG", "", "COSMIC", "1"),
        ("", "", "oncogene", "Y", "", "0"),
    ])
    roles, fusion = annot.load_cancer_genes(path)
    assert roles == {"ABL1": "oncogene", "X": "TSG"}
    assert fusion["ABL1"] == ["BCR", "ETV6"] and fusion["BCR"] == ["ABL1"] and not fusion["X"]
    assert annot.is_known_fusion("ETV6", "ABL1", fusion)
    assert not annot.is_known_fusion("ETV6", "BCR", fusion)


def test_column_order_and_extra_columns(tmp_path):
    path = write_list(tmp_path / "c.tsv", [("a", "ABL1", "", "oncogene", "BCR")],
                      header=("note", "Gene_symbol", "fusion", "Role_in_cancer", "extra"))
    roles, fusion = annot.load_cancer_genes(path)
    assert roles == {"ABL1": "oncogene"} and fusion["ABL1"] == []


def test_missing_column(tmp_path):
    path = write_list(tmp_path / "c.tsv", [("ABL1", "oncogene")], header=("Gene_symbol", "Role_in_cancer"))
    with pytest.raises(ValueError, match="fusion"):
        annot.load_cancer_genes(path)


def test_bundled_list():
    roles, fusion = annot.load_cancer_genes(BUNDLED)
    assert len(roles) == 974
    assert roles["ERBB2"] == "oncogene" and roles["PTEN"] == "TSG"
    assert annot.is_known_fusion("BCR", "ABL1", fusion) and annot.is_known_fusion("ABL1", "BCR", fusion)
    assert "SEPTIN9" in fusion and not any(g.endswith("-Sep") for g in roles)
