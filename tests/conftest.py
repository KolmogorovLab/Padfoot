import os
import sys
from collections import defaultdict

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from padfoot import annot          # noqa: E402
import toy_genome                  # noqa: E402
import severus_vcf                 # noqa: E402


@pytest.fixture(scope="session")
def v01_annotation(tmp_path_factory):
    path = tmp_path_factory.mktemp("annot") / "toy.gff3.gz"
    return str(toy_genome.write_v01_annotation(path))


@pytest.fixture
def parse_svs(tmp_path):
    """Write Severus-style records to a VCF and parse them with get_bps."""
    def _parse(records):
        path = severus_vcf.write_vcf(tmp_path / "severus.vcf", records)
        return annot.get_bps(str(path))
    return _parse


@pytest.fixture
def classify(v01_annotation, parse_svs):
    """Run SV-level classification on the toy genome; return {SV id: SV}."""
    def _classify(records):
        genes, exon_pos = annot.get_genes(v01_annotation)
        svs = parse_svs(records)
        by_gene = defaultdict(list)
        annot.annot_SVS(genes, exon_pos, svs, by_gene)
        return {sv.vcf_id: sv for sv in svs}
    return _classify
