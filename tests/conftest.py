import os
import sys
from collections import defaultdict

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from padfoot import annot                                      # noqa: E402
from padfoot.annotation import build_annotation, load_annotation  # noqa: E402
import toy_genome                                              # noqa: E402
import severus_vcf                                             # noqa: E402


@pytest.fixture(scope="session")
def toy_table(tmp_path_factory):
    """The toy GFF3 converted to Padfoot's compact annotation table."""
    d = tmp_path_factory.mktemp("annot")
    gff = toy_genome.write_gff3(d / "toy.gff3")
    out = d / "toy_genes.tsv.gz"
    build_annotation(str(gff), str(out))
    return str(out)


@pytest.fixture(scope="session")
def toy_annotation(toy_table):
    return load_annotation(toy_table)


@pytest.fixture
def parse_svs(tmp_path):
    """Write Severus-style records to a VCF and parse them with get_bps."""
    def _parse(records):
        path = severus_vcf.write_vcf(tmp_path / "severus.vcf", records)
        return annot.get_bps(str(path))
    return _parse


@pytest.fixture
def classify(toy_annotation, parse_svs):
    """Run SV-level classification on the toy genome; return {SV id: SV}."""
    def _classify(records):
        svs = parse_svs(records)
        annot.annot_SVS(toy_annotation, svs, defaultdict(list))
        return {sv.vcf_id: sv for sv in svs}
    return _classify
