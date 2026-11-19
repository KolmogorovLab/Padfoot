"""Command-line checks in main.py."""

import logging

import pytest

from padfoot.main import _check_genome, GENOME_CHR1


def fai(tmp_path, chr1, name="chr1"):
    path = tmp_path / "ref.fa.fai"
    path.write_text(f"{name}\t{chr1}\t6\t60\t61\nchr2\t1000\t100\t60\t61\n")
    return str(path)


@pytest.mark.parametrize("genome,chr1,name,warns", [
    ("hg38", GENOME_CHR1["hg38"], "chr1", False),
    ("hg38", GENOME_CHR1["hg38"], "1", False),
    ("hg38", GENOME_CHR1["chm13"], "chr1", True),
    ("mm10", GENOME_CHR1["hg38"], "chr1", True),
    ("hg38", 12345, "chr1", True),
    ("hg38", 12345, "chrA", True),       # no chr1
    ("canfam", 12345, "chr1", True),     # unknown genome
])
def test_check_genome(tmp_path, caplog, genome, chr1, name, warns):
    with caplog.at_level(logging.INFO):
        _check_genome(genome, fai(tmp_path, chr1, name))
    assert any(r.levelno == logging.WARNING for r in caplog.records) == warns
