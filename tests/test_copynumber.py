"""Wakhan integer profile reader and gene CN fields."""

from collections import defaultdict

import pytest

from padfoot import annot
from padfoot.copynumber import (read_wakhan_profile, gene_cn, gene_cn_fields, coverage_gaps,
                                GENE_CN_COLUMNS)
from wakhan_profile import write_profile, COLUMNS

M = 1_000_000

# 100 kb blocks on chr1; total CN 2 (1+1) covers most autosomal length, so modal CN = 2.
ROWS = [
    ("chr1", 0, 1 * M - 1, 1, 1),
    ("chr1", 1_000_000, 1_099_999, 0, 0),             # HOMDEL
    ("chr1", 1_100_000, 1_199_999, 1, 0),             # LOH with loss
    ("chr1", 1_200_000, 1_299_999, 2, 0),             # LOH copy-neutral
    ("chr1", 1_300_000, 1_399_999, 3, 0),             # LOH with gain
    ("chr1", 1_400_000, 1_499_999, 2, 1),             # total 3: GAIN
    ("chr1", 1_500_000, 1_599_999, 4, 2),             # total 6 = 3 x modal: still GAIN
    ("chr1", 1_600_000, 1_699_999, 5, 2),             # total 7 > 3 x modal: GAIN_HIGH
    ("chr1", 1_700_000, 1_799_999, (26, 0.0), 1),     # amplicon with confidence 0 (Wakhan < 0.5)
    ("chr1", 1_800_000, 1_899_999, None, 1),          # HP1 zero-filled by the merge step
    ("chr1", 1_900_000, 99 * M - 1, 1, 1),
    ("chr2", 0, 10 * M - 1, 2, 2),
    ("chr2", 10_500_000, 20 * M - 1, 2, 2),           # gap 10.0-10.5 Mb
    ("chrX", 0, 1 * M - 1, 1, 0),
    ("chrX", 1 * M, 200 * M - 1, 2, 2),               # longest, but not autosomal
    ("chrY", 0, 1 * M - 1, 1, 0),
]


@pytest.fixture
def profile(tmp_path):
    return read_wakhan_profile(str(write_profile(tmp_path / "integer_profile.bed", ROWS)))


def block(k):
    """A gene inside the k-th 100 kb block after 1 Mb on chr1."""
    s = 1 * M + k * 100_000
    return ("chr1", s + 10_000, s + 20_000)


# --- reader --------------------------------------------------------------------

def test_reads_profile(tmp_path):
    p = read_wakhan_profile(str(write_profile(tmp_path / "integer_profile.bed", ROWS)))
    assert len(p) == len(ROWS) and p.path.endswith("integer_profile.bed")


def test_columns_come_from_header(tmp_path):
    shuffled = COLUMNS[:3] + COLUMNS[6:9] + COLUMNS[3:6] + COLUMNS[9:]
    p = read_wakhan_profile(str(write_profile(tmp_path / "p.bed", ROWS, columns=shuffled)))
    seg = p.overlapping(*block(4))[0]
    assert (seg.hp1, seg.hp2) == (2, 1)


def test_old_per_haplotype_bed_is_rejected(tmp_path):
    path = tmp_path / "old_HP_1.bed"
    path.write_text("#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tsvs_breakpoints_ids\n"
                    "chr1\t0\t100\t8.0\t1.0\t0.9\t\n")
    with pytest.raises(ValueError, match="missing columns"):
        read_wakhan_profile(str(path))


def test_zero_filled_haplotype_is_missing(profile):
    seg = profile.overlapping("chr1", 1_850_000, 1_850_000)[0]
    assert seg.hp1 is None and seg.hp2 == 1 and seg.total is None


def test_confidence_zero_alone_is_not_missing(profile):
    seg = profile.overlapping("chr1", 1_750_000, 1_750_000)[0]
    assert seg.hp1 == 26


def test_modal_cn_is_autosomal(profile):
    assert profile.modal_cn == 2


def test_chrom_without_prefix(profile):
    assert profile.overlapping("1", 500, 600)


# --- coverage check --------------------------------------------------------------

def test_coverage_gaps(profile):
    lengths = {"chr1": 99 * M, "chr2": 20 * M + 50_000, "chrX": 200 * M, "chrY": 1 * M, "chr3": 5 * M}
    absent, gaps = coverage_gaps(profile, lengths, ["chr1", "chr2", "chr3", "chrX", "chrY"])
    assert absent == ["chr3"]
    assert ("chr2", 10 * M, 10_500_000 - 1) in gaps           # internal gap
    assert ("chr1", 1_800_000, 1_899_999) in gaps             # zero-filled fragment
    assert not any(g[0] == "chr2" and g[1] >= 20 * M for g in gaps)   # end within tolerance


def test_chromosome_end_gap_beyond_tolerance(profile):
    _, gaps = coverage_gaps(profile, {"chr1": 100 * M}, ["chr1"])
    assert ("chr1", 99 * M, 100 * M) in gaps


# --- gene CN status --------------------------------------------------------------

@pytest.mark.parametrize("k,status,total", [
    (-5, "NEUTRAL", (2, 2)),
    (0, "HOMDEL", (0, 0)),
    (1, "LOSS", (1, 1)),
    (2, "NEUTRAL", (2, 2)),
    (3, "GAIN", (3, 3)),
    (4, "GAIN", (3, 3)),
    (5, "GAIN", (6, 6)),
    (6, "GAIN_HIGH", (7, 7)),
    (7, "GAIN_HIGH", (27, 27)),
])
def test_cn_status(profile, k, status, total):
    g = gene_cn(profile, *block(k))
    assert (g.status, g.total, g.modal, g.boundary, g.missing) == (status, total, 2, False, False)


def test_gain_high_multiplier(profile):
    assert gene_cn(profile, *block(5), gain_high_multiplier=2.5).status == "GAIN_HIGH"


@pytest.mark.parametrize("k,loh_type", [(1, "with_loss"), (2, "copy_neutral"), (3, "with_gain")])
def test_loh(profile, k, loh_type):
    g = gene_cn(profile, *block(k))
    assert g.loh and g.loh_type == loh_type


def test_homdel_is_not_loh(profile):
    g = gene_cn(profile, *block(0))
    assert not g.loh and g.loh_type == ""


def test_haplotype_min_max(profile):
    g = gene_cn(profile, "chr1", 1_350_000, 1_450_000)        # (3,0) then (2,1)
    assert (g.hp1, g.hp2, g.total) == ((2, 3), (0, 1), (3, 3))
    assert g.boundary


@pytest.mark.parametrize("start,end,status", [
    (1_050_000, 1_650_000, "HOMDEL"),       # HOMDEL over GAIN_HIGH
    (1_450_000, 1_650_000, "GAIN_HIGH"),    # GAIN_HIGH from the maximum
    (1_150_000, 1_450_000, "LOSS"),         # LOSS over GAIN
    (1_450_000, 1_550_000, "GAIN"),
])
def test_cn_boundary_in_gene(profile, start, end, status):
    g = gene_cn(profile, "chr1", start, end)
    assert g.status == status and g.boundary


def test_loh_type_uses_lowest_total(profile):
    g = gene_cn(profile, "chr1", 1_150_000, 1_350_000)        # LOH with loss, copy-neutral, with gain
    assert g.loh_type == "with_loss"


def test_partly_missing_gene(profile):
    g = gene_cn(profile, "chr1", 1_850_000, 1_950_000)
    assert g.missing and g.status == "NEUTRAL" and g.total == (2, 2) and g.hp2 == (1, 1)


def test_fully_missing_gene(profile):
    g = gene_cn(profile, "chr1", 1_810_000, 1_820_000)
    assert g.missing and g.status == "" and g.total == (None, None) and not g.loh


def test_gene_in_gap(profile):
    g = gene_cn(profile, "chr2", 10_100_000, 10_200_000)
    assert g.missing and g.status == ""


def test_gene_on_unprofiled_chromosome(profile):
    g = gene_cn(profile, "chr3", 100, 200)
    assert g.missing and g.status == ""


@pytest.mark.parametrize("chrom,sex,loh", [
    ("chrX", "female", True),
    ("chrX", "unknown", True),
    ("chrX", "male", False),
    ("chrY", "male", False),
    ("chrY", "female", False),
    ("chr1", "male", True),
])
def test_sex_chromosome_loh(profile, chrom, sex, loh):
    start = 1_100_000 if chrom == "chr1" else 10_000
    g = gene_cn(profile, chrom, start + 10_000, start + 20_000, sex=sex)
    assert g.loh == loh
    if chrom == "chrX" and sex == "male":
        assert g.status == "LOSS"                              # CN status is still called


def test_fields(profile):
    vals = dict(zip(GENE_CN_COLUMNS, gene_cn_fields(gene_cn(profile, *block(1)))))
    assert vals == {"HP1_CN_min": "1", "HP1_CN_max": "1", "HP2_CN_min": "0", "HP2_CN_max": "0",
                    "total_CN_min": "1", "total_CN_max": "1", "modal_CN": "2", "cn_status": "LOSS",
                    "LOH": "yes", "LOH_type": "with_loss", "cn_boundary_in_gene": "no", "cn_missing": "no"}
    assert gene_cn_fields(None) == [""] * len(GENE_CN_COLUMNS)


# --- genes in the toy annotation -----------------------------------------------

TOY_ROWS = [
    ("chr1", 0, 199_999, 1, 1),
    ("chr1", 200_000, 299_999, 0, 0),            # GENEB (200,000-219,999)
    ("chr1", 300_000, 312_000, 1, 1),            # GENEC (300,000-309,999)
    ("chr1", 312_001, 999_999, 2, 1),            # GENED (315,000-324,999): GAIN
    ("chr2", 0, 999_999, 1, 1),
]


def test_by_gene_cn(toy_annotation, tmp_path):
    profile = read_wakhan_profile(str(write_profile(tmp_path / "integer_profile.bed", TOY_ROWS)))
    by_gene = defaultdict(list)
    annot.annot_gene_cn(toy_annotation, profile, by_gene)
    assert set(by_gene) == {g.name for g in toy_annotation.genes}
    assert by_gene["GENEB"].cn.status == "HOMDEL"
    assert by_gene["GENED"].cn.status == "GAIN"
    assert by_gene["GENEA"].cn.status == "NEUTRAL"
    annot.output_genes(by_gene, str(tmp_path))
    rows = [l.rstrip("\n").split("\t") for l in open(tmp_path / "by_gene.tsv")]
    assert rows[0][-len(GENE_CN_COLUMNS):] == GENE_CN_COLUMNS
    assert all(len(r) == len(rows[0]) for r in rows)
    b = dict(zip(rows[0], next(r for r in rows if r[0] == "GENEB")))
    assert (b["cn_status"], b["total_CN_min"], b["HP1 Impact"]) == ("HOMDEL", "0", "")
