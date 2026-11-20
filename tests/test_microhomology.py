"""Microhomology at SV junctions, read from the reference with pysam."""

import random

import pysam

from padfoot import annot

MH = "GATTACAG"


def write_ref(path, lower=False):
    rnd = random.Random(1)
    seq = [rnd.choice("ACGT") for _ in range(2_000)]
    for pos in (500, 1_500):                  # the same 8-mer starts at both breakends (1-based 500 and 1500)
        seq[pos - 1:pos - 1 + len(MH)] = MH.lower() if lower else MH
    with open(path, "w") as f:
        f.write(">chr1\n" + "".join(seq) + "\n")
    pysam.faidx(str(path))
    return str(path)


def run(parse_svs, ref, records):
    svs = parse_svs(records)
    annot.get_microhomology(svs, ref)
    return {sv.vcf_id: sv.microh for sv in svs}


DEL = dict(id="D", chrom="chr1", pos=500, end=1_500, svtype="DEL", strands="+-")


def test_shared_sequence_is_found(parse_svs, tmp_path):
    mh = run(parse_svs, write_ref(tmp_path / "ref.fa"), [DEL])["D"]
    assert len(mh) == 1 and mh[0].startswith(MH)        # random flanks may extend it by chance


def test_no_shared_sequence(parse_svs, tmp_path):
    mh = run(parse_svs, write_ref(tmp_path / "ref.fa"), [dict(DEL, pos=800, end=1_200)])["D"]
    assert not mh or len(mh[0]) < len(MH)


def test_soft_masked_reference(parse_svs, tmp_path):
    upper = run(parse_svs, write_ref(tmp_path / "a.fa"), [DEL])["D"]
    lower = run(parse_svs, write_ref(tmp_path / "b.fa", lower=True), [DEL])["D"]
    assert lower == upper


def test_contig_missing_from_reference(parse_svs, tmp_path, caplog):
    records = [dict(id="F", chrom="chr1", pos=500, mate=("chrZ", 1_000), svtype="BND", strands="+-"), DEL]
    mh = run(parse_svs, write_ref(tmp_path / "ref.fa"), records)
    assert mh["F_1"] == [] and mh["D"]                # later SVs are unaffected
    assert "chrZ" in caplog.text


def test_insertions_are_skipped(parse_svs, tmp_path):
    mh = run(parse_svs, write_ref(tmp_path / "ref.fa"),
             [dict(id="I", chrom="chr1", pos=500, svtype="INS", svlen=60, alt="A" * 61)])
    assert mh["I"] == []


def test_breakend_near_contig_end(parse_svs, tmp_path, caplog):
    records = [dict(id="E", chrom="chr1", pos=500, end=1_990, svtype="DEL", strands="+-"),
               dict(id="S", chrom="chr1", pos=10, end=1_500, svtype="DEL", strands="+-"), DEL]
    mh = run(parse_svs, write_ref(tmp_path / "ref.fa"), records)
    assert mh["E"] == [] and mh["S"] == [] and mh["D"]
    assert "E, S" in caplog.text
