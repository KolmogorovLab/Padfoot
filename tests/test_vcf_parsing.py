"""Severus VCF parsing (annot.get_bps)."""


def by_id(svs):
    return {sv.vcf_id: sv for sv in svs}


def test_del_uses_end(parse_svs):
    sv = by_id(parse_svs([dict(id="DEL1", chrom="chr1", pos=1000, end=5000, svlen=4000, svtype="DEL")]))["DEL1"]
    assert sv.bp_1 == ("chr1", 1000)
    assert sv.bp_2 == ("chr1", 5000)
    assert (sv.direction_1, sv.direction_2) == ("+", "-")


def test_ins_has_single_position(parse_svs):
    sv = by_id(parse_svs([dict(id="INS1", chrom="chr1", pos=1000, svlen=50, svtype="INS", alt="A" * 50)]))["INS1"]
    assert sv.bp_1 == sv.bp_2 == ("chr1", 1000)
    assert sv.ins_seq == "A" * 50


def test_sbnd_has_single_position(parse_svs):
    sv = by_id(parse_svs([dict(id="SB1", chrom="chr1", pos=1000, svtype="sBND")]))["SB1"]
    assert sv.is_single
    assert sv.bp_1 == sv.bp_2


def test_bnd_pair_becomes_one_sv(parse_svs):
    svs = parse_svs([dict(id="BND1", chrom="chr1", pos=1000, mate=("chr2", 2000), svtype="BND", strands="+-")])
    assert [sv.vcf_id for sv in svs] == ["BND1_1"]
    sv = svs[0]
    assert sv.bp_2 == ("chr2", 2000)
    assert (sv.direction_1, sv.direction_2) == ("+", "-")


def test_bnd_haplotypes_come_from_first_record(parse_svs):
    # Severus writes HP=<hp at bp1>|<hp at bp2> in the _1 record (and the
    # reverse in the _2 record), so the _1 record alone is enough.
    sv = parse_svs([dict(id="BND1", chrom="chr1", pos=1000, mate=("chr2", 2000), svtype="BND", hp="2|0")])[0]
    assert (sv.hp1, sv.hp2) == (2, 0)


def test_cluster_and_detailed_type_are_kept(parse_svs):
    sv = parse_svs([dict(id="BND1", chrom="chr1", pos=1000, mate=("chr2", 2000), svtype="BND",
                         cluster="severus_7", detailed_type="complex_inv")])[0]
    assert sv.cluster_id == "severus_7"
    assert sv.detailed_type == "complex_inv"
