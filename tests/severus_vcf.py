"""
Minimal writer for Severus-style somatic VCFs (header copied from Severus 1.7).

Each record is a dict with keys:
  id, chrom, pos, svtype            (required)
  end, svlen, strands, hp, mate     (mate = (chrom, pos) for BND)
  vaf, dv, cluster, detailed_type, alt, vntr
BND records are written as a _1/_2 pair, like Severus does; HP in the _2 record
is given in mate order (bp2|bp1).
"""

HEADER = """##fileformat=VCFv4.2
##source=Severus_v1.7
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Reciprocal Inversion">
##ALT=<ID=BND,Description="Breakend">
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="SV with precise breakpoints coordinates and length">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the SV">
##INFO=<ID=STRANDS,Number=1,Type=String,Description="Breakpoint strandedness">
##INFO=<ID=DETAILED_TYPE,Number=1,Type=String,Description="Detailed type of the SV">
##INFO=<ID=HP,Number=1,Type=String,Description="Matching haplotype ID for phased SVs">
##INFO=<ID=CLUSTERID,Number=1,Type=String,Description="Cluster ID in breakpoint_graph">
##INFO=<ID=INSSEQ,Number=1,Type=String,Description="Insertion sequence between breakpoints">
##INFO=<ID=MATE_ID,Number=1,Type=String,Description="MATE ID for breakends">
##INFO=<ID=INSIDE_VNTR,Number=1,Type=String,Description="True if an indel is inside a VNTR">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thaplotagged
"""


def _bnd_alt(strands, mate):
    """VCF breakend ALT for the first record of a pair (Severus convention)."""
    c, p = mate
    return {"+-": f"N[{c}:{p}[", "++": f"N]{c}:{p}]",
            "-+": f"]{c}:{p}]N", "--": f"[{c}:{p}[N"}[strands]


def _line(chrom, pos, vid, alt, info, vaf, dv):
    fmt = f"0/1:{dv}:{dv}:{vaf}"
    return "\t".join([chrom, str(pos), vid, "N", alt, "60", "PASS", ";".join(info), "GT:DR:DV:VAF", fmt])


def write_vcf(path, records):
    lines = []
    for r in records:
        svtype = r["svtype"]
        vaf, dv = r.get("vaf", 0.4), r.get("dv", 10)
        info = ["PRECISE", f"SVTYPE={svtype}"]
        if "svlen" in r:
            info.append(f"SVLEN={r['svlen']}")
        if "end" in r:
            info.append(f"END={r['end']}")
        strands = r.get("strands", {"DEL": "+-", "DUP": "-+", "INV": "++", "INS": None,
                                    "sBND": "-", "BND": "+-"}[svtype])
        if strands:
            info.append(f"STRANDS={strands}")
        if "detailed_type" in r:
            info.append(f"DETAILED_TYPE={r['detailed_type']}")
        if r.get("vntr"):
            info.append("INSIDE_VNTR=TRUE")
        hp = r.get("hp", "0|0")                         # None: no HP field (unphased in Severus)
        cluster = [f"CLUSTERID={r['cluster']}"] if "cluster" in r else []

        if svtype == "BND":
            mchrom, mpos = r["mate"]
            hp1 = [f"HP={hp}"] if hp else []
            hp2 = ["HP={1}|{0}".format(*hp.split("|"))] if hp else []
            lines.append(_line(r["chrom"], r["pos"], r["id"] + "_1", _bnd_alt(strands, r["mate"]),
                               info + [f"MATE_ID={r['id']}_2"] + hp1 + cluster, vaf, dv))
            lines.append(_line(mchrom, mpos, r["id"] + "_2", _bnd_alt(strands, (r["chrom"], r["pos"])),
                               info + [f"MATE_ID={r['id']}_1"] + hp2 + cluster, vaf, dv))
        else:
            alt = r.get("alt", "N" if svtype == "sBND" else f"<{svtype}>")
            lines.append(_line(r["chrom"], r["pos"], r["id"], alt, info + ([f"HP={hp}"] if hp else []) + cluster, vaf, dv))

    with open(path, "w") as f:
        f.write(HEADER)
        f.write("\n".join(lines) + "\n")
    return path
