"""
Minimal writer for Wakhan >= 0.5.0 ``integer_profile.bed`` files.

Rows are (chrom, start, end, hp1, hp2) with closed coordinates, as Wakhan
writes them. A haplotype CN may be:
  int                 a normal call (coverage > 0 unless CN is 0, confidence 0.95)
  None                zero-filled by Wakhan's merge step (coverage, CN, confidence all 0)
  (cn, confidence)    a call with the given confidence
"""

COLUMNS = ["chr", "start", "end", "hp1_coverage", "hp1_copynumber_state", "hp1_confidence",
           "hp2_coverage", "hp2_copynumber_state", "hp2_confidence", "svs_breakpoints_ids"]

HEADER = """#chr: chromosome number
#start: start address for CN segment
#end: end address for CN segment
#hp1_coverage / hp2_coverage: median coverage for this segment, per haplotype
#hp1_copynumber_state / hp2_copynumber_state: detected copy number state (integer/fraction), per haplotype
#hp1_confidence / hp2_confidence: confidence score, per haplotype
#svs_breakpoints_ids: corresponding structural variations (breakpoints) IDs from VCF file (union across haplotypes)
"""


def _hp(value):
    if value is None:
        return {"coverage": 0, "copynumber_state": 0, "confidence": 0}
    cn, conf = value if isinstance(value, tuple) else (value, 0.95)
    return {"coverage": round(cn * 8.0 + 0.9, 2), "copynumber_state": float(cn), "confidence": conf}


def write_profile(path, rows, columns=COLUMNS):
    with open(path, "w") as f:
        f.write(HEADER)
        f.write("#" + "\t".join(columns) + "\n")
        for chrom, start, end, hp1, hp2 in rows:
            v = {"chr": chrom, "start": start, "end": end, "svs_breakpoints_ids": ""}
            for hp, val in (("hp1", hp1), ("hp2", hp2)):
                v.update({f"{hp}_{k}": x for k, x in _hp(val).items()})
            f.write("\t".join(str(v[c]) for c in columns) + "\n")
    return path
