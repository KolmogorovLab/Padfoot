"""
Copy number from Wakhan (>= 0.5.0) and gene-level CN status.

Input: the solution's ``integer_profile.bed``. Each row is one fragment with
HP1 and HP2 on shared boundaries; the column layout comes from the
``#chr\tstart\t...`` header line. Fragments are closed intervals
[start, end] (adjacent rows meet at end + 1), as Wakhan writes them.

Wakhan's merge step fills a fragment that only one haplotype covers with
coverage 0, CN 0 and confidence 0 for the other haplotype. Such a haplotype
value is treated as missing, not as CN 0. Confidence alone is not used:
real segments can have confidence 0 (in Wakhan < 0.5 every state above the
top centre did, e.g. the ERBB2 amplicon).

Gene CN status is measured against the sample's modal CN: the integer total
CN covering the most autosomal length.
"""

import bisect
import logging
from collections import defaultdict

logger = logging.getLogger()

HAPLOTYPE_COLUMNS = ("hp1_copynumber_state", "hp2_copynumber_state")
GAIN_HIGH_MULTIPLIER = 3
END_TOLERANCE = 100_000          # profiles stop at the last full bin of a chromosome
SEXES = ("female", "male", "unknown")


class Segment:
    __slots__ = ("chrom", "start", "end", "hp1", "hp2")

    def __init__(self, chrom, start, end, hp1, hp2):
        self.chrom, self.start, self.end = chrom, start, end
        self.hp1, self.hp2 = hp1, hp2

    @property
    def total(self):
        if self.hp1 is None or self.hp2 is None:
            return None
        return self.hp1 + self.hp2

    @property
    def state(self):
        return (self.hp1, self.hp2)


class CNProfile:
    def __init__(self, segments, path=""):
        self.path = path
        self.by_chrom = defaultdict(list)
        for s in sorted(segments, key=lambda s: (s.chrom, s.start)):
            self.by_chrom[s.chrom].append(s)
        self._starts = {c: [s.start for s in segs] for c, segs in self.by_chrom.items()}
        self.modal_cn = modal_cn(self)

    def __len__(self):
        return sum(len(v) for v in self.by_chrom.values())

    def _chrom(self, chrom):
        if chrom in self.by_chrom:
            return chrom
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        return alt if alt in self.by_chrom else None

    def overlapping(self, chrom, start, end):
        """Segments overlapping the closed interval [start, end]."""
        c = self._chrom(chrom)
        if c is None:
            return []
        segs, starts = self.by_chrom[c], self._starts[c]
        i = max(bisect.bisect_right(starts, start) - 1, 0)
        out = []
        while i < len(segs) and segs[i].start <= end:
            if segs[i].end >= start:
                out.append(segs[i])
            i += 1
        return out


def _header_columns(path):
    cols = None
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("#chr\t"):
                cols = line[1:].rstrip("\n").split("\t")
    if cols is None:
        raise ValueError(f"{path}: no '#chr<TAB>start<TAB>...' header line")
    missing = [c for c in ("chr", "start", "end") + HAPLOTYPE_COLUMNS
               + ("hp1_coverage", "hp1_confidence", "hp2_coverage", "hp2_confidence") if c not in cols]
    if missing:
        raise ValueError(f"{path}: not a Wakhan >= 0.5.0 phased integer profile; missing columns {missing}")
    return cols


def _haplotype_cn(row, col, hp):
    cov, cn, conf = (float(row[col[f"hp{hp}_{k}"]]) for k in ("coverage", "copynumber_state", "confidence"))
    if cov == 0 and cn == 0 and conf == 0:
        return None                                    # zero-filled by Wakhan's merge step
    return int(round(cn))


def read_wakhan_profile(path):
    """Read a Wakhan integer_profile.bed into a CNProfile."""
    names = _header_columns(path)
    col = {c: i for i, c in enumerate(names)}
    segments = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            segments.append(Segment(row[col["chr"]], int(row[col["start"]]), int(row[col["end"]]),
                                    _haplotype_cn(row, col, 1), _haplotype_cn(row, col, 2)))
    return CNProfile(segments, path)


def is_autosome(chrom):
    return chrom.removeprefix("chr").isdigit()


def is_chrom(chrom, name):
    return chrom.removeprefix("chr") == name


def modal_cn(profile):
    """Integer total CN covering the most autosomal length (missing fragments ignored)."""
    length = defaultdict(int)
    for chrom, segs in profile.by_chrom.items():
        if not is_autosome(chrom):
            continue
        for s in segs:
            if s.total is not None:
                length[s.total] += s.end - s.start + 1
    if not length:
        return None
    return max(length, key=lambda cn: (length[cn], -cn))


def coverage_gaps(profile, chrom_lengths, chroms=None, end_tolerance=END_TOLERANCE):
    """Regions without a CN call for both haplotypes.

    Returns (absent, gaps): chromosomes in ``chroms`` with no profile at all,
    and (chrom, start, end) gaps inside the profiled chromosomes, including
    zero-filled fragments. Gaps at a chromosome end shorter than
    ``end_tolerance`` are ignored."""
    chroms = list(chroms) if chroms is not None else list(profile.by_chrom)
    absent, gaps = [], []
    for chrom in chroms:
        c = profile._chrom(chrom)
        if c is None:
            absent.append(chrom)
            continue
        length = chrom_lengths.get(c, chrom_lengths.get(chrom))
        prev_end = 0
        for s in profile.by_chrom[c]:
            if s.start > prev_end + 1 and not (prev_end == 0 and s.start <= end_tolerance):
                gaps.append((c, prev_end + 1, s.start - 1))
            if s.total is None:
                gaps.append((c, s.start, s.end))
            prev_end = max(prev_end, s.end)
        if length and length - prev_end > end_tolerance:
            gaps.append((c, prev_end + 1, length))
    return absent, _merge(gaps)


def _merge(gaps):
    out = []
    for g in sorted(gaps):
        if out and out[-1][0] == g[0] and g[1] <= out[-1][2] + 1:
            out[-1] = (g[0], out[-1][1], max(out[-1][2], g[2]))
        else:
            out.append(g)
    return out


class GeneCN:
    """CN fields of one gene body (see the docstring of gene_cn)."""
    __slots__ = ("hp1", "hp2", "total", "modal", "status", "loh", "loh_type", "boundary", "missing")

    def __init__(self):
        self.hp1 = self.hp2 = self.total = (None, None)            # (min, max)
        self.modal = None
        self.status = self.loh_type = ""
        self.loh = self.boundary = self.missing = False


GENE_CN_COLUMNS = ["HP1_CN_min", "HP1_CN_max", "HP2_CN_min", "HP2_CN_max", "total_CN_min", "total_CN_max",
                   "modal_CN", "cn_status", "LOH", "LOH_type", "cn_boundary_in_gene", "cn_missing"]


def gene_cn_fields(cn):
    """Values for GENE_CN_COLUMNS ('' when there is no CN)."""
    if cn is None:
        return [""] * len(GENE_CN_COLUMNS)
    fmt = lambda v: "" if v is None else str(v)
    vals = [fmt(v) for pair in (cn.hp1, cn.hp2, cn.total) for v in pair]
    return vals + [fmt(cn.modal), cn.status, _flag(cn.loh), cn.loh_type, _flag(cn.boundary), _flag(cn.missing)]


def _flag(b):
    return "yes" if b else "no"


def _min_max(values):
    values = [v for v in values if v is not None]
    return (min(values), max(values)) if values else (None, None)


def gene_cn(profile, chrom, start, end, sex="unknown", gain_high_multiplier=GAIN_HIGH_MULTIPLIER):
    """CN of the gene body [start, end].

    * HP1/HP2/total CN: min and max over the fragments overlapping the gene.
    * cn_status against the modal CN M: HOMDEL (min total 0), GAIN_HIGH
      (max total > multiplier x M), LOSS (0 < min total < M), GAIN (max total
      > M), else NEUTRAL; checked in this order, so losses use the minimum
      and gains the maximum when a CN breakpoint lies in the gene.
    * LOH: a fragment with one haplotype at 0 and total > 0; LOH_type
      (copy_neutral / with_loss / with_gain) from the lowest total CN among
      the LOH fragments. Not called on chrY, nor on chrX of male lines.
    * cn_boundary_in_gene: the gene overlaps fragments with different CN.
    * cn_missing: part of the gene has no CN call for both haplotypes."""
    res = GeneCN()
    res.modal = profile.modal_cn
    segs = profile.overlapping(chrom, start, end)
    covered = sum(min(s.end, end) - max(s.start, start) + 1 for s in segs if s.total is not None)
    res.missing = covered < end - start + 1
    res.hp1 = _min_max(s.hp1 for s in segs)
    res.hp2 = _min_max(s.hp2 for s in segs)
    called = [s for s in segs if s.total is not None]
    if not called:
        return res
    res.total = _min_max(s.total for s in called)
    res.boundary = len({s.state for s in called}) > 1

    lo, hi = res.total
    modal = res.modal
    if lo == 0:
        res.status = "HOMDEL"
    elif modal is not None and hi > gain_high_multiplier * modal:
        res.status = "GAIN_HIGH"
    elif modal is not None and lo < modal:
        res.status = "LOSS"
    elif modal is not None and hi > modal:
        res.status = "GAIN"
    else:
        res.status = "NEUTRAL"

    no_loh = is_chrom(chrom, "Y") or (is_chrom(chrom, "X") and sex == "male")
    loh = [s for s in called if not no_loh and s.total > 0 and 0 in (s.hp1, s.hp2)]
    if loh:
        res.loh = True
        t = min(s.total for s in loh)
        if modal is not None:
            res.loh_type = "copy_neutral" if t == modal else "with_loss" if t < modal else "with_gain"
    return res
