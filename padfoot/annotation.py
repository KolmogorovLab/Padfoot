"""
Gene annotation for Padfoot.

build_annotation() converts a GFF3 (GENCODE, or Ensembl-style) into Padfoot's
compact annotation table; load_annotation() reads that table and answers
"which genes does this breakend hit, and where".

Compact table (tab-separated, optionally gzipped), one row per transcript:

    #padfoot_annotation  format=1  source=...  md5=...  built=...  transcripts=...
    #chrom start end strand gene_id gene transcript_id tags naming cds_start cds_end exon_starts exon_ends

Coordinates are 1-based and inclusive, as in GFF3. Exons are listed in genomic
order. `naming` is 1 for the transcript used for exon/intron numbers.

Transcript set: protein-coding (incl. protein_coding_LoF) transcripts of protein-coding genes tagged
MANE_Select, Ensembl_canonical or CCDS. One naming transcript per gene:
MANE_Select, else Ensembl_canonical, else the longest CCDS transcript, else the
longest protein-coding transcript (which is then the only transcript kept).
"""

import bisect
import datetime
import gzip
import hashlib
import logging
import os
from collections import defaultdict

from intervaltree import IntervalTree

logger = logging.getLogger()

PADDING = 10_000
FORMAT_HEADER = "#padfoot_annotation"
FORMAT_VERSION = "1"
SELECT_TAGS = ("MANE_Select", "Ensembl_canonical", "CCDS")
# transcript types kept for protein-coding genes; GENCODE uses protein_coding_LoF
# for genes whose reference allele is a known loss-of-function variant (GSTT2, CASP12, ...)
CODING_TRANSCRIPT_TYPES = ("protein_coding", "protein_coding_LoF")
COLUMNS = ["#chrom", "start", "end", "strand", "gene_id", "gene", "transcript_id", "tags", "naming",
           "cds_start", "cds_end", "exon_starts", "exon_ends"]


# --------------------------------------------------------------------------
# Building the compact table from GFF3
# --------------------------------------------------------------------------

def _open_text(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def _attrs(field):
    out = {}
    for kv in field.rstrip("\n").rstrip(";").split(";"):
        k, sep, v = kv.partition("=")
        if sep:
            out[k.strip()] = v
    return out


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def is_compact_annotation(path):
    with _open_text(path) as f:
        return f.readline().startswith(FORMAT_HEADER)


def _choose_transcripts(txs):
    """Return (selected transcripts, naming transcript) for one gene."""
    span = lambda t: t["end"] - t["start"]                      # noqa: E731
    by_tag = lambda tag: [t for t in txs if tag in t["tags"]]    # noqa: E731
    naming = None
    for tag in ("MANE_Select", "Ensembl_canonical"):
        if by_tag(tag):
            naming = by_tag(tag)[0]
            break
    if naming is None and by_tag("CCDS"):
        naming = max(by_tag("CCDS"), key=span)
    if naming is None:
        naming = max(txs, key=span)
    selected = [t for t in txs if t["tags"] & set(SELECT_TAGS)]
    if naming not in selected:
        selected.append(naming)
    return selected, naming


def build_annotation(gff3_path, out_path):
    """Convert a GFF3 to the compact annotation table. Returns summary counts."""
    # pass 1: genes and protein-coding transcripts
    genes = {}
    txs = {}
    chrom_order = []
    with _open_text(gff3_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9:
                continue
            ftype = c[2]
            if ftype == "gene":
                a = _attrs(c[8])
                gid = a.get("ID", "")
                genes[gid] = (a.get("gene_name") or a.get("Name") or gid,
                              a.get("gene_type") or a.get("biotype"),
                              a.get("gene_id") or gid.replace("gene:", ""))
            elif ftype in ("transcript", "mRNA"):
                a = _attrs(c[8])
                gname, gtype, gid = genes.get(a.get("Parent"), (None, None, a.get("Parent")))
                gname = a.get("gene_name") or gname
                gtype = a.get("gene_type") or gtype
                gid = a.get("gene_id") or gid
                ttype = a.get("transcript_type") or a.get("biotype") or ("protein_coding" if ftype == "mRNA" else None)
                if gtype != "protein_coding" or ttype not in CODING_TRANSCRIPT_TYPES:
                    continue
                tags = set(filter(None, a.get("tag", "").split(",")))
                if a.get("ccdsid"):
                    tags.add("CCDS")
                if c[0] not in chrom_order:
                    chrom_order.append(c[0])
                txs[a["ID"]] = dict(chrom=c[0], start=int(c[3]), end=int(c[4]), strand=c[6], gene_id=gid,
                                    gene=gname or gid, transcript_id=a.get("transcript_id") or a["ID"].replace("transcript:", ""),
                                    tags=tags, exons=[], cds=[])

    by_gene = defaultdict(list)
    for tid, t in txs.items():
        t["id"] = tid
        by_gene[(t["chrom"], t["gene_id"])].append(t)
    selected_ids = set()
    for gtx in by_gene.values():
        selected, naming = _choose_transcripts(gtx)
        for t in selected:
            t["naming"] = t is naming
            selected_ids.add(t["id"])

    # pass 2: exons and CDS of the selected transcripts
    with _open_text(gff3_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] not in ("exon", "CDS"):
                continue
            i = c[8].find("Parent=")
            if i < 0:
                continue
            parents = c[8][i + 7:].split(";", 1)[0].rstrip("\n").split(",")
            for p in parents:
                if p in selected_ids:
                    txs[p]["exons" if c[2] == "exon" else "cds"].append((int(c[3]), int(c[4])))

    rows = []
    for tid in selected_ids:
        t = txs[tid]
        if not t["exons"]:
            continue
        exons = sorted(t["exons"])
        cds_s = min(s for s, _ in t["cds"]) if t["cds"] else "."
        cds_e = max(e for _, e in t["cds"]) if t["cds"] else "."
        rows.append((chrom_order.index(t["chrom"]), t["start"], t["gene_id"], t["transcript_id"], [
            t["chrom"], t["start"], t["end"], t["strand"], t["gene_id"], t["gene"], t["transcript_id"],
            ",".join(sorted(t["tags"] & set(SELECT_TAGS))) or ".", int(t["naming"]), cds_s, cds_e,
            ",".join(str(s) for s, _ in exons), ",".join(str(e) for _, e in exons)]))
    rows.sort(key=lambda r: r[:4])

    meta = [FORMAT_HEADER, f"format={FORMAT_VERSION}", f"source={os.path.basename(gff3_path)}",
            f"md5={_md5(gff3_path)}", f"built={datetime.date.today().isoformat()}",
            "transcripts=" + "|".join(SELECT_TAGS) + "(protein_coding)"]
    opener = gzip.open if str(out_path).endswith(".gz") else open
    with opener(out_path, "wt") as out:
        out.write("\t".join(meta) + "\n")
        out.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            out.write("\t".join(map(str, r[4])) + "\n")
    n_genes = len({(r[4][0], r[4][4]) for r in rows})
    return dict(genes=n_genes, transcripts=len(rows))


# --------------------------------------------------------------------------
# Reading and lookup
# --------------------------------------------------------------------------

class Transcript(object):
    __slots__ = ("transcript_id", "tags", "start", "end", "cds_start", "cds_end", "exons", "exon_starts")

    def __init__(self, transcript_id, tags, start, end, cds_start, cds_end, exons):
        self.transcript_id = transcript_id
        self.tags = tags
        self.start = start
        self.end = end
        self.cds_start = cds_start
        self.cds_end = cds_end
        self.exons = exons                      # genomic order
        self.exon_starts = [s for s, _ in exons]

    def exon_index(self, pos):
        """Index (genomic order) of the exon containing pos, or None."""
        i = bisect.bisect_right(self.exon_starts, pos) - 1
        if i >= 0 and pos <= self.exons[i][1]:
            return i
        return None


class Gene(object):
    __slots__ = ("gene_id", "name", "chrom", "strand", "start", "end", "naming", "transcripts")

    def __init__(self, gene_id, name, chrom, strand):
        self.gene_id = gene_id
        self.name = name
        self.chrom = chrom
        self.strand = strand
        self.start = None
        self.end = None
        self.naming = None
        self.transcripts = []

    @property
    def key(self):
        return (self.chrom, self.gene_id)


class BPHit(object):
    """Where one breakend falls relative to one gene (its naming transcript)."""
    __slots__ = ("gene", "pos", "region", "number", "side", "coding", "alt_isoform_exon")

    def __init__(self, gene, pos, region, number=None, side=None, coding=None):
        self.gene = gene
        self.pos = pos
        self.region = region                    # 'exon' | 'intron' | 'flank'
        self.number = number                    # exon / intron number, transcript order
        self.side = side                        # 'upstream' | 'downstream' for flanks
        self.coding = coding                    # exon overlaps the CDS
        self.alt_isoform_exon = False           # exonic only in a non-naming transcript

    @property
    def in_body(self):
        return self.gene.start <= self.pos <= self.gene.end

    @property
    def label(self):
        if self.region == "flank":
            return "promoter/utr"
        return f"{self.region}{self.number}"

    @property
    def key(self):
        """Identity used to compare the two breakends of an SV (flank side included)."""
        return (self.gene.key, self.label, self.side)

    def as_tuple(self):
        return (self.gene.name, self.label, self.gene.strand)


class Annotation(object):
    def __init__(self, genes, meta, padding=PADDING):
        self.genes = genes
        self.meta = meta
        self.padding = padding
        self.trees = defaultdict(IntervalTree)
        for g in genes:
            self.trees[g.chrom].addi(g.start - padding, g.end + padding + 1, g)

    def _chrom(self, chrom):
        if chrom in self.trees:
            return chrom
        if "chr" + chrom in self.trees:
            return "chr" + chrom
        return None

    def lookup(self, chrom, pos):
        """All genes whose padded interval contains pos, as BPHits (sorted by gene start)."""
        chrom = self._chrom(chrom)
        if chrom is None:
            return []
        genes = sorted((iv.data for iv in self.trees[chrom][pos]), key=lambda g: (g.start, g.end, g.gene_id))
        return [locate(g, pos) for g in genes]

    def genes_by_chrom(self):
        out = defaultdict(list)
        for g in self.genes:
            out[g.chrom].append(g)
        for gs in out.values():
            gs.sort(key=lambda g: (g.start, g.end))
        return out


def locate(gene, pos):
    t = gene.naming
    plus = gene.strand == "+"
    if pos < t.start or pos > t.end:
        hit = BPHit(gene, pos, "flank", side="upstream" if (pos < t.start) == plus else "downstream")
    else:
        n = len(t.exons)
        i = t.exon_index(pos)
        if i is not None:
            s, e = t.exons[i]
            coding = t.cds_start is not None and s <= t.cds_end and e >= t.cds_start
            hit = BPHit(gene, pos, "exon", number=i + 1 if plus else n - i, coding=coding)
        else:
            i = bisect.bisect_right(t.exon_starts, pos) - 1      # exon to the left
            hit = BPHit(gene, pos, "intron", number=i + 1 if plus else n - i - 1)
    if hit.region != "exon":
        hit.alt_isoform_exon = any(tr.exon_index(pos) is not None for tr in gene.transcripts if tr is not t)
    return hit


def load_annotation(path, padding=PADDING):
    meta = {}
    genes = {}
    with _open_text(path) as f:
        first = f.readline()
        if not first.startswith(FORMAT_HEADER):
            raise ValueError(f"{path} is not a Padfoot annotation table (run scripts/build_annotation.py)")
        for field in first.rstrip("\n").split("\t")[1:]:
            k, _, v = field.partition("=")
            meta[k] = v
        for line in f:
            if line.startswith("#"):
                continue
            (chrom, start, end, strand, gene_id, name, tid, tags, naming,
             cds_s, cds_e, ex_s, ex_e) = line.rstrip("\n").split("\t")
            key = (chrom, gene_id)
            g = genes.get(key)
            if g is None:
                g = genes[key] = Gene(gene_id, name, chrom, strand)
            exons = list(zip(map(int, ex_s.split(",")), map(int, ex_e.split(","))))
            t = Transcript(tid, set(tags.split(",")) - {"."}, int(start), int(end),
                           None if cds_s == "." else int(cds_s), None if cds_e == "." else int(cds_e), exons)
            g.transcripts.append(t)
            if naming == "1":
                g.naming = t
    for g in genes.values():
        g.start = min(t.start for t in g.transcripts)
        g.end = max(t.end for t in g.transcripts)
        if g.naming is None:
            g.naming = g.transcripts[0]
    logger.info("Loaded annotation %s: %d genes, %d transcripts (source %s, built %s)", path, len(genes),
                sum(len(g.transcripts) for g in genes.values()), meta.get("source"), meta.get("built"))
    return Annotation(list(genes.values()), meta, padding)
