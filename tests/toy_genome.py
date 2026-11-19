"""
Synthetic gene models used by the tests, written as a GENCODE-style GFF3.

Coordinates are 1-based, inclusive (as in GFF3).

Layout (chr1 unless stated):

  GENEA  +  100,000-128,999   13 exons of 200 bp, every 2,400 bp
                               exon k = [100000 + (k-1)*2400, +199]
            GENEA-201  MANE_Select, Ensembl_canonical, CCDS  (naming)
            GENEA-202  CCDS; extra exon 106,500-106,699 in intron 3
            GENEA-203  untagged (not kept); extra exon 110,500-110,699 in intron 5
  GENEB  -  200,000-219,999    3 exons (exon1 is the right-most one)
  GENEC  +  300,000-309,999    2 exons   } padded intervals (+-10 kb) overlap
  GENED  +  315,000-324,999    3 exons   } exon1 of GENED is entirely 5' UTR
  GENEE  +  400,000-449,999    3 exons
  GENEF  -  420,000-421,999    2 exons, nested in intron 1 of GENEE
  GENEG  +  chr2:100,000-109,999   3 exons, no tags (naming falls back to the
            longest transcript); fusion partner on another chromosome
  LNC1      chr1:500,000-510,000 lncRNA (not kept)
"""

PADDING = 10_000

_A = [(100_000 + k * 2400, 100_000 + k * 2400 + 199) for k in range(13)]

MANE = "basic,Ensembl_canonical,MANE_Select,CCDS"

# gene -> (chrom, strand, gene_type, [(transcript, tags, exons in genomic order, (cds_start, cds_end))])
GENES = {
    "GENEA": ("chr1", "+", "protein_coding", [
        ("GENEA-201", MANE, _A, (100_100, 128_899)),
        ("GENEA-202", "basic,CCDS", sorted(_A + [(106_500, 106_699)]), (100_100, 128_899)),
        ("GENEA-203", "basic", sorted(_A + [(110_500, 110_699)]), (100_100, 128_899)),
    ]),
    "GENEB": ("chr1", "-", "protein_coding", [
        ("GENEB-201", MANE, [(200_000, 200_499), (210_000, 210_199), (219_500, 219_999)], (200_100, 219_899)),
    ]),
    "GENEC": ("chr1", "+", "protein_coding", [
        ("GENEC-201", MANE, [(300_000, 300_199), (309_800, 309_999)], (300_050, 309_949)),
    ]),
    "GENED": ("chr1", "+", "protein_coding", [
        ("GENED-201", MANE, [(315_000, 315_199), (320_000, 320_199), (324_800, 324_999)], (320_000, 324_949)),
    ]),
    "GENEE": ("chr1", "+", "protein_coding", [
        ("GENEE-201", MANE, [(400_000, 400_199), (425_000, 425_199), (449_800, 449_999)], (400_100, 449_899)),
    ]),
    "GENEF": ("chr1", "-", "protein_coding", [
        ("GENEF-201", MANE, [(420_000, 420_199), (421_800, 421_999)], (420_100, 421_899)),
    ]),
    "GENEG": ("chr2", "+", "protein_coding", [
        ("GENEG-201", "basic", [(100_000, 100_199), (104_000, 104_199), (109_800, 109_999)], (100_050, 109_949)),
        ("GENEG-202", "basic", [(100_000, 100_199), (104_000, 104_199)], (100_050, 104_150)),
    ]),
    "LNC1": ("chr1", "+", "lncRNA", [
        ("LNC1-201", "basic,Ensembl_canonical", [(500_000, 500_500), (509_500, 510_000)], None),
    ]),
}


def exon(gene, k, transcript=0):
    """Coordinates of exon k (transcript order, 1-based) of a toy transcript."""
    _, strand, _, txs = GENES[gene]
    exons = txs[transcript][2]
    ordered = exons if strand == "+" else exons[::-1]
    return ordered[k - 1]


def write_gff3(path):
    lines = ["##gff-version 3"]
    for gi, (gene, (chrom, strand, gtype, txs)) in enumerate(GENES.items()):
        gid = f"ENSG{gi:011d}.1"
        gstart = min(e[0][0] for _, _, e, _ in txs)
        gend = max(e[-1][1] for _, _, e, _ in txs)
        common = f"gene_id={gid};gene_type={gtype};gene_name={gene}"
        lines.append(f"{chrom}\tTOY\tgene\t{gstart}\t{gend}\t.\t{strand}\t.\tID={gid};{common}")
        for ti, (tname, tags, exons, cds) in enumerate(txs):
            tid = f"ENST{gi:06d}{ti:05d}.1"
            ttype = gtype
            tattr = f"{common};transcript_id={tid};transcript_type={ttype};transcript_name={tname};tag={tags}"
            lines.append(f"{chrom}\tTOY\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\t"
                         f"ID={tid};Parent={gid};{tattr}")
            ordered = exons if strand == "+" else exons[::-1]
            for k, (s, e) in enumerate(ordered, 1):
                lines.append(f"{chrom}\tTOY\texon\t{s}\t{e}\t.\t{strand}\t.\t"
                             f"ID=exon:{tid}:{k};Parent={tid};{tattr};exon_number={k}")
                if cds and s <= cds[1] and e >= cds[0]:
                    cs, ce = max(s, cds[0]), min(e, cds[1])
                    lines.append(f"{chrom}\tTOY\tCDS\t{cs}\t{ce}\t.\t{strand}\t0\t"
                                 f"ID=CDS:{tid};Parent={tid};{tattr};exon_number={k}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path
