"""
Synthetic gene models used by the tests.

Coordinates are 1-based, inclusive (as in GFF3). Each gene has one transcript
for now; the v0.2 annotation tests will add alternative transcripts.

Layout (chr1 unless stated):

  GENEA  +  100,000-128,999   13 exons of 200 bp, every 2,400 bp
                               exon k = [100000 + (k-1)*2400, +199]
  GENEB  -  200,000-219,999    3 exons (exon1 is the right-most one)
  GENEC  +  300,000-309,999    2 exons   } padded intervals (+-10 kb) overlap
  GENED  +  315,000-324,999    3 exons   }
  GENEG  +  chr2:100,000-109,999   3 exons (fusion partner on another chromosome)
"""

import gzip

PADDING = 10_000


def _genea_exons():
    return [(100_000 + k * 2400, 100_000 + k * 2400 + 199) for k in range(13)]


# gene -> (chrom, strand, [exons in genomic order], five_prime_UTR, three_prime_UTR)
GENES = {
    "GENEA": ("chr1", "+", _genea_exons(), (100_000, 100_099), (128_900, 128_999)),
    "GENEB": ("chr1", "-", [(200_000, 200_499), (210_000, 210_199), (219_500, 219_999)],
              (219_900, 219_999), (200_000, 200_099)),
    "GENEC": ("chr1", "+", [(300_000, 300_199), (309_800, 309_999)],
              (300_000, 300_049), (309_950, 309_999)),
    "GENED": ("chr1", "+", [(315_000, 315_199), (320_000, 320_199), (324_800, 324_999)],
              (315_000, 315_049), (324_950, 324_999)),
    "GENEG": ("chr2", "+", [(100_000, 100_199), (104_000, 104_199), (109_800, 109_999)],
              (100_000, 100_049), (109_950, 109_999)),
}


def exon(gene, k):
    """Coordinates of exon k (transcript order, 1-based) of a toy gene."""
    chrom, strand, exons, _, _ = GENES[gene]
    ordered = exons if strand == "+" else exons[::-1]
    return ordered[k - 1]


def write_v01_annotation(path):
    """Write the toy genes in the compact format bundled with Padfoot v0.1
    (beds/hg38.gff3.gz): chrom, gene, strand, type, start, end."""
    rows = []
    for gene, (chrom, strand, exons, utr5, utr3) in GENES.items():
        start, end = exons[0][0], exons[-1][1]
        rows.append((chrom, gene, strand, "gene", start, end))
        rows.append((chrom, gene, strand, f"{gene}-201", start, end))
        ordered = exons if strand == "+" else exons[::-1]
        for k, (s, e) in enumerate(ordered, 1):
            rows.append((chrom, gene, strand, f"exon{k}", s, e))
        rows.append((chrom, gene, strand, "five_prime_UTR", *utr5))
        rows.append((chrom, gene, strand, "three_prime_UTR", *utr3))
    with gzip.open(path, "wt") as f:
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    return path
