#!/usr/bin/env python3
"""
Build Padfoot's compact gene annotation from a GENCODE (or Ensembl-style) GFF3.

    python scripts/build_annotation.py gencode.v50.basic.annotation.gff3.gz beds/hg38_genes.tsv.gz

The bundled tables were made from the GENCODE "basic" GFF3 files listed in
beds/README.md. Padfoot --gff accepts either a GFF3 (converted on the fly with
the same code) or a table made by this script.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from padfoot.annotation import build_annotation   # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gff3", help="input GFF3 (plain or gzipped)")
    ap.add_argument("out", help="output table (.tsv or .tsv.gz)")
    args = ap.parse_args()
    counts = build_annotation(args.gff3, args.out)
    print(f"{args.out}: {counts['genes']} genes, {counts['transcripts']} transcripts")


if __name__ == "__main__":
    main()
