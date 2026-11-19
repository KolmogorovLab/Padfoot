# Bundled annotation files

## Gene annotation (`<genome>_genes.tsv.gz`)

Built with `scripts/build_annotation.py` from the GENCODE "basic" GFF3. The
source file name, its md5 and the build date are recorded in the first line of
each table.

| File | Source | md5 of source | Genes | Transcripts |
|---|---|---|---|---|
| `hg38_genes.tsv.gz` | [GENCODE v50](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gff3.gz) `gencode.v50.basic.annotation.gff3.gz` | `72a7c4f0a34f139d6c8b9dd7a866637a` | 19,790 | 113,613 |
| `mm10_genes.tsv.gz` | [GENCODE vM25](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gff3.gz) `gencode.vM25.basic.annotation.gff3.gz` (last mm10 release) | `bc9ff6b8c68f5c7eb2c39f5243734ddd` | 21,766 | 34,504 |

**Genes:** `gene_type=protein_coding`.

**Transcripts kept:** `protein_coding` or `protein_coding_LoF` transcripts
tagged `MANE_Select`, `Ensembl_canonical` or `CCDS`.

**Naming transcript** (exon and intron numbers), one per gene:
1. MANE Select;
2. otherwise Ensembl canonical;
3. otherwise the longest CCDS transcript;
4. otherwise the longest protein-coding transcript, which is then the only transcript kept.

GENCODE vM25 has no MANE or Ensembl canonical tags, so mouse genes are named
after their longest CCDS transcript.

To rebuild:

```
python scripts/build_annotation.py gencode.v50.basic.annotation.gff3.gz beds/hg38_genes.tsv.gz
```

## Cancer gene list (`cancer_genes.tsv`)

Human only; used by default with `--specie human`. Replace it with
`--cancer-genes <file>`. Padfoot reads `Gene_symbol`, `Role_in_cancer`
(oncogene / TSG / oncogene,TSG) and `fusion` (comma-separated partner
symbols); other columns are kept for reference. Both are labels only and
never change a tier. A fusion pair matches when it is listed under either
partner.

- **Source:** the lab's curated list, `aysegokce/PadfootFiles/cancer_genes.tsv`
  at commit 059a3a0 (md5 `429b50dc806bae0058e6acc14b1f423e`). *TODO: document
  how the list was compiled.*
- **Cleaned** on 2026-10-30 against the HGNC complete set: symbols updated to
  current HGNC names, Excel-mangled SEPTIN symbols restored, duplicate rows
  merged, non-gene rows dropped, fusion partners corrected and placeholders
  removed. The original symbol is in `original_symbol`.
- **Contents:** 1,422 genes, 974 with a role; 447 distinct fusion pairs
  (763 listings; 131 pairs are listed under one partner only).
- **Not matched by the protein-coding annotation:** non-coding entries (IG/TR
  loci, MALAT1, TERC, PVT1, pseudogenes …) and three unresolved partner
  names (CEP88, C12orf9, ATL).

md5 of the bundled file: `3adae88d7d6a04dadd36f59fc16708d8`.

## RepeatMasker (`<genome>_rm.bed.gz`)

Unchanged from Padfoot v0.1. *TODO: record the source and version.*
