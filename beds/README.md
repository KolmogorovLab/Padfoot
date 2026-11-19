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

## RepeatMasker (`<genome>_rm.bed.gz`)

Unchanged from Padfoot v0.1. *TODO: record the source and version.*
