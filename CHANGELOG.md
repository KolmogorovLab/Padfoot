# Changelog

## 0.2.0

### SV classification (same classes, bugs fixed)
- Breakend gene lookup uses an interval tree over the ±10 kb padded genes; v0.1
  assigned no gene when padded genes overlapped.
- The padding on both sides is labelled `promoter/utr` (the right-hand side was
  labelled as an intron).
- Exon/intron numbers are taken from the naming transcript; fixes the
  `len(exons)` bug (exon 4–5 labelled `promoter/utr`) and truncated intron numbers.
- Frame label uses `length % 3`, with SVLEN for insertions.
- Fusions need both breakends in gene bodies; partners are reported 5′::3′.
- Known-fusion lookup (`oncogenic_fusion`) works, in either orientation.
- One primary gene per breakend; other overlapping genes are reported separately
  and all of them get the SV in the gene tables.
- Whole genes inside a DEL/DUP/INV (≤ 1 Mb) are labelled `whole_gene_*`.

### Tiers
- `HIGH` / `MODERATE` from class × breakend location × SV type, replacing
  `HIGH`, `HIGH_ONCO`, `MODERATE`/`LOW`. Cancer-gene status is a label only.
- The `Score` column is removed from both tables.

### Copy number
- Reads the Wakhan (≥ 0.5.0) `integer_profile.bed` (`--wakhan-profile`); the
  Wakhan VCF and per-haplotype BEDs are no longer supported.
- Gene CN status (`HOMDEL`, `LOSS`, `NEUTRAL`, `GAIN`, `GAIN_HIGH`) relative to
  the sample's modal CN; LOH and LOH type; CN boundary and missing-CN flags.
  Zero-filled haplotypes are treated as missing. `--sex` controls chrX LOH.
- CN no longer assigns tiers (`ONC_AMP`, `TSG_DEL` removed).

### Annotation inputs
- Gene annotation rebuilt from GENCODE v50 (hg38) and vM25 (mm10) with
  `scripts/build_annotation.py`; transcripts: MANE Select ∪ Ensembl canonical ∪ CCDS.
  `--gff` converts a user GFF3 with the same rules.
- Bundled, cleaned cancer gene list for hg38; `--cancer-genes` to replace it or
  `none` to turn it off.
- Sources and md5 sums in `beds/README.md`; md5 sums are logged.

### Outputs
- `annotated_svs.tsv` and `by_gene.tsv` keep their other v0.1 columns and gain
  new ones (see README). In `by_gene.tsv` the v0.1 `HP1`, `HP1 CN`, `HP2`,
  `HP2 CN` columns are replaced by the CN block, and rows are per gene ID.
- New `gene_sv_pairs.tsv`: one row per gene and SV, with the class and tier for that gene.

### Command line
- New: `--wakhan-profile`, `--sex`, `--cancer-genes`, `--whole-gene-max-sv`,
  `--gain-high-multiplier`, `--repeatmasker-species`.
- Removed: `--wakhan-vcf`, `--specie`. RepeatMasker on inserted sequences runs
  only with `--repeatmasker-species`.
- Checks that `--genome` matches the reference and that `bedtools` is installed.

### Other
- Microhomology reads the reference with pysam (much faster; SVs on contigs
  missing from the reference no longer shift the results of later SVs).
- `cluster_cns.py` removed.
- pytest suite (`tests/`).

## 0.1

Initial release.
