# Input files

Download the source data listed in the [README](../README.md#required-inputs). Keep them outside this repository and arrange them as follows:

```text
<DATA_ROOT>/
├── Gene_ref/gencode.v49.basic.annotation.gtf.gz
└── SV_association/NFE/qtbig.A/<TRAIT>_adjAgeSexYobPC_InvNorm.txt.gz
```

`<TRAIT>` names are listed in [`metadata/traits.tsv`](../metadata/traits.tsv).

## Required columns

- SV summary statistics (`.txt.gz`): `Chrom`, `Pos`, `Name`, `effectAllele`, `effectAlleleFreq`, `pval`, `Beta`, `SE`, `N`, `info`.

`Pos` is 1-based. SV `Name` must retain its original identifier for exact-variant recurrence. The scripts convert positions to 0-based coordinates for GENCODE gene-body membership; promoter, CDS, and UTR annotations do not extend that membership window.
