# Input files

Download the three source datasets listed in the [README](../README.md#required-inputs). Keep them outside this repository and arrange them as follows:

```text
<DATA_ROOT>/
├── Gene_ref/gencode.v49.basic.annotation.gtf.gz
└── SV_association/NFE/qtbig.A/<TRAIT>_adjAgeSexYobPC_InvNorm.txt.gz

<GWAS_DIR>/
└── <SOURCE>/<GCST_ID>.tsv.gz
```

`<TRAIT>` names are listed in [`metadata/traits.tsv`](../metadata/traits.tsv). The exact trait-to-`GCST_ID` mapping is `TRAIT_TO_GCST` in [`scripts/06_snv_acat.py`](../scripts/06_snv_acat.py). `<GWAS_DIR>` is passed to that script with `--gwas-dir`.

## Required columns

- SV summary statistics (`.txt.gz`): `Chrom`, `Pos`, `Name`, `effectAllele`, `effectAlleleFreq`, `pval`, `Beta`, `SE`, `N`, `info`.
- Single-variant summary statistics (`.tsv.gz`): `chromosome`, `base_pair_location`, `effect_allele`, `other_allele`, `p_value`, `effect_allele_frequency`, `n`.

`Pos` and `base_pair_location` are 1-based. SV `Name` must retain its original identifier for exact-variant recurrence. The scripts convert positions to 0-based coordinates for GENCODE gene-body membership; promoter, CDS, and UTR annotations do not extend that membership window.
