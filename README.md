# SV association architecture

Code for the analyses underlying a study of rare structural-variant (SV) gene-level associations in UK Biobank whole-genome sequencing. This is a research reproduction pipeline, not a copy of the working notebooks or a general-purpose package. It includes statistical analyses, a fixed trait mapping, and an input schema; it omits restricted association data, generated results, figure layout, and Supplementary Table formatting.

The primary analysis uses SV MAF <5% and 72 non-ratio traits. The fixed manuscript set comprises 404 significant gene–trait associations, 193 genes, and 169 exact lead SVs. Rare-SNV convergence used matched GENCODE gene-body membership: 204/363 evaluable SV associations (56.2%) had rare-SNV gene-level support, versus 4.86869/363 under the matched null (41.90-fold). These are reference checks, not bundled output.

## Environment and inputs

Use Python 3.14 with [`requirements.txt`](requirements.txt). A full run requires authorized access to UK Biobank SV summary statistics, exact phenotype-matched WGS SNV summary statistics, and the GENCODE v49 basic GRCh38 GTF. See [`examples/input_schema.md`](examples/input_schema.md) for the private directory layout and required columns. No access-controlled data or per-trait result files are distributed here. Keep generated outputs in a private directory outside this repository.

Most scripts take `--project-root <PRIVATE_ROOT>`; Steps 32 and 33 require `SVARCH_PROJECT_ROOT=<PRIVATE_ROOT>`. The rare-SNV batch script takes explicit input, output, and scratch paths. `metadata/traits.tsv` fixes trait category, module, family, and ratio status. The exact trait-to-GCST mapping is fixed in the rare-SNV batch script.

## Run order

This is the dependency order for an authorized reproduction; inspecting the code does not require running it. `--help` lists optional flags. Several later scripts perform large permutation analyses.

| Order | Code | Purpose / principal private output |
| --- | --- | --- |
| 1 | `scripts/01_prepare_reference.py --gtf <GTF> --out <PRIVATE_ROOT>/results/gene_table.parquet` | GENCODE v49 gene-body table |
| 2 | `scripts/02_03_run_sv_gene_acat.py --project-root <PRIVATE_ROOT>` | SV–gene overlap, annotations, MAF weights, SV-type-stratified gene-level ACAT; `results/<TRAIT>/{sv_weighted,acat_gene}.parquet` |
| 3 | `scripts/03_primary_associations.py --project-root <PRIVATE_ROOT>` | Significant associations and exact lead-SV recurrence; `results/sv_pleiotropy/sv_pleiotropy_master.non_ratio_primary.csv` |
| 4 | `scripts/03_qc_caution.py --project-root <PRIVATE_ROOT>` | Raw-summary QC flags for recurrent exact lead SVs; updates recurrent summary |
| 5 | `scripts/04_lead_sv_decomposition.py --project-root <PRIVATE_ROOT>` | Lead removal and driver classes; `driver_decomposition/` |
| 6 | `scripts/05_recurrence_enrichment.py --project-root <PRIVATE_ROOT>` | Trait-level and family-collapsed module/category enrichment; `category_enrichment/` |
| 7 | `scripts/06_locus_audit.py --project-root <PRIVATE_ROOT>` | Gene-target ambiguity, origin flags, exclusion sensitivity; `target_ambiguity/`, `somatic_immune_audit/` |
| 8 | `scripts/07_recurrent_sv_features.py --project-root <PRIVATE_ROOT>` | Recurrent-SV feature comparisons and joint model; `recurrent_feature_architecture/` |
| 9 | `scripts/27_full_rare_snp_acat.py --gwas-dir <GWAS_DIR> --gene-table <GENE_TABLE> --gtf <GTF> --out-dir <SNV_OUT> --scratch-dir <SCRATCH> --preflight-only` | Check rare-SNV inputs; rerun without `--preflight-only` for gene-level ACAT. Place complete output under `<PRIVATE_ROOT>/results_server/step27_full_rare_snv_gene_body_all_traits/` |
| 10 | `scripts/08_rare_snv_convergence.py --project-root <PRIVATE_ROOT>` | Gene-level SV/SNV convergence and matched null; `full_sv_rare_snv_convergence_gene_body/` |
| 11 | `scripts/09_representative_loci.py --project-root <PRIVATE_ROOT>` | Predefined complete/partial/absent representative-locus selection; `candidate_locus_dossiers_gene_body/` |
| 12 | `SVARCH_PROJECT_ROOT=<PRIVATE_ROOT> python scripts/32_maf1_sensitivity.py` | SV MAF <1% sensitivity; `maf1_sensitivity/` |
| 13 | `SVARCH_PROJECT_ROOT=<PRIVATE_ROOT> python scripts/33_lead_criticality_locality_audit.py` | Random-removal and conditional residual-SV analyses; `lead_criticality_locality_audit/` |

`src/svarch/acat.py` contains the shared ACAT implementation. Scripts write result CSV/parquet files but do not render figures or format a Supplementary Tables workbook. The code-to-claim/output map is in [`docs/provenance.md`](docs/provenance.md).

## Interpretation

Lead-SV removal measures statistical concentration, not causality. Conditional residual-SV overlap and proximity were not enriched; same-SV-type concordance was enriched. Rare-SNV support indicates gene-level convergence, not colocalization or a shared causal mechanism. Absence of support does not establish SV specificity or SNV independence. The method is **SV-type-stratified gene-level ACAT**, not ACAT-O.

This release reorganizes code without changing the fixed analysis design. The complete restricted-data pipeline has not been rerun from this distribution copy; focused equivalence and static checks are documented in the release history.
