# SV association architecture

Analysis code for *Gene-level analysis of structural variants reveals a lead-variant-centered signal architecture and phenotype-structured recurrence*.

## Set up the environment

We recommend Conda. From the repository directory, create a Python 3.14 environment and install the pinned packages in [`requirements.txt`](requirements.txt):

```bash
conda create -n sv-architecture python=3.14 pip -y
conda activate sv-architecture
python -m pip install -r requirements.txt
```

Run the analysis commands below with this environment active. Use `conda deactivate` when finished.

## Inputs

Download these inputs from their original providers:

| Input | Source | Expected location |
| --- | --- | --- |
| UK Biobank WGS SV association summary statistics | [deCODE summary data](https://www.decode.com/summarydata/) — UK Biobank Whole-Genome Sequencing Consortium entry | `<DATA_ROOT>/SV_association/NFE/qtbig.A/` |
| Phenotype-matched single-variant summary statistics | [GWAS Catalog](https://www.ebi.ac.uk/gwas/downloads/summary-statistics); exact trait-to-GCST IDs in [`scripts/06_snv_acat.py`](scripts/06_snv_acat.py) | `<GWAS_DIR>/<SOURCE>/<GCST_ID>.tsv.gz` |
| GENCODE v49 basic gene annotation (GRCh38) | [GENCODE human release 49](https://www.gencodegenes.org/human/release_49.html) | `<DATA_ROOT>/Gene_ref/gencode.v49.basic.annotation.gtf.gz` |

[`metadata/traits.tsv`](metadata/traits.tsv) records the analyzed traits and phenotype groupings. See [`examples/input_schema.md`](examples/input_schema.md) for filenames and required columns. Neither source summary statistics nor generated results are redistributed here.

Most scripts take `--project-root <DATA_ROOT>`; the rare-SNV batch script takes explicit input, output, and scratch paths. Keep downloaded inputs and generated outputs outside this repository.

## Run order

This is the dependency order for a reproduction; inspecting the code does not require running it. `--help` lists optional flags. Several later scripts perform large permutation analyses.

| Order | Code | Purpose / principal local output |
| --- | --- | --- |
| 1 | `scripts/01_reference.py --gtf <GTF> --out <DATA_ROOT>/results/gene_table.parquet` | GENCODE v49 gene-body table |
| 2 | `scripts/02_sv_acat.py --project-root <DATA_ROOT>` | SV–gene overlap, annotations, MAF weights, SV-type-stratified gene-level ACAT; `results/<TRAIT>/{sv_weighted,acat_gene}.parquet` |
| 3 | `scripts/03_lead_architecture.py --project-root <DATA_ROOT>` | Significant associations, lead-removal classes, random-removal and conditional residual-SV audit; `sv_pleiotropy_master.non_ratio_primary.csv`, `driver_decomposition/`, `lead_criticality_locality_audit/` |
| 4 | `scripts/04_recurrence.py --project-root <DATA_ROOT>` | Recurrent-lead-SV QC, phenotype enrichment, target/origin audit; `category_enrichment/`, `target_ambiguity/`, `somatic_immune_audit/` |
| 5 | `scripts/05_sv_features.py --project-root <DATA_ROOT>` | Recurrent-SV feature comparisons and joint model; `recurrent_feature_architecture/` |
| 6 | `scripts/06_snv_acat.py --gwas-dir <GWAS_DIR> --gene-table <GENE_TABLE> --gtf <GTF> --out-dir <SNV_OUT> --scratch-dir <SCRATCH> --preflight-only` | Check rare-SNV inputs; rerun without `--preflight-only` for gene-level ACAT. Place complete output under `<DATA_ROOT>/results_server/step27_full_rare_snv_gene_body_all_traits/` |
| 7 | `scripts/07_snv_convergence.py --project-root <DATA_ROOT>` | Gene-level SV/SNV convergence, matched null and representative loci; `full_sv_rare_snv_convergence_gene_body/`, `candidate_locus_dossiers_gene_body/` |
| 8 | `scripts/08_maf1_sensitivity.py --project-root <DATA_ROOT>` | SV MAF <1% sensitivity; `maf1_sensitivity/` |

`src/svarch/acat.py` contains the shared ACAT implementation. Scripts write result CSV/parquet files but do not render figures or format a Supplementary Tables workbook. The code-to-claim/output map is in [`docs/provenance.md`](docs/provenance.md).
Legacy output names beginning `step27_`, `step32_`, or `step33_` remain unchanged for manuscript compatibility; they do not denote additional scripts.

## Interpretation

Lead-SV removal measures statistical concentration, not causality. Conditional residual-SV overlap and proximity were not enriched; same-SV-type concordance was enriched. Rare-SNV support indicates gene-level convergence, not colocalization or a shared causal mechanism. Absence of support does not establish SV specificity or SNV independence. The method is **SV-type-stratified gene-level ACAT**, not ACAT-O.

This release reorganizes code without changing the fixed analysis design. The complete pipeline has not been rerun from this distribution copy; focused equivalence and static checks are documented in the release history.
