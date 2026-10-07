# Manuscript provenance

All output paths are relative to the private project root (`results/sv_pleiotropy/` unless otherwise stated). No output data are included here.

| Manuscript result | Source code | Principal output |
| --- | --- | --- |
| Gene-level SV landscape; 404 associations, 193 genes, 169 exact lead SVs | `01_reference.py`, `02_sv_acat.py`, `03_lead_architecture.py` | `sv_pleiotropy_master.non_ratio_primary.csv`, `lead_sv_level_pleiotropy_summary.non_ratio_primary.csv` |
| Lead-SV removal classes (8/339/57); Fig. 1B | `03_lead_architecture.py` | `driver_decomposition/lead_sv_driver_decomposition.non_ratio_primary.csv` |
| Random-removal and conditional residual-SV null; Fig. 1C, Supplementary Table S10 source | `03_lead_architecture.py` | `lead_criticality_locality_audit/step33_criticality_summary.csv`, `step33_locality_summary.csv` |
| Exact lead-SV recurrence and phenotype structure; Figs. 1D and 2 | `03_lead_architecture.py`, `04_recurrence.py` | `candidate_recurrent_lead_svs.non_ratio_primary.csv`, `category_enrichment/recurrent_lead_sv_category_assignment.csv`, `lead_sv_redundancy_sensitivity.csv` |
| Target ambiguity, blood-derived candidates, exclusion sensitivity | `04_recurrence.py` | `target_ambiguity/`, `somatic_immune_audit/` |
| Recurrent-SV features and joint model; Fig. 3 | `05_sv_features.py` | `recurrent_feature_architecture/primary_univariate_tests.csv`, `multivariable_model_diagnostics.csv` |
| MAF <1% sensitivity; Supplementary Table S9 source | `06_maf1_sensitivity.py` | `maf1_sensitivity/step32_primary_comparison.csv` |

Figure composition and Supplementary Table workbook formatting are not included here.
