"""MAF <1% sensitivity analysis for the rare-SV architecture results.

This script reuses the cached MAF <5% ``sv_weighted.parquet`` tables, applies
MAF <1%, and reruns the same SV-type-stratified gene-level ACAT. It then
reconstructs lead-SV dependence, exact lead-SV recurrence, individual-trait
category enrichment, and trait-family-collapsed enrichment. Existing primary
results are read only and are never overwritten.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from itertools import product
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from acat import run_acat_gene


PROJECT_ROOT = RESULTS_DIR = PRIMARY_DIR = OUT_DIR = None

MAF_CUTOFF = 0.01

RATIO_TRAITS = {
    "Monocyte_to_lymphocyte_ratio",
    "Neutrophil_to_lymphocyte_ratio",
    "Non_high_density_lipoprotein_Apolipoprotein_B_ratio",
    "Platelet_to_lymphocyte_ratio",
    "Waist_hip_ratio",
}

# Full 77-trait mapping used in Step 16. Ratio traits are excluded below.
TRAIT_CATEGORY = {
    "Platelet_count": "Platelet", "Platelet_crit": "Platelet",
    "Platelet_distribution_width": "Platelet",
    "Mean_platelet_thrombocyte_volume": "Platelet",
    "Platelet_to_lymphocyte_ratio": "Platelet",
    "Red_blood_cell_erythrocyte_count": "RBC",
    "Red_blood_cell_erythrocyte_distribution_width": "RBC",
    "Haemoglobin_concentration": "RBC", "Haematocrit_percentage": "RBC",
    "Mean_corpuscular_volume": "RBC", "Mean_corpuscular_haemoglobin": "RBC",
    "Mean_corpuscular_haemoglobin_concentration": "RBC",
    "Mean_sphered_cell_volume": "RBC",
    "Reticulocyte_count": "Reticulocyte",
    "Reticulocyte_percentage": "Reticulocyte",
    "Immature_reticulocyte_fraction": "Reticulocyte",
    "High_light_scatter_reticulocyte_count": "Reticulocyte",
    "High_light_scatter_reticulocyte_percentage": "Reticulocyte",
    "Mean_reticulocyte_volume": "Reticulocyte",
    "White_blood_cell_leukocyte_count": "WBC", "Lymphocyte_count": "WBC",
    "Lymphocyte_percentage": "WBC", "Monocyte_count": "WBC",
    "Monocyte_percentage": "WBC", "Monocyte_to_lymphocyte_ratio": "WBC",
    "Neutrophill_count": "WBC", "Neutrophill_percentage": "WBC",
    "Neutrophil_to_lymphocyte_ratio": "WBC", "Eosinophill_count": "WBC",
    "Eosinophill_percentage": "WBC", "Basophill_count": "WBC",
    "Basophill_percentage": "WBC",
    "High_density_lipoprotein": "Lipid",
    "Low_density_lipoprotein_calculated": "Lipid",
    "Low_density_lipoprotein_measured": "Lipid",
    "Cholesterol_total": "Lipid", "Triglycerides": "Lipid",
    "Non_high_density_lipoprotein": "Lipid", "Apolipoprotein_A": "Lipid",
    "Apolipoprotein_B": "Lipid", "Remnant_cholesterol": "Lipid",
    "Non_high_density_lipoprotein_Apolipoprotein_B_ratio": "Lipid",
    "Alanine_aminotransferase": "Liver",
    "Aspartate_aminotransferase": "Liver",
    "Gamma_glutamyltransferase": "Liver", "Total_bilirubin": "Liver",
    "Direct_bilirubin": "Liver", "Albumin": "Liver",
    "Total_protein": "Liver", "Alkaline_phosphatase": "Liver",
    "Creatinine": "Renal", "Cystatin_C": "Renal", "Urea": "Renal",
    "Urate": "Renal", "Phosphate": "Renal", "Calcium": "Renal",
    "Glucose": "Glycemic", "Glycated_haemoglobin_HbA1c": "Glycemic",
    "Body_mass_index": "Anthropometric", "Weight": "Anthropometric",
    "Hip_circumference": "Anthropometric",
    "Waist_circumference": "Anthropometric",
    "Waist_hip_ratio": "Anthropometric", "Seated_height": "Anthropometric",
    "Sitting_height": "Anthropometric",
    "Systolic_blood_pressure": "Cardiovascular",
    "Diastolic_blood_pressure": "Cardiovascular",
    "Mean_arterial_pressure": "Cardiovascular",
    "Pulse_pressure": "Cardiovascular", "Pulse_rate": "Cardiovascular",
    "IGF_1": "Hormonal", "Testosterone": "Hormonal",
    "Oestradiol": "Hormonal",
    "Sex_hormon_binding_globulin_SHBG": "Hormonal", "Vitamin_D": "Hormonal",
    "Rheumatoid_factor": "Inflammatory",
    "C_reactive_protein": "Inflammatory",
}


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(len(values), dtype=float)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def fisher_freeman_halton_2xk(table: np.ndarray) -> float:
    """Two-sided conditional exact test for a 2-by-k contingency table."""
    table = np.asarray(table, dtype=int)
    if table.ndim != 2 or table.shape[0] != 2:
        raise ValueError("Expected a 2-by-k contingency table")
    table = table[:, table.sum(axis=0) > 0]
    column_totals = table.sum(axis=0)
    first_row_total = int(table[0].sum())
    denominator = comb(int(table.sum()), first_row_total)

    def conditional_probability(first_row: tuple[int, ...] | np.ndarray) -> float:
        numerator = 1
        for column_total, count in zip(column_totals, first_row):
            numerator *= comb(int(column_total), int(count))
        return numerator / denominator

    observed_probability = conditional_probability(table[0])
    tolerance = max(1e-15, observed_probability * 1e-12)
    p_value = 0.0
    for first_row in product(*(range(int(total) + 1) for total in column_totals)):
        if sum(first_row) != first_row_total:
            continue
        probability = conditional_probability(first_row)
        if probability <= observed_probability + tolerance:
            p_value += probability
    return min(float(p_value), 1.0)


def semicolon_join(values: pd.Series) -> str:
    return "; ".join(sorted(pd.Series(values).dropna().astype(str).unique()))


def annotation_layer(row: pd.Series) -> str:
    if bool(row["is_cds"]) or bool(row["is_utr"]) or bool(row["is_promoter"]):
        return "Functional"
    return "Intronic"


def trait_family_map() -> dict[str, str]:
    path = Path(__file__).resolve().parents[1] / "metadata/traits.tsv"
    mapping = pd.read_csv(path, sep="\t")
    return dict(zip(mapping["trait"], mapping["family"]))


def analyze_trait(path: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    trait = path.parent.name
    sv = pd.read_parquet(path)
    sv = sv.loc[sv["maf"].lt(MAF_CUTOFF)].copy()
    if sv.empty:
        return pd.DataFrame(), {
            "trait": trait, "n_sv_gene_pairs": 0, "n_unique_svs": 0,
            "n_genes_tested": 0, "bonferroni_threshold": np.nan,
            "n_significant_associations": 0,
        }

    # Preserve the original Beta(1,25) weights already stored in the cache.
    sv["w_final"] = sv["w_maf"]
    _, gene = run_acat_gene(sv, min_sv=1)
    threshold = 0.05 / len(gene)
    significant = gene.loc[gene["p_acat_o"].lt(threshold)].copy()

    summary = {
        "trait": trait,
        "n_sv_gene_pairs": len(sv),
        "n_unique_svs": sv["sv_id"].nunique(),
        "n_genes_tested": len(gene),
        "bonferroni_threshold": threshold,
        "n_significant_associations": len(significant),
    }
    if significant.empty:
        return pd.DataFrame(), summary

    rows: list[dict[str, object]] = []
    for gene_row in significant.itertuples(index=False):
        gene_sv = sv.loc[sv["gene_id"].eq(gene_row.gene_id)].copy()
        gene_sv = gene_sv.sort_values(["pval", "maf", "sv_id"])
        gene_sv = gene_sv.drop_duplicates("sv_id", keep="first")
        lead = gene_sv.iloc[0]
        remaining = gene_sv.loc[gene_sv["sv_id"].ne(lead["sv_id"])].copy()

        residual_p = np.nan
        residual_significant = False
        residual = None
        if remaining.empty:
            driver_class = "single-SV-only"
        else:
            _, residual_gene = run_acat_gene(remaining, min_sv=1)
            if not residual_gene.empty:
                residual_p = float(residual_gene.iloc[0]["p_acat_o"])
                residual_significant = residual_p < threshold
            driver_class = (
                "multi-SV-supported" if residual_significant
                else "lead-SV-dependent"
            )
            residual = remaining.iloc[0]

        record = {
            "trait": trait,
            "trait_category": TRAIT_CATEGORY[trait],
            "gene_id": gene_row.gene_id,
            "gene_name": gene_row.gene_name,
            "p_acat_o": float(gene_row.p_acat_o),
            "sv_bonf_threshold": threshold,
            "n_strata": int(gene_row.n_strata),
            "n_sv_total": len(gene_sv),
            "lead_sv_id": lead["sv_id"],
            "lead_sv_chrom": lead["chrom"],
            "lead_sv_start": int(lead["sv_start"]),
            "lead_sv_end": int(lead["sv_end"]),
            "lead_sv_type": lead["sv_type"],
            "lead_sv_size": int(lead["svsize"]),
            "lead_sv_maf": float(lead["maf"]),
            "lead_sv_p": float(lead["pval"]),
            "lead_sv_beta": float(lead["Beta"]),
            "lead_sv_se": float(lead["SE"]),
            "lead_sv_layer": annotation_layer(lead),
            "driver_class": driver_class,
            "residual_p_acat": residual_p,
        }
        if residual is not None:
            distance = max(
                0,
                int(lead["sv_start"]) - int(residual["sv_end"]),
                int(residual["sv_start"]) - int(lead["sv_end"]),
            )
            record.update({
                "residual_sv_id": residual["sv_id"],
                "residual_sv_type": residual["sv_type"],
                "residual_distance_bp": distance,
                "residual_overlaps_lead": distance == 0,
                "residual_same_type": residual["sv_type"] == lead["sv_type"],
            })
        rows.append(record)

    return pd.DataFrame(rows), summary


def enrichment_table(
    lead_trait: pd.DataFrame,
    background: pd.DataFrame,
    unit_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    recurrent_ids = (
        lead_trait.groupby("lead_sv_id")[unit_column].nunique().loc[lambda x: x >= 2].index
    )
    recurrent = lead_trait.loc[lead_trait["lead_sv_id"].isin(recurrent_ids)].copy()
    categories = sorted(background["trait_category"].unique())
    n_background = background[unit_column].nunique()
    background_counts = background.groupby("trait_category")[unit_column].nunique()

    rows: list[dict[str, object]] = []
    for lead_sv_id, sub in recurrent.groupby("lead_sv_id"):
        sub = sub.drop_duplicates(unit_column)
        n_units = sub[unit_column].nunique()
        counts = Counter(sub["trait_category"])
        for category in categories:
            observed = counts.get(category, 0)
            category_total = int(background_counts.get(category, 0))
            p_value = float(
                hypergeom.sf(observed - 1, n_background, category_total, n_units)
            ) if observed else 1.0
            rows.append({
                "lead_sv_id": lead_sv_id,
                "unit": unit_column,
                "category": category,
                "observed_count": observed,
                "n_units": n_units,
                "background_category_count": category_total,
                "background_total": n_background,
                "hypergeom_p": p_value,
            })
    tests = pd.DataFrame(rows)
    if len(tests):
        tests["category_fdr"] = bh_adjust(tests["hypergeom_p"].to_numpy())
        tests["category_fdr_lt_0_05"] = tests["category_fdr"].le(0.05)
        assignment = (
            tests.groupby("lead_sv_id", as_index=False)
            .agg(
                n_significant_categories=("category_fdr_lt_0_05", "sum"),
                min_category_fdr=("category_fdr", "min"),
            )
        )
        assignment["is_category_enriched"] = assignment[
            "n_significant_categories"
        ].gt(0)
    else:
        assignment = pd.DataFrame()
    return tests, assignment


def recurrent_feature_comparison(lead_summary: pd.DataFrame) -> pd.DataFrame:
    recurrent = lead_summary.loc[lead_summary["is_recurrent"]].copy()
    enriched = recurrent.loc[recurrent["trait_category_enriched"]]
    other = recurrent.loc[~recurrent["trait_category_enriched"]]
    rows: list[dict[str, object]] = []

    for feature in ["lead_sv_size", "lead_sv_maf"]:
        statistic, p_value = mannwhitneyu(
            enriched[feature], other[feature], alternative="two-sided"
        )
        rows.append({
            "feature": feature,
            "test": "two-sided Mann-Whitney U",
            "n_category_enriched": len(enriched),
            "n_other_recurrent": len(other),
            "category_enriched_summary": float(enriched[feature].median()),
            "other_recurrent_summary": float(other[feature].median()),
            "summary_type": "median",
            "statistic": float(statistic),
            "p_value": float(p_value),
        })

    for feature in ["lead_sv_layer", "lead_sv_type"]:
        table = pd.crosstab(recurrent["trait_category_enriched"], recurrent[feature])
        p_value = fisher_freeman_halton_2xk(table.to_numpy())
        rows.append({
            "feature": feature,
            "test": "two-sided Fisher-Freeman-Halton exact",
            "n_category_enriched": len(enriched),
            "n_other_recurrent": len(other),
            "category_enriched_summary": "; ".join(
                f"{key}:{value}" for key, value in enriched[feature].value_counts().sort_index().items()
            ),
            "other_recurrent_summary": "; ".join(
                f"{key}:{value}" for key, value in other[feature].value_counts().sort_index().items()
            ),
            "summary_type": "counts",
            "statistic": np.nan,
            "p_value": p_value,
        })

    result = pd.DataFrame(rows)
    result["adjusted_p"] = bh_adjust(result["p_value"].to_numpy())
    result["adjusted_p_lt_0_05"] = result["adjusted_p"].le(0.05)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()
    global PROJECT_ROOT, RESULTS_DIR, PRIMARY_DIR, OUT_DIR
    PROJECT_ROOT = args.project_root.expanduser().resolve()
    RESULTS_DIR = PROJECT_ROOT / "results"
    PRIMARY_DIR = RESULTS_DIR / "combined"
    OUT_DIR = PRIMARY_DIR / "maf1_sensitivity"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    trait_paths = sorted((RESULTS_DIR / "traits").glob("*/sv_weighted.parquet"))
    trait_paths = [p for p in trait_paths if p.parent.name not in RATIO_TRAITS]
    if len(trait_paths) != 72:
        raise ValueError(f"Expected 72 non-ratio trait tables, found {len(trait_paths)}")
    missing_categories = sorted({p.parent.name for p in trait_paths} - set(TRAIT_CATEGORY))
    if missing_categories:
        raise ValueError(f"Missing trait categories: {missing_categories}")

    association_chunks: list[pd.DataFrame] = []
    trait_rows: list[dict[str, object]] = []
    for index, path in enumerate(trait_paths, start=1):
        associations, summary = analyze_trait(path)
        trait_rows.append(summary)
        if not associations.empty:
            association_chunks.append(associations)
        print(
            f"[{index:02d}/72] {path.parent.name}: "
            f"{summary['n_significant_associations']} associations"
        )

    associations = pd.concat(association_chunks, ignore_index=True)
    trait_summary = pd.DataFrame(trait_rows)

    lead_summary = (
        associations.groupby("lead_sv_id", as_index=False)
        .agg(
            genes=("gene_name", semicolon_join),
            n_genes=("gene_name", "nunique"),
            n_traits=("trait", "nunique"),
            traits=("trait", semicolon_join),
            categories=("trait_category", semicolon_join),
            lead_sv_type=("lead_sv_type", "first"),
            lead_sv_layer=("lead_sv_layer", semicolon_join),
            lead_sv_size=("lead_sv_size", "first"),
            lead_sv_maf=("lead_sv_maf", "median"),
            min_lead_sv_p=("lead_sv_p", "min"),
            min_gene_p=("p_acat_o", "min"),
        )
    )
    lead_summary["is_recurrent"] = lead_summary["n_traits"].ge(2)

    lead_trait = associations[
        ["lead_sv_id", "trait", "trait_category"]
    ].drop_duplicates()
    hit_background = associations[
        ["trait", "trait_category"]
    ].drop_duplicates()
    category_tests, category_assignment = enrichment_table(
        lead_trait, hit_background, "trait"
    )

    family_map = trait_family_map()
    lead_family = lead_trait.copy()
    lead_family["trait_family"] = lead_family["trait"].map(family_map)
    family_background = hit_background.copy()
    family_background["trait_family"] = family_background["trait"].map(family_map)
    if lead_family["trait_family"].isna().any():
        raise ValueError("Missing trait-family assignments")
    lead_family = lead_family.drop_duplicates(
        ["lead_sv_id", "trait_family", "trait_category"]
    )
    family_background = family_background.drop_duplicates(
        ["trait_family", "trait_category"]
    )
    family_tests, family_assignment = enrichment_table(
        lead_family, family_background, "trait_family"
    )

    lead_summary = lead_summary.merge(
        category_assignment.rename(columns={
            "n_significant_categories": "n_significant_trait_categories",
            "min_category_fdr": "min_trait_category_fdr",
            "is_category_enriched": "trait_category_enriched",
        }),
        on="lead_sv_id", how="left", validate="one_to_one",
    ).merge(
        family_assignment.rename(columns={
            "n_significant_categories": "n_significant_family_categories",
            "min_category_fdr": "min_family_category_fdr",
            "is_category_enriched": "family_category_enriched",
        }),
        on="lead_sv_id", how="left", validate="one_to_one",
    )
    for column in ["trait_category_enriched", "family_category_enriched"]:
        lead_summary[column] = lead_summary[column].eq(True)

    primary_associations = pd.read_csv(
        PRIMARY_DIR / "gene_trait_associations.csv"
    )
    primary_leads = pd.read_csv(
        PRIMARY_DIR / "lead_sv_summary.csv"
    )
    primary_category = pd.read_csv(
        PRIMARY_DIR / "category_enrichment/category_assignment.csv"
    )

    primary_lead_trait = primary_associations[
        ["lead_sv_id", "trait", "trait_category"]
    ].drop_duplicates()
    primary_hit_background = primary_associations[
        ["trait", "trait_category"]
    ].drop_duplicates()
    primary_lead_family = primary_lead_trait.copy()
    primary_lead_family["trait_family"] = primary_lead_family["trait"].map(family_map)
    primary_lead_family = primary_lead_family.drop_duplicates(
        ["lead_sv_id", "trait_family", "trait_category"]
    )
    primary_family_background = primary_hit_background.copy()
    primary_family_background["trait_family"] = primary_family_background["trait"].map(
        family_map
    )
    primary_family_background = primary_family_background.drop_duplicates(
        ["trait_family", "trait_category"]
    )
    _, primary_family_assignment = enrichment_table(
        primary_lead_family, primary_family_background, "trait_family"
    )

    maf1_keys = set(map(tuple, associations[["trait", "gene_id"]].to_numpy()))
    primary_keys = set(map(tuple, primary_associations[["trait", "gene_id"]].to_numpy()))
    maf1_lead_ids = set(lead_summary["lead_sv_id"])
    primary_lead_ids = set(primary_leads["lead_sv_id"])
    maf1_recurrent_ids = set(lead_summary.loc[lead_summary["is_recurrent"], "lead_sv_id"])
    primary_recurrent_ids = set(
        primary_leads.loc[primary_leads["is_recurrent_non_ratio"].astype(bool), "lead_sv_id"]
    )
    primary_enriched_ids = set(
        primary_category.loc[
            primary_category["n_significant_categories"].gt(0), "lead_sv_id"
        ]
    )
    maf1_enriched_ids = set(
        lead_summary.loc[lead_summary["trait_category_enriched"], "lead_sv_id"]
    )
    maf1_family_ids = set(
        lead_summary.loc[lead_summary["family_category_enriched"], "lead_sv_id"]
    )
    primary_family_ids = set(
        primary_family_assignment.loc[
            primary_family_assignment["is_category_enriched"], "lead_sv_id"
        ]
    )

    feature_comparison = recurrent_feature_comparison(lead_summary)

    driver_counts = associations["driver_class"].value_counts()
    lead_centered = int(
        associations["driver_class"].isin(
            ["single-SV-only", "lead-SV-dependent"]
        ).sum()
    )
    multi = associations.loc[associations["driver_class"].eq("multi-SV-supported")]

    comparison = pd.DataFrame([
        {"metric": "gene_trait_associations", "maf_lt_5_percent": len(primary_keys),
         "maf_lt_1_percent": len(maf1_keys), "overlap": len(primary_keys & maf1_keys)},
        {"metric": "traits_with_hits", "maf_lt_5_percent": primary_associations["trait"].nunique(),
         "maf_lt_1_percent": associations["trait"].nunique(), "overlap": np.nan},
        {"metric": "unique_genes", "maf_lt_5_percent": primary_associations["gene_id"].nunique(),
         "maf_lt_1_percent": associations["gene_id"].nunique(), "overlap": len(
             set(primary_associations["gene_id"]) & set(associations["gene_id"])
         )},
        {"metric": "unique_exact_lead_svs", "maf_lt_5_percent": len(primary_lead_ids),
         "maf_lt_1_percent": len(maf1_lead_ids), "overlap": len(primary_lead_ids & maf1_lead_ids)},
        {"metric": "recurrent_exact_lead_svs", "maf_lt_5_percent": len(primary_recurrent_ids),
         "maf_lt_1_percent": len(maf1_recurrent_ids), "overlap": len(primary_recurrent_ids & maf1_recurrent_ids)},
        {"metric": "trait_category_enriched_recurrent_svs", "maf_lt_5_percent": len(primary_enriched_ids),
         "maf_lt_1_percent": len(maf1_enriched_ids), "overlap": len(primary_enriched_ids & maf1_enriched_ids)},
        {"metric": "family_category_enriched_recurrent_svs",
         "maf_lt_5_percent": len(primary_family_ids),
         "maf_lt_1_percent": len(maf1_family_ids),
         "overlap": len(primary_family_ids & maf1_family_ids)},
        {"metric": "single_sv_only", "maf_lt_5_percent": 8,
         "maf_lt_1_percent": int(driver_counts.get("single-SV-only", 0)), "overlap": np.nan},
        {"metric": "lead_sv_dependent", "maf_lt_5_percent": 339,
         "maf_lt_1_percent": int(driver_counts.get("lead-SV-dependent", 0)), "overlap": np.nan},
        {"metric": "multi_sv_supported", "maf_lt_5_percent": 57,
         "maf_lt_1_percent": int(driver_counts.get("multi-SV-supported", 0)), "overlap": np.nan},
        {"metric": "lead_centered_fraction", "maf_lt_5_percent": 347 / 404,
         "maf_lt_1_percent": lead_centered / len(associations), "overlap": np.nan},
        {"metric": "multi_supported_within_250kb", "maf_lt_5_percent": 57,
         "maf_lt_1_percent": int(multi["residual_distance_bp"].le(250_000).sum()), "overlap": np.nan},
        {"metric": "multi_supported_same_type_within_50kb", "maf_lt_5_percent": 51,
         "maf_lt_1_percent": int((multi["residual_same_type"] & multi["residual_distance_bp"].le(50_000)).sum()),
         "overlap": np.nan},
    ])

    qc = pd.DataFrame([
        {"check": "72 non-ratio trait tables", "value": len(trait_paths),
         "expected": 72, "passed": len(trait_paths) == 72},
        {"check": "all analyzed SVs have MAF <1%", "value": float(associations["lead_sv_maf"].max()),
         "expected": "<0.01", "passed": associations["lead_sv_maf"].lt(0.01).all()},
        {"check": "unique association keys", "value": len(maf1_keys),
         "expected": len(associations), "passed": len(maf1_keys) == len(associations)},
        {"check": "driver classes sum to associations", "value": int(driver_counts.sum()),
         "expected": len(associations), "passed": int(driver_counts.sum()) == len(associations)},
        {"check": "lead SVs are represented in associations", "value": len(maf1_lead_ids),
         "expected": associations["lead_sv_id"].nunique(),
         "passed": len(maf1_lead_ids) == associations["lead_sv_id"].nunique()},
        {"check": "MAF <1% associations are a subset of primary associations",
         "value": len(maf1_keys - primary_keys), "expected": 0,
         "passed": len(maf1_keys - primary_keys) == 0},
    ])

    outputs = {
        "trait_summary.csv": trait_summary,
        "gene_trait_associations.csv": associations,
        "lead_sv_summary.csv": lead_summary,
        "trait_category_enrichment.csv": category_tests,
        "family_category_enrichment.csv": family_tests,
        "recurrent_feature_comparison.csv": feature_comparison,
        "primary_comparison.csv": comparison,
        "qc.csv": qc,
    }
    for filename, table in outputs.items():
        table.to_csv(OUT_DIR / filename, index=False)

    print("\nPrimary comparison")
    print(comparison.to_string(index=False))
    print("\nQC")
    print(qc.to_string(index=False))
    if not qc["passed"].all():
        failed = qc.loc[~qc["passed"], "check"].tolist()
        raise AssertionError(f"Step 32 QC failed: {failed}")


if __name__ == "__main__":
    main()
