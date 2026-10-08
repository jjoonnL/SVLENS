"""Recurrent exact lead-SV feature comparisons and joint model."""

from __future__ import annotations

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()

    # Inputs and output location
    import numpy as np
    import pandas as pd

    PROJECT_ROOT = args.project_root.expanduser().resolve()

    RESULT_DIR = PROJECT_ROOT / "results" / "main"
    OUT_DIR = PROJECT_ROOT / "results" / "work" / "sv_features"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    lead_path = RESULT_DIR / "lead_sv_summary.csv"
    if not lead_path.exists():
        raise FileNotFoundError(lead_path)

    # Analysis block 3
    lead = pd.read_csv(lead_path)
    recurrent = lead.loc[lead["is_recurrent_non_ratio"].astype(bool)].copy()
    category_assignment = recurrent
    origin_audit = lead

    input_summary = pd.DataFrame(
        {
            "table": ["lead", "recurrent", "category_assignment", "origin_audit"],
            "rows": [len(lead), len(recurrent), len(category_assignment), len(origin_audit)],
            "unique_lead_sv_ids": [
                lead["lead_sv_id"].nunique(),
                recurrent["lead_sv_id"].nunique(),
                category_assignment["lead_sv_id"].nunique(),
                origin_audit["lead_sv_id"].nunique(),
            ],
        }
    )

    # Analysis block 5
    for table_name, table in {
        "lead": lead,
        "recurrent": recurrent,
        "category_assignment": category_assignment,
        "origin_audit": origin_audit,
    }.items():
        if table["lead_sv_id"].duplicated().any():
            raise ValueError(f"{table_name} contains duplicated lead_sv_id values")

    recurrent_ids = set(recurrent["lead_sv_id"])
    assignment_ids = set(category_assignment["lead_sv_id"])

    if len(recurrent_ids) != 64:
        raise ValueError(f"Expected 64 recurrent lead SVs, found {len(recurrent_ids)}")
    if recurrent_ids != assignment_ids:
        only_recurrent = sorted(recurrent_ids - assignment_ids)
        only_assignment = sorted(assignment_ids - recurrent_ids)
        raise ValueError(
            "Recurrent/category-assignment ID mismatch: "
            f"only recurrent={only_recurrent}; only assignment={only_assignment}"
        )

    source_compare = recurrent[[
        "lead_sv_id", "lead_sv_type", "lead_sv_layer", "lead_sv_size", "lead_sv_maf"
    ]].merge(
        category_assignment[[
            "lead_sv_id", "lead_sv_type", "lead_sv_layer", "lead_sv_size", "lead_sv_maf"
        ]],
        on="lead_sv_id",
        how="inner",
        validate="one_to_one",
        suffixes=("_recurrent", "_assignment"),
    )

    for variable in ["lead_sv_type", "lead_sv_layer"]:
        left = source_compare[f"{variable}_recurrent"].astype(str)
        right = source_compare[f"{variable}_assignment"].astype(str)
        if not left.equals(right):
            raise ValueError(f"Source mismatch for {variable}")

    for variable in ["lead_sv_size", "lead_sv_maf"]:
        left = pd.to_numeric(source_compare[f"{variable}_recurrent"], errors="coerce")
        right = pd.to_numeric(source_compare[f"{variable}_assignment"], errors="coerce")
        if not np.allclose(left, right, rtol=1e-12, atol=0, equal_nan=True):
            raise ValueError(f"Source mismatch for {variable}")

    print("Source-table QC passed: 64 recurrent lead SVs with concordant annotations.")

    # Analysis block 7
    lead_columns = [
        "lead_sv_id",
        "genes",
        "n_genes",
        "n_traits_non_ratio",
        "traits_non_ratio",
        "fine_categories_non_ratio",
        "broad_modules_non_ratio",
        "lead_sv_chrom",
        "lead_sv_pos",
        "lead_sv_type",
        "lead_sv_layer",
        "lead_sv_size",
        "lead_sv_maf",
    ]
    assignment_columns = [
        "lead_sv_id",
        "n_significant_categories",
        "significant_categories",
        "significant_category_counts",
        "plot_group",
        "min_category_fdr",
        "best_category_by_fdr",
        "observed_categories",
    ]
    origin_columns = [
        "lead_sv_id",
        "within_immune_locus",
        "known_somatic_label",
        "origin_class",
        "blood_derived_candidate",
    ]

    master = (
        lead.loc[lead["lead_sv_id"].isin(recurrent_ids), lead_columns]
        .merge(
            category_assignment[assignment_columns],
            on="lead_sv_id",
            how="left",
            validate="one_to_one",
        )
        .merge(
            origin_audit[origin_columns],
            on="lead_sv_id",
            how="left",
            validate="one_to_one",
        )
        .rename(
            columns={
                "n_traits_non_ratio": "n_traits",
                "traits_non_ratio": "traits",
                "fine_categories_non_ratio": "fine_categories",
                "broad_modules_non_ratio": "broad_modules",
                "lead_sv_chrom": "sv_chrom",
                "lead_sv_pos": "sv_pos",
                "lead_sv_type": "sv_type",
                "lead_sv_layer": "annotation_layer",
                "lead_sv_size": "sv_size_bp",
                "lead_sv_maf": "maf",
            }
        )
    )

    master["is_category_enriched"] = master["n_significant_categories"].gt(0)
    master["analysis_group"] = np.where(
        master["is_category_enriched"],
        "Category-enriched",
        "Non-enriched",
    )
    master["has_functional_annotation"] = master["annotation_layer"].str.contains(
        "Functional", regex=False, na=False
    )
    master["has_intronic_annotation"] = master["annotation_layer"].str.contains(
        "Intronic", regex=False, na=False
    )
    master["origin_audit_flag"] = master["blood_derived_candidate"].astype(bool)
    master["log10_sv_size_bp"] = np.log10(master["sv_size_bp"])
    master["log10_maf"] = np.log10(master["maf"])

    master = master.sort_values(
        ["is_category_enriched", "n_traits", "lead_sv_id"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


    # Analysis block 9
    primary_features = ["sv_size_bp", "maf", "annotation_layer", "sv_type"]
    group_counts = master["analysis_group"].value_counts().to_dict()
    allowed_sv_types = {"DEL", "DUP", "INS"}
    allowed_layers = {"Functional", "Functional; Intronic", "Intronic"}

    qc_rows = [
        {"check": "rows", "observed": len(master), "expected": 64, "passed": len(master) == 64},
        {
            "check": "unique lead_sv_id",
            "observed": master["lead_sv_id"].nunique(),
            "expected": 64,
            "passed": master["lead_sv_id"].nunique() == 64,
        },
        {
            "check": "category-enriched count",
            "observed": group_counts.get("Category-enriched", 0),
            "expected": 22,
            "passed": group_counts.get("Category-enriched", 0) == 22,
        },
        {
            "check": "non-enriched count",
            "observed": group_counts.get("Non-enriched", 0),
            "expected": 42,
            "passed": group_counts.get("Non-enriched", 0) == 42,
        },
        {
            "check": "missing primary features",
            "observed": int(master[primary_features].isna().sum().sum()),
            "expected": 0,
            "passed": master[primary_features].notna().all().all(),
        },
        {
            "check": "non-positive SV size",
            "observed": int(master["sv_size_bp"].le(0).sum()),
            "expected": 0,
            "passed": master["sv_size_bp"].gt(0).all(),
        },
        {
            "check": "non-positive MAF",
            "observed": int(master["maf"].le(0).sum()),
            "expected": 0,
            "passed": master["maf"].gt(0).all(),
        },
        {
            "check": "unexpected SV types",
            "observed": "; ".join(sorted(set(master["sv_type"]) - allowed_sv_types)),
            "expected": "",
            "passed": set(master["sv_type"]).issubset(allowed_sv_types),
        },
        {
            "check": "unexpected annotation layers",
            "observed": "; ".join(sorted(set(master["annotation_layer"]) - allowed_layers)),
            "expected": "",
            "passed": set(master["annotation_layer"]).issubset(allowed_layers),
        },
        {
            "check": "origin-audit flags",
            "observed": int(master["origin_audit_flag"].sum()),
            "expected": 3,
            "passed": int(master["origin_audit_flag"].sum()) == 3,
        },
        {
            "check": "flagged SVs are category-enriched",
            "observed": int(master.loc[master["origin_audit_flag"], "is_category_enriched"].sum()),
            "expected": 3,
            "passed": master.loc[master["origin_audit_flag"], "is_category_enriched"].all(),
        },
    ]
    qc = pd.DataFrame(qc_rows)

    if not qc["passed"].all():
        failed = qc.loc[~qc["passed"], "check"].tolist()
        raise AssertionError(f"Master-table QC failed: {failed}")

    print("Master-table QC passed.")

    # Analysis block 11
    descriptive_summary = (
        master.groupby("analysis_group", sort=False)
        .agg(
            n=("lead_sv_id", "size"),
            functional_n=("has_functional_annotation", "sum"),
            intronic_n=("has_intronic_annotation", "sum"),
            median_size_bp=("sv_size_bp", "median"),
            median_maf=("maf", "median"),
            origin_audit_flag_n=("origin_audit_flag", "sum"),
        )
        .reset_index()
    )

    # Analysis block 13
    master_path = OUT_DIR / "recurrent_lead_sv_feature_master.csv"
    master.to_csv(master_path, index=False)
    print(f"Saved: {master_path}")

    # Analysis block 15
    from itertools import product
    from math import comb

    from scipy.stats import chi2_contingency, mannwhitneyu

    GROUP_ORDER = ["Category-enriched", "Non-enriched"]
    N_BOOTSTRAP = 20_000
    BOOTSTRAP_SEED = 20260714


    def bh_adjust(p_values):
        """Return Benjamini-Hochberg adjusted P values in the original order."""
        p_values = np.asarray(p_values, dtype=float)
        if np.isnan(p_values).any():
            raise ValueError("BH adjustment requires non-missing P values")

        order = np.argsort(p_values)
        ranked = p_values[order]
        ranks = np.arange(1, len(ranked) + 1)
        adjusted_ranked = ranked * len(ranked) / ranks
        adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
        adjusted_ranked = np.clip(adjusted_ranked, 0, 1)

        adjusted = np.empty_like(adjusted_ranked)
        adjusted[order] = adjusted_ranked
        return adjusted


    def rank_biserial_from_u(u_statistic, n_enriched, n_non_enriched):
        """Positive values indicate higher ranks in the category-enriched group."""
        return 2 * u_statistic / (n_enriched * n_non_enriched) - 1


    def bootstrap_continuous_effects(enriched, non_enriched, n_bootstrap, seed):
        enriched = np.asarray(enriched, dtype=float)
        non_enriched = np.asarray(non_enriched, dtype=float)
        rng = np.random.default_rng(seed)

        rank_biserial = np.empty(n_bootstrap, dtype=float)
        median_ratio = np.empty(n_bootstrap, dtype=float)

        for index in range(n_bootstrap):
            enriched_boot = rng.choice(enriched, size=len(enriched), replace=True)
            non_enriched_boot = rng.choice(
                non_enriched, size=len(non_enriched), replace=True
            )
            u_statistic = mannwhitneyu(
                enriched_boot, non_enriched_boot, alternative="two-sided", method="auto"
            ).statistic
            rank_biserial[index] = rank_biserial_from_u(
                u_statistic, len(enriched_boot), len(non_enriched_boot)
            )
            median_ratio[index] = np.median(enriched_boot) / np.median(non_enriched_boot)

        return {
            "rank_biserial_ci_low": np.quantile(rank_biserial, 0.025),
            "rank_biserial_ci_high": np.quantile(rank_biserial, 0.975),
            "median_ratio_ci_low": np.quantile(median_ratio, 0.025),
            "median_ratio_ci_high": np.quantile(median_ratio, 0.975),
        }


    def fisher_freeman_halton_2xk(table):
        """Two-sided conditional exact test for a 2-by-k contingency table."""
        table = np.asarray(table, dtype=int)
        if table.ndim != 2 or table.shape[0] != 2:
            raise ValueError("Expected a 2-by-k contingency table")
        if (table < 0).any() or not np.equal(table, table.astype(int)).all():
            raise ValueError("Contingency-table counts must be non-negative integers")

        column_totals = table.sum(axis=0)
        if (column_totals == 0).any():
            raise ValueError("Every category must be observed at least once")

        first_row_total = int(table[0].sum())
        grand_total = int(table.sum())
        denominator = comb(grand_total, first_row_total)

        def conditional_probability(first_row):
            numerator = 1
            for column_total, count in zip(column_totals, first_row):
                numerator *= comb(int(column_total), int(count))
            return numerator / denominator

        observed_probability = conditional_probability(table[0])
        possible_probabilities = []
        ranges = [range(int(column_total) + 1) for column_total in column_totals]
        for first_row in product(*ranges):
            if sum(first_row) == first_row_total:
                possible_probabilities.append(conditional_probability(first_row))

        tolerance = max(1e-15, observed_probability * 1e-12)
        p_value = sum(
            probability
            for probability in possible_probabilities
            if probability <= observed_probability + tolerance
        )
        return min(float(p_value), 1.0)


    def bias_corrected_cramers_v(table):
        """Small-sample bias-corrected Cramer's V (Bergsma correction)."""
        table = np.asarray(table, dtype=int)
        table = table[:, table.sum(axis=0) > 0]
        if min(table.shape) < 2:
            return 0.0

        n_total = table.sum()
        n_rows, n_columns = table.shape
        phi_squared = chi2_contingency(table, correction=False).statistic / n_total
        phi_squared_corrected = max(
            0.0,
            phi_squared - ((n_columns - 1) * (n_rows - 1)) / (n_total - 1),
        )
        rows_corrected = n_rows - (n_rows - 1) ** 2 / (n_total - 1)
        columns_corrected = n_columns - (n_columns - 1) ** 2 / (n_total - 1)
        denominator = min(rows_corrected - 1, columns_corrected - 1)
        if denominator <= 0:
            return 0.0
        return float(np.sqrt(phi_squared_corrected / denominator))


    def bootstrap_bias_corrected_cramers_v(
        enriched, non_enriched, levels, n_bootstrap, seed
    ):
        level_to_index = {level: index for index, level in enumerate(levels)}
        enriched_index = np.array([level_to_index[value] for value in enriched], dtype=int)
        non_enriched_index = np.array(
            [level_to_index[value] for value in non_enriched], dtype=int
        )
        rng = np.random.default_rng(seed)
        bootstrapped_v = np.empty(n_bootstrap, dtype=float)

        for index in range(n_bootstrap):
            enriched_boot = rng.choice(
                enriched_index, size=len(enriched_index), replace=True
            )
            non_enriched_boot = rng.choice(
                non_enriched_index, size=len(non_enriched_index), replace=True
            )
            table = np.vstack(
                [
                    np.bincount(enriched_boot, minlength=len(levels)),
                    np.bincount(non_enriched_boot, minlength=len(levels)),
                ]
            )
            bootstrapped_v[index] = bias_corrected_cramers_v(table)

        return (
            np.quantile(bootstrapped_v, 0.025),
            np.quantile(bootstrapped_v, 0.975),
        )

    # Analysis block 17
    continuous_specs = {
        "sv_size_bp": "SV size",
        "maf": "MAF",
    }

    continuous_summary_rows = []
    primary_test_rows = []

    for feature_index, (feature, feature_label) in enumerate(continuous_specs.items()):
        enriched = master.loc[master["is_category_enriched"], feature].to_numpy(dtype=float)
        non_enriched = master.loc[~master["is_category_enriched"], feature].to_numpy(
            dtype=float
        )

        for group_name, values in [
            ("Category-enriched", enriched),
            ("Non-enriched", non_enriched),
        ]:
            continuous_summary_rows.append(
                {
                    "feature": feature,
                    "feature_label": feature_label,
                    "analysis_group": group_name,
                    "n": len(values),
                    "median": np.median(values),
                    "q25": np.quantile(values, 0.25),
                    "q75": np.quantile(values, 0.75),
                    "minimum": np.min(values),
                    "maximum": np.max(values),
                }
            )

        mann_whitney = mannwhitneyu(
            enriched, non_enriched, alternative="two-sided", method="auto"
        )
        rank_biserial = rank_biserial_from_u(
            mann_whitney.statistic, len(enriched), len(non_enriched)
        )
        bootstrap = bootstrap_continuous_effects(
            enriched,
            non_enriched,
            n_bootstrap=N_BOOTSTRAP,
            seed=BOOTSTRAP_SEED + feature_index,
        )

        primary_test_rows.append(
            {
                "feature": feature,
                "feature_label": feature_label,
                "feature_type": "continuous",
                "test": "Mann-Whitney U, two-sided",
                "n_category_enriched": len(enriched),
                "n_non_enriched": len(non_enriched),
                "test_statistic": mann_whitney.statistic,
                "effect_name": "rank-biserial correlation",
                "effect": rank_biserial,
                "effect_ci_low": bootstrap["rank_biserial_ci_low"],
                "effect_ci_high": bootstrap["rank_biserial_ci_high"],
                "median_ratio": np.median(enriched) / np.median(non_enriched),
                "median_ratio_ci_low": bootstrap["median_ratio_ci_low"],
                "median_ratio_ci_high": bootstrap["median_ratio_ci_high"],
                "p_value": mann_whitney.pvalue,
            }
        )

    continuous_summary = pd.DataFrame(continuous_summary_rows)

    # Analysis block 19
    categorical_specs = {
        "annotation_layer": {
            "label": "Annotation layer",
            "levels": ["Functional", "Functional; Intronic", "Intronic"],
        },
        "sv_type": {
            "label": "SV type",
            "levels": ["DEL", "DUP", "INS"],
        },
    }

    categorical_composition_rows = []

    for feature_index, (feature, spec) in enumerate(categorical_specs.items()):
        levels = spec["levels"]
        table = (
            pd.crosstab(master["analysis_group"], master[feature])
            .reindex(index=GROUP_ORDER, columns=levels, fill_value=0)
            .astype(int)
        )

        for group_name in GROUP_ORDER:
            group_total = int(table.loc[group_name].sum())
            for level in levels:
                count = int(table.loc[group_name, level])
                categorical_composition_rows.append(
                    {
                        "feature": feature,
                        "feature_label": spec["label"],
                        "analysis_group": group_name,
                        "level": level,
                        "count": count,
                        "total": group_total,
                        "percent": 100 * count / group_total,
                    }
                )

        enriched_values = master.loc[master["is_category_enriched"], feature].to_numpy()
        non_enriched_values = master.loc[~master["is_category_enriched"], feature].to_numpy()
        effect_ci_low, effect_ci_high = bootstrap_bias_corrected_cramers_v(
            enriched_values,
            non_enriched_values,
            levels=levels,
            n_bootstrap=N_BOOTSTRAP,
            seed=BOOTSTRAP_SEED + 100 + feature_index,
        )

        primary_test_rows.append(
            {
                "feature": feature,
                "feature_label": spec["label"],
                "feature_type": "categorical",
                "test": "Fisher-Freeman-Halton exact, two-sided",
                "n_category_enriched": int(table.loc["Category-enriched"].sum()),
                "n_non_enriched": int(table.loc["Non-enriched"].sum()),
                "test_statistic": np.nan,
                "effect_name": "bias-corrected Cramer's V",
                "effect": bias_corrected_cramers_v(table.to_numpy()),
                "effect_ci_low": effect_ci_low,
                "effect_ci_high": effect_ci_high,
                "median_ratio": np.nan,
                "median_ratio_ci_low": np.nan,
                "median_ratio_ci_high": np.nan,
                "p_value": fisher_freeman_halton_2xk(table.to_numpy()),
            }
        )

    categorical_composition = pd.DataFrame(categorical_composition_rows)

    # Analysis block 21
    primary_tests = pd.DataFrame(primary_test_rows)
    expected_features = ["sv_size_bp", "maf", "annotation_layer", "sv_type"]

    if primary_tests["feature"].tolist() != expected_features:
        raise AssertionError(
            "Primary-test family changed: "
            f"expected {expected_features}, observed {primary_tests['feature'].tolist()}"
        )

    primary_tests["adj_p_value_bh"] = bh_adjust(primary_tests["p_value"])
    primary_tests["bh_significant_0_05"] = primary_tests["adj_p_value_bh"].lt(0.05)

    continuous_summary_path = OUT_DIR / "primary_continuous_summary.csv"
    categorical_composition_path = OUT_DIR / "primary_categorical_composition.csv"

    continuous_summary.to_csv(continuous_summary_path, index=False)
    categorical_composition.to_csv(categorical_composition_path, index=False)

    print(f"Saved: {continuous_summary_path}")
    print(f"Saved: {categorical_composition_path}")

    # Analysis block 23
    expected_excluded_genes = {"CALR", "IGLL5", "TARP"}
    excluded_sv_details = master.loc[
        master["origin_audit_flag"],
        [
            "lead_sv_id",
            "genes",
            "analysis_group",
            "sv_type",
            "annotation_layer",
            "sv_size_bp",
            "maf",
            "origin_class",
        ],
    ].copy()

    if len(excluded_sv_details) != 3:
        raise AssertionError(f"Expected 3 excluded SVs, found {len(excluded_sv_details)}")
    if set(excluded_sv_details["genes"]) != expected_excluded_genes:
        raise AssertionError(
            "Origin-exclusion set changed: "
            f"observed {sorted(excluded_sv_details['genes'])}"
        )
    if not excluded_sv_details["analysis_group"].eq("Category-enriched").all():
        raise AssertionError("All three excluded SVs must be category-enriched")

    sensitivity_master = master.loc[~master["origin_audit_flag"]].copy()
    sensitivity_group_counts = sensitivity_master["analysis_group"].value_counts().to_dict()

    if len(sensitivity_master) != 61:
        raise AssertionError(f"Expected 61 retained recurrent SVs, found {len(sensitivity_master)}")
    if sensitivity_group_counts != {"Non-enriched": 42, "Category-enriched": 19}:
        raise AssertionError(
            f"Expected the 19/42 sensitivity split, observed {sensitivity_group_counts}"
        )
    if sensitivity_master["origin_audit_flag"].any():
        raise AssertionError("Origin-audit flags remain after exclusion")

    print("Sensitivity set QC passed: 19 category-enriched vs 42 non-enriched.")

    # Analysis block 25
    sensitivity_continuous_rows = []
    sensitivity_test_rows = []

    for feature_index, (feature, feature_label) in enumerate(continuous_specs.items()):
        enriched = sensitivity_master.loc[
            sensitivity_master["is_category_enriched"], feature
        ].to_numpy(dtype=float)
        non_enriched = sensitivity_master.loc[
            ~sensitivity_master["is_category_enriched"], feature
        ].to_numpy(dtype=float)

        for group_name, values in [
            ("Category-enriched", enriched),
            ("Non-enriched", non_enriched),
        ]:
            sensitivity_continuous_rows.append(
                {
                    "feature": feature,
                    "feature_label": feature_label,
                    "analysis_group": group_name,
                    "n": len(values),
                    "median": np.median(values),
                    "q25": np.quantile(values, 0.25),
                    "q75": np.quantile(values, 0.75),
                    "minimum": np.min(values),
                    "maximum": np.max(values),
                }
            )

        mann_whitney = mannwhitneyu(
            enriched, non_enriched, alternative="two-sided", method="auto"
        )
        rank_biserial = rank_biserial_from_u(
            mann_whitney.statistic, len(enriched), len(non_enriched)
        )
        bootstrap = bootstrap_continuous_effects(
            enriched,
            non_enriched,
            n_bootstrap=N_BOOTSTRAP,
            seed=BOOTSTRAP_SEED + 1_000 + feature_index,
        )

        sensitivity_test_rows.append(
            {
                "feature": feature,
                "feature_label": feature_label,
                "feature_type": "continuous",
                "test": "Mann-Whitney U, two-sided",
                "n_category_enriched": len(enriched),
                "n_non_enriched": len(non_enriched),
                "test_statistic": mann_whitney.statistic,
                "effect_name": "rank-biserial correlation",
                "effect": rank_biserial,
                "effect_ci_low": bootstrap["rank_biserial_ci_low"],
                "effect_ci_high": bootstrap["rank_biserial_ci_high"],
                "median_ratio": np.median(enriched) / np.median(non_enriched),
                "median_ratio_ci_low": bootstrap["median_ratio_ci_low"],
                "median_ratio_ci_high": bootstrap["median_ratio_ci_high"],
                "p_value": mann_whitney.pvalue,
            }
        )

    sensitivity_continuous_summary = pd.DataFrame(sensitivity_continuous_rows)

    # Analysis block 27
    sensitivity_categorical_rows = []

    for feature_index, (feature, spec) in enumerate(categorical_specs.items()):
        levels = spec["levels"]
        table = (
            pd.crosstab(sensitivity_master["analysis_group"], sensitivity_master[feature])
            .reindex(index=GROUP_ORDER, columns=levels, fill_value=0)
            .astype(int)
        )

        for group_name in GROUP_ORDER:
            group_total = int(table.loc[group_name].sum())
            for level in levels:
                count = int(table.loc[group_name, level])
                sensitivity_categorical_rows.append(
                    {
                        "feature": feature,
                        "feature_label": spec["label"],
                        "analysis_group": group_name,
                        "level": level,
                        "count": count,
                        "total": group_total,
                        "percent": 100 * count / group_total,
                    }
                )

        enriched_values = sensitivity_master.loc[
            sensitivity_master["is_category_enriched"], feature
        ].to_numpy()
        non_enriched_values = sensitivity_master.loc[
            ~sensitivity_master["is_category_enriched"], feature
        ].to_numpy()
        effect_ci_low, effect_ci_high = bootstrap_bias_corrected_cramers_v(
            enriched_values,
            non_enriched_values,
            levels=levels,
            n_bootstrap=N_BOOTSTRAP,
            seed=BOOTSTRAP_SEED + 1_100 + feature_index,
        )
        exact_test_table = table.loc[:, table.sum(axis=0).gt(0)]

        sensitivity_test_rows.append(
            {
                "feature": feature,
                "feature_label": spec["label"],
                "feature_type": "categorical",
                "test": "Fisher-Freeman-Halton exact, two-sided",
                "n_category_enriched": int(table.loc["Category-enriched"].sum()),
                "n_non_enriched": int(table.loc["Non-enriched"].sum()),
                "test_statistic": np.nan,
                "effect_name": "bias-corrected Cramer's V",
                "effect": bias_corrected_cramers_v(table.to_numpy()),
                "effect_ci_low": effect_ci_low,
                "effect_ci_high": effect_ci_high,
                "median_ratio": np.nan,
                "median_ratio_ci_low": np.nan,
                "median_ratio_ci_high": np.nan,
                "p_value": fisher_freeman_halton_2xk(exact_test_table.to_numpy()),
            }
        )

    sensitivity_categorical_composition = pd.DataFrame(sensitivity_categorical_rows)

    # Analysis block 29
    sensitivity_tests = pd.DataFrame(sensitivity_test_rows)
    if sensitivity_tests["feature"].tolist() != expected_features:
        raise AssertionError(
            "Sensitivity-test family changed: "
            f"expected {expected_features}, observed {sensitivity_tests['feature'].tolist()}"
        )

    sensitivity_tests["adj_p_value_bh"] = bh_adjust(sensitivity_tests["p_value"])
    sensitivity_tests["bh_significant_0_05"] = sensitivity_tests[
        "adj_p_value_bh"
    ].lt(0.05)

    sensitivity_continuous_path = OUT_DIR / "origin_exclusion_continuous_summary.csv"
    sensitivity_categorical_path = OUT_DIR / "origin_exclusion_categorical_composition.csv"
    sensitivity_continuous_summary.to_csv(sensitivity_continuous_path, index=False)
    sensitivity_categorical_composition.to_csv(sensitivity_categorical_path, index=False)
    print(f"Saved: {sensitivity_continuous_path}")
    print(f"Saved: {sensitivity_categorical_path}")

    # Analysis block 32
    import warnings

    import statsmodels.api as sm
    from scipy.stats import chi2

    scale_source = master.copy()
    scale_source["log1p_n_traits"] = np.log1p(scale_source["n_traits"])

    scale_specs = {
        "log10_sv_size_bp": "z_log10_sv_size_bp",
        "log10_maf": "z_log10_maf",
        "log1p_n_traits": "z_log1p_n_traits",
    }
    scaling_rows = []
    for source_variable, standardized_variable in scale_specs.items():
        mean = scale_source[source_variable].mean()
        standard_deviation = scale_source[source_variable].std(ddof=0)
        if not np.isfinite(standard_deviation) or standard_deviation <= 0:
            raise ValueError(f"Cannot standardize {source_variable}")
        scaling_rows.append(
            {
                "source_variable": source_variable,
                "standardized_variable": standardized_variable,
                "mean_full_set": mean,
                "sd_full_set": standard_deviation,
            }
        )

    predictor_scaling = pd.DataFrame(scaling_rows)

    # Analysis block 33
    def prepare_model_data(data):
        model_data = data.copy()
        model_data["log1p_n_traits"] = np.log1p(model_data["n_traits"])

        for row in predictor_scaling.itertuples(index=False):
            model_data[row.standardized_variable] = (
                model_data[row.source_variable] - row.mean_full_set
            ) / row.sd_full_set

        model_data["functional_overlap"] = model_data[
            "has_functional_annotation"
        ].astype(int)
        model_data["category_enriched_outcome"] = model_data[
            "is_category_enriched"
        ].astype(int)
        return model_data


    TERM_LABELS = {
        "z_log10_sv_size_bp": "SV size (per 1 SD log10)",
        "z_log10_maf": "MAF (per 1 SD log10)",
        "functional_overlap": "Functional overlap vs intronic only",
        "z_log1p_n_traits": "Trait count (per 1 SD log1p)",
    }


    def run_logistic_model(data, model_name, predictors):
        required_columns = ["category_enriched_outcome"] + list(predictors)
        if data[required_columns].isna().any().any():
            missing = data[required_columns].isna().sum()
            raise ValueError(f"Missing model values in {model_name}:\n{missing}")

        outcome = data["category_enriched_outcome"].astype(int)
        design = sm.add_constant(data[list(predictors)].astype(float), has_constant="add")

        with warnings.catch_warnings(record=True) as caught_warnings:
            warnings.simplefilter("always")
            result = sm.GLM(
                outcome, design, family=sm.families.Binomial()
            ).fit(maxiter=200)

        warning_text = " | ".join(str(item.message) for item in caught_warnings)
        confidence_interval = result.conf_int(alpha=0.05)
        coefficient_rows = []
        for term in predictors:
            coefficient_rows.append(
                {
                    "model": model_name,
                    "term": term,
                    "term_label": TERM_LABELS[term],
                    "beta": result.params[term],
                    "standard_error": result.bse[term],
                    "odds_ratio": np.exp(result.params[term]),
                    "or_ci_low": np.exp(confidence_interval.loc[term, 0]),
                    "or_ci_high": np.exp(confidence_interval.loc[term, 1]),
                    "p_value": result.pvalues[term],
                }
            )

        finite_estimates = bool(
            np.isfinite(result.params).all()
            and np.isfinite(result.bse).all()
            and np.isfinite(confidence_interval.to_numpy()).all()
        )
        separation_warning = "separation" in warning_text.lower()
        events = int(outcome.sum())
        n_predictors = len(predictors)
        likelihood_ratio = result.null_deviance - result.deviance
        likelihood_ratio_p = chi2.sf(likelihood_ratio, df=n_predictors)

        diagnostic_row = {
            "model": model_name,
            "n": len(data),
            "events": events,
            "non_events": len(data) - events,
            "n_predictors": n_predictors,
            "events_per_predictor": events / n_predictors,
            "converged": bool(result.converged),
            "finite_estimates": finite_estimates,
            "separation_warning": separation_warning,
            "warning_text": warning_text,
            "condition_number": np.linalg.cond(design.to_numpy()),
            "deviance": result.deviance,
            "null_deviance": result.null_deviance,
            "likelihood_ratio_statistic": likelihood_ratio,
            "likelihood_ratio_df": n_predictors,
            "likelihood_ratio_p_value": likelihood_ratio_p,
            "aic": result.aic,
            "model_usable": bool(
                result.converged and finite_estimates and not separation_warning
            ),
        }

        return (
            pd.DataFrame(coefficient_rows),
            pd.DataFrame([diagnostic_row]),
        )

    # Analysis block 35
    full_model_data = prepare_model_data(master)
    exclusion_model_data = prepare_model_data(sensitivity_master)

    primary_predictors = [
        "z_log10_sv_size_bp",
        "z_log10_maf",
        "functional_overlap",
    ]
    trait_count_sensitivity_predictors = primary_predictors + ["z_log1p_n_traits"]

    model_specs = [
        ("Full primary", full_model_data, primary_predictors),
        ("Exclusion primary", exclusion_model_data, primary_predictors),
        (
            "Full + trait-count sensitivity",
            full_model_data,
            trait_count_sensitivity_predictors,
        ),
        (
            "Exclusion + trait-count sensitivity",
            exclusion_model_data,
            trait_count_sensitivity_predictors,
        ),
    ]

    coefficient_tables = []
    diagnostic_tables = []

    for model_name, model_data, predictors in model_specs:
        coefficients, diagnostics = run_logistic_model(
            model_data, model_name, predictors
        )
        coefficient_tables.append(coefficients)
        diagnostic_tables.append(diagnostics)

    multivariable_coefficients = pd.concat(coefficient_tables, ignore_index=True)
    multivariable_diagnostics = pd.concat(diagnostic_tables, ignore_index=True)

    coefficients_path = OUT_DIR / "multivariable_logistic_coefficients.csv"
    diagnostics_path = OUT_DIR / "multivariable_model_diagnostics.csv"
    multivariable_coefficients.to_csv(coefficients_path, index=False)
    multivariable_diagnostics.to_csv(diagnostics_path, index=False)
    print(f"Saved: {coefficients_path}")
    print(f"Saved: {diagnostics_path}")

    feature_tests = pd.concat([
        primary_tests.assign(analysis_set="All recurrent lead SVs"),
        sensitivity_tests.assign(analysis_set="Exclude flagged lead SVs"),
    ], ignore_index=True)
    feature_tests.to_csv(RESULT_DIR / "sv_feature_tests.csv", index=False)


if __name__ == "__main__":
    main()
