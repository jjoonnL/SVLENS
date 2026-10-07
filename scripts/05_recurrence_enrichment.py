"""Test phenotype-module and category enrichment of recurrent exact lead SVs."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


N_PERMUTATIONS = 20_000
RANDOM_SEED = 20260626


def join_unique(values):
    return "; ".join(sorted(pd.Series(values).dropna().astype(str).unique()))


def bh_adjust(values):
    p = np.asarray(values, dtype=float)
    result = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    order = np.argsort(p[valid])
    ranked = p[valid][order]
    if not len(ranked):
        return result
    q = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    result[np.flatnonzero(valid)[order]] = np.clip(q, 0, 1)
    return result


def concentration(labels):
    counts = pd.Series(labels).dropna().astype(str).value_counts()
    n = int(counts.sum())
    proportions = counts / n
    return {
        "best": counts.index[0],
        "purity": float(counts.iloc[0] / n),
        "simpson": float((proportions ** 2).sum()),
        "counts": "; ".join(f"{key}:{value}" for key, value in counts.items()),
    }


def module_scenario(name, lead_traits, background, unit_col, rng, n_permutations):
    fine = background["trait_category"].to_numpy()
    broad = background["broad_module"].to_numpy()
    null_cache = {}

    def null_for(n):
        if n not in null_cache:
            draws = {
                "fine_purity": np.empty(n_permutations),
                "fine_simpson": np.empty(n_permutations),
                "broad_purity": np.empty(n_permutations),
                "broad_simpson": np.empty(n_permutations),
            }
            for index in range(n_permutations):
                selected = rng.choice(len(background), size=n, replace=False)
                fine_draw = concentration(fine[selected])
                broad_draw = concentration(broad[selected])
                for metric in ("purity", "simpson"):
                    draws[f"fine_{metric}"][index] = fine_draw[metric]
                    draws[f"broad_{metric}"][index] = broad_draw[metric]
            null_cache[n] = draws
        return null_cache[n]

    rows = []
    for lead_sv_id, group in lead_traits.groupby("lead_sv_id"):
        units = group[[unit_col, "trait_category", "broad_module"]].drop_duplicates()
        n = units[unit_col].nunique()
        if n < 2:
            continue
        fine_obs = concentration(units["trait_category"])
        broad_obs = concentration(units["broad_module"])
        null = null_for(n)
        rows.append({
            "scenario": name, "lead_sv_id": lead_sv_id,
            "genes": join_unique(group["gene_name"]), "n_units": n,
            "unit_col": unit_col, "units": join_unique(units[unit_col]),
            "traits": join_unique(group["trait"]),
            "fine_category_counts": fine_obs["counts"],
            "broad_module_counts": broad_obs["counts"],
            "best_fine_category": fine_obs["best"],
            "fine_purity": fine_obs["purity"],
            "fine_simpson": fine_obs["simpson"],
            "best_broad_module": broad_obs["best"],
            "broad_purity": broad_obs["purity"],
            "broad_simpson": broad_obs["simpson"],
            **{
                f"{layer}_{metric}_emp_p": (
                    1 + np.sum(null[f"{layer}_{metric}"] >= observed[metric])
                ) / (n_permutations + 1)
                for layer, observed in (("fine", fine_obs), ("broad", broad_obs))
                for metric in ("purity", "simpson")
            },
        })
    result = pd.DataFrame(rows)
    for layer in ("fine", "broad"):
        for metric in ("purity", "simpson"):
            result[f"{layer}_{metric}_fdr"] = bh_adjust(result[f"{layer}_{metric}_emp_p"])
        result[f"{layer}_fdr_lt_0_05"] = result[f"{layer}_purity_fdr"].lt(0.05)
    result["any_module_fdr_lt_0_05"] = (
        result["fine_fdr_lt_0_05"] | result["broad_fdr_lt_0_05"]
    )
    return result


def category_pairs(lead, background, unit_col, scenario):
    counts = background["trait_category"].value_counts()
    total = int(background[unit_col].nunique())
    rows = []
    for record in lead.itertuples(index=False):
        observed = {}
        for item in record.fine_category_counts.split("; "):
            name, count = item.rsplit(":", 1)
            observed[name] = int(count)
        for category, available in counts.items():
            count = observed.get(category, 0)
            expected = record.n_units * available / total
            rows.append({
                "scenario": scenario, "lead_sv_id": record.lead_sv_id,
                "genes": record.genes, "n_units": record.n_units,
                "category": category, "observed_count": count,
                "background_category_count": int(available),
                "background_total_units": total, "expected_count": expected,
                "fold_enrichment": count / expected if expected else np.nan,
                "hypergeom_p": float(hypergeom.sf(count - 1, total, available, record.n_units))
                if count else 1.0,
            })
    pairs = pd.DataFrame(rows)
    pairs["category_fdr"] = bh_adjust(pairs["hypergeom_p"])
    pairs["category_fdr_lt_0_05"] = pairs["category_fdr"].lt(0.05)
    return pairs


def category_assignment(pairs, recurrence):
    lead = recurrence.set_index("lead_sv_id")
    rows = []
    for lead_sv_id, group in pairs.groupby("lead_sv_id", sort=False):
        significant = group.loc[group["category_fdr_lt_0_05"]].sort_values(
            ["category_fdr", "hypergeom_p", "category"]
        )
        observed = group.loc[group["observed_count"].gt(0)].sort_values(
            ["observed_count", "category_fdr"], ascending=[False, True]
        )
        first = lead.loc[lead_sv_id]
        rows.append({
            "lead_sv_id": lead_sv_id,
            "genes": first["genes"],
            "n_traits_non_ratio": int(first["n_traits_non_ratio"]),
            "n_significant_categories": len(significant),
            "significant_categories": join_unique(significant["category"]),
            "significant_category_counts": "; ".join(
                f"{r.category}:{r.observed_count}" for r in significant.itertuples()
            ),
            "plot_group": (
                "Not enriched" if significant.empty else
                significant.iloc[0]["category"] if len(significant) == 1 else "Mixed"
            ),
            "min_category_fdr": float(group["category_fdr"].min()),
            "best_category_by_fdr": group.sort_values(
                ["category_fdr", "hypergeom_p"]
            ).iloc[0]["category"],
            "observed_categories": "; ".join(
                f"{r.category}:{r.observed_count}" for r in observed.itertuples()
            ),
            "lead_sv_type": first["lead_sv_type"],
            "lead_sv_layer": first["lead_sv_layer"],
            "lead_sv_size": first["lead_sv_size"],
            "lead_sv_maf": first["lead_sv_maf"],
        })
    return pd.DataFrame(rows).sort_values(
        ["n_significant_categories", "n_traits_non_ratio", "min_category_fdr"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def run(master, metadata, n_permutations=N_PERMUTATIONS):
    traits = master[["trait", "trait_category", "broad_module"]].drop_duplicates().sort_values("trait")
    traits = traits.merge(metadata[["trait", "family", "is_ratio"]], on="trait", validate="one_to_one")
    traits["is_ratio"] = traits["is_ratio"].astype(bool)
    non_ratio_traits = traits.loc[~traits["is_ratio"]].copy()
    family_background = (
        non_ratio_traits[["family", "trait_category", "broad_module"]]
        .drop_duplicates()
        .sort_values("family")
    )
    lead_counts = master.groupby("lead_sv_id")["trait"].nunique()
    recurrent_ids = set(lead_counts.loc[lead_counts.ge(2)].index)
    lead_traits = master.loc[master["lead_sv_id"].isin(recurrent_ids)].copy()
    lead_traits = lead_traits.merge(
        traits[["trait", "family", "is_ratio"]], on="trait", validate="many_to_one"
    )
    lead_traits = lead_traits.loc[~lead_traits["is_ratio"]]

    rng = np.random.default_rng(RANDOM_SEED)
    individual = module_scenario(
        "no_ratio", lead_traits,
        non_ratio_traits[["trait", "trait_category", "broad_module"]],
        "trait", rng, n_permutations,
    )
    family_units = (
        lead_traits.sort_values(["lead_sv_id", "family", "trait"])
        .drop_duplicates(["lead_sv_id", "family"])
    )
    collapsed = module_scenario(
        "family_no_ratio", family_units,
        family_background,
        "family", rng, n_permutations,
    )
    module = pd.concat([individual, collapsed], ignore_index=True)
    individual_pairs = category_pairs(individual, non_ratio_traits, "trait", "no_ratio")
    family_pairs = category_pairs(collapsed, family_background, "family", "family_no_ratio")
    pairs = pd.concat([individual_pairs, family_pairs], ignore_index=True)
    return module, pairs, traits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--traits", type=Path, default=Path(__file__).resolve().parents[1] / "metadata/traits.tsv")
    args = parser.parse_args()
    base = args.project_root / "results/sv_pleiotropy"
    master = pd.read_csv(base / "sv_pleiotropy_master.csv")
    metadata = pd.read_csv(args.traits, sep="\t")
    module, pairs, traits = run(master, metadata)
    recurrence_path = base / "lead_sv_level_pleiotropy_summary.non_ratio_primary.csv"
    recurrence = pd.read_csv(recurrence_path)
    individual_module = module.loc[module["scenario"].eq("no_ratio")].copy()
    module_columns = [
        "lead_sv_id", "n_units", "fine_category_counts", "broad_module_counts",
        "best_broad_module", "broad_purity", "broad_purity_emp_p", "broad_purity_fdr",
        "best_fine_category", "fine_purity", "fine_purity_emp_p", "fine_purity_fdr",
        "broad_fdr_lt_0_05", "fine_fdr_lt_0_05", "any_module_fdr_lt_0_05",
    ]
    recurrence = recurrence.merge(
        individual_module[module_columns].rename(columns={"n_units": "n_units_no_ratio_permutation"}),
        on="lead_sv_id", how="left", validate="one_to_one",
    )
    recurrence.to_csv(recurrence_path, index=False)
    recurrence.loc[recurrence["is_recurrent_non_ratio"]].to_csv(
        base / "candidate_recurrent_lead_svs.non_ratio_primary.csv", index=False
    )
    output = base / "category_enrichment"
    output.mkdir(parents=True, exist_ok=True)
    module.to_csv(base / "lead_sv_redundancy_sensitivity.csv", index=False)
    pairs.to_csv(output / "recurrent_lead_sv_category_enrichment_sensitivity.csv", index=False)
    individual_pairs = pairs.loc[pairs["scenario"].eq("no_ratio")].copy()
    assignment = category_assignment(individual_pairs, recurrence)
    individual_pairs.rename(columns={
        "n_units": "n_traits_non_ratio", "background_category_count": "background_trait_count",
        "background_total_units": "background_total_traits",
    }).to_csv(output / "recurrent_lead_sv_category_enrichment.csv", index=False)
    assignment.to_csv(output / "recurrent_lead_sv_category_assignment.csv", index=False)
    phenotype = traits.rename(columns={
        "trait_category": "trait_category", "family": "collapsed_family",
        "is_ratio": "is_ratio_trait",
    })
    phenotype["used_in_non_ratio_primary"] = ~phenotype["is_ratio_trait"]
    phenotype_dir = base / "phenotype_grouping"
    phenotype_dir.mkdir(parents=True, exist_ok=True)
    phenotype.to_csv(phenotype_dir / "phenotype_family_mapping.csv", index=False)
    print(f"{len(module)} module rows and {len(pairs)} category-pair rows")


if __name__ == "__main__":
    main()
