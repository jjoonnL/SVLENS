"""Classify SV gene–trait associations by lead-SV removal."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from svarch.acat import run_acat_gene


P_FLOOR = 1e-300


def safe_logp(value):
    return np.nan if pd.isna(value) else -np.log10(max(float(value), P_FLOOR))


def summarize_gene_svs(gene_sv):
    valid = gene_sv.loc[gene_sv["w_final"].notna()]
    if valid.empty:
        return pd.DataFrame()
    per_sv = (
        valid.groupby("sv_id", dropna=False)
        .agg(
            sv_type=("sv_type", "first"), svsize=("svsize", "first"),
            maf=("maf", "first"), min_pval=("pval", "min"),
            n_rows=("sv_id", "size"),
        )
        .reset_index()
        .sort_values(["min_pval", "sv_id"])
        .reset_index(drop=True)
    )
    per_sv["rank_by_pval"] = np.arange(1, len(per_sv) + 1)
    return per_sv


def without_lead_acat(gene_sv, lead_sv_id, gene_id, gene_name):
    remaining = gene_sv.loc[gene_sv["sv_id"].ne(lead_sv_id)]
    if remaining.empty or remaining["w_final"].notna().sum() == 0:
        return np.nan, 0, 0
    _, genes = run_acat_gene(remaining)
    if genes.empty:
        return np.nan, 0, 0
    match = genes.loc[genes["gene_id"].eq(gene_id)]
    if match.empty:
        match = genes.loc[genes["gene_name"].eq(gene_name)]
    if match.empty:
        return np.nan, 0, int(remaining["w_final"].notna().sum())
    row = match.iloc[0]
    return float(row["p_acat_o"]), int(row["n_strata"]), int(row["n_sv_total"])


def decompose_hit(hit, weighted):
    gene_sv = weighted.loc[weighted["gene_id"].eq(hit["gene_id"])].copy()
    if gene_sv.empty:
        gene_sv = weighted.loc[weighted["gene_name"].eq(hit["gene_name"])].copy()
    row = hit.to_dict()
    row["status"] = "ok"
    row["n_sv_rows_original"] = len(gene_sv)
    row["n_sv_valid_original"] = int(gene_sv["w_final"].notna().sum()) if not gene_sv.empty else 0
    if gene_sv.empty:
        row["status"] = "gene_missing_in_weighted"
        row["driver_class"] = "lead_missing_or_error"
        return row

    per_sv = summarize_gene_svs(gene_sv)
    if per_sv.empty:
        row["status"] = "no_valid_sv_for_gene"
        row["driver_class"] = "lead_missing_or_error"
        return row
    row["n_unique_sv_original"] = len(per_sv)
    best = per_sv.iloc[0]
    second = per_sv.iloc[1] if len(per_sv) > 1 else None
    row["best_sv_id"] = best["sv_id"]
    row["best_sv_pval"] = float(best["min_pval"])
    row["best_sv_type"] = best["sv_type"]
    row["best_sv_maf"] = float(best["maf"])
    row["second_sv_id"] = second["sv_id"] if second is not None else np.nan
    row["second_sv_pval"] = float(second["min_pval"]) if second is not None else np.nan
    row["second_sv_type"] = second["sv_type"] if second is not None else np.nan
    row["second_sv_maf"] = float(second["maf"]) if second is not None else np.nan
    lead_rows = per_sv.loc[per_sv["sv_id"].eq(hit["lead_sv_id"])]
    if lead_rows.empty:
        row["status"] = "lead_missing_in_weighted"
        row["driver_class"] = "lead_missing_or_error"
        return row

    lead = lead_rows.iloc[0]
    row["lead_rank_by_pval"] = int(lead["rank_by_pval"])
    row["lead_pval_weighted"] = float(lead["min_pval"])
    row["lead_maf_weighted"] = float(lead["maf"])
    row["lead_sv_type_weighted"] = lead["sv_type"]
    row["lead_sv_size_weighted"] = int(lead["svsize"])
    row["lead_neglog10_pval"] = safe_logp(lead["min_pval"])
    row["second_neglog10_pval"] = safe_logp(row["second_sv_pval"])
    row["lead_second_log10_margin"] = (
        row["lead_neglog10_pval"] - row["second_neglog10_pval"]
        if pd.notna(row["second_neglog10_pval"]) else np.nan
    )
    p_without, n_strata, n_sv = without_lead_acat(
        gene_sv, hit["lead_sv_id"], hit["gene_id"], hit["gene_name"]
    )
    row["p_acat_o_without_lead"] = p_without
    row["n_strata_without_lead"] = n_strata
    row["n_sv_without_lead"] = n_sv
    threshold = float(hit["sv_bonf_threshold"])
    row["original_bonf_significant"] = float(hit["p_acat_o"]) < threshold
    row["without_lead_bonf_significant"] = bool(pd.notna(p_without) and p_without < threshold)
    row["original_neglog10_p"] = safe_logp(hit["p_acat_o"])
    row["without_lead_neglog10_p"] = safe_logp(p_without)
    row["delta_neglog10_original_minus_without"] = (
        row["original_neglog10_p"] - row["without_lead_neglog10_p"]
        if pd.notna(row["without_lead_neglog10_p"]) else np.nan
    )
    row["lead_not_best_by_pval"] = row["lead_rank_by_pval"] != 1
    row["low_margin_to_second_sv"] = bool(
        pd.notna(row["lead_second_log10_margin"])
        and row["lead_second_log10_margin"] < 1.0
    )
    row["driver_class"] = (
        "single_sv_only" if row["n_sv_valid_original"] <= 1 or pd.isna(p_without)
        else "multi_sv_supported" if row["without_lead_bonf_significant"]
        else "lead_anchored"
    )
    return row


def run(master, results_dir):
    required = {"trait", "gene_id", "gene_name", "p_acat_o", "sv_bonf_threshold", "lead_sv_id"}
    missing = required - set(master.columns)
    if missing:
        raise ValueError(f"Association master is missing columns: {sorted(missing)}")
    rows = []
    for trait, hits in master.groupby("trait", sort=True):
        path = results_dir / trait / "sv_weighted.parquet"
        weighted = pd.read_parquet(path)
        rows.extend(decompose_hit(hit, weighted) for _, hit in hits.iterrows())
    return pd.DataFrame(rows)


def semicolon_join(values):
    return "; ".join(sorted(pd.Series(values).dropna().astype(str).unique()))


def summarize_by_lead(detail):
    summary = (
        detail.groupby(["lead_sv_id", "gene_name"], dropna=False)
        .agg(
            n_gene_trait_hits=("trait", "nunique"),
            traits=("trait", semicolon_join),
            driver_classes=("driver_class", semicolon_join),
            n_single_sv_only=("driver_class", lambda x: int(x.eq("single_sv_only").sum())),
            n_lead_anchored=("driver_class", lambda x: int(x.eq("lead_anchored").sum())),
            n_multi_sv_supported=("driver_class", lambda x: int(x.eq("multi_sv_supported").sum())),
            n_lead_missing_or_error=("driver_class", lambda x: int(x.eq("lead_missing_or_error").sum())),
            median_n_unique_sv_original=("n_unique_sv_original", "median"),
            median_lead_rank_by_pval=("lead_rank_by_pval", "median"),
            median_lead_second_log10_margin=("lead_second_log10_margin", "median"),
            median_delta_neglog10_original_minus_without=("delta_neglog10_original_minus_without", "median"),
            min_without_lead_p=("p_acat_o_without_lead", "min"),
            any_without_lead_bonf=("without_lead_bonf_significant", "max"),
            any_lead_not_best=("lead_not_best_by_pval", "max"),
            any_low_margin_to_second=("low_margin_to_second_sv", "max"),
            lead_sv_maf=("lead_sv_maf", "first"),
            lead_sv_type=("lead_sv_type", "first"),
            lead_sv_layer=("lead_sv_layer", semicolon_join),
            support_classes=("support_class_final", semicolon_join),
        )
        .reset_index()
    )
    summary["frac_lead_anchored_or_single"] = (
        (summary["n_single_sv_only"] + summary["n_lead_anchored"])
        / summary["n_gene_trait_hits"]
    )
    summary["frac_multi_sv_supported"] = summary["n_multi_sv_supported"] / summary["n_gene_trait_hits"]
    return summary.sort_values(
        ["n_gene_trait_hits", "frac_lead_anchored_or_single"], ascending=[False, False]
    ).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()
    results_dir = args.project_root / "results"
    output = results_dir / "sv_pleiotropy" / "driver_decomposition"
    master = pd.read_csv(results_dir / "sv_pleiotropy" / "sv_pleiotropy_master.non_ratio_primary.csv")
    detail = run(master, results_dir)
    if len(detail) != len(master) or detail["driver_class"].eq("lead_missing_or_error").any():
        raise AssertionError("Lead-SV decomposition is incomplete")
    output.mkdir(parents=True, exist_ok=True)
    detail.to_csv(output / "lead_sv_driver_decomposition.non_ratio_primary.csv", index=False)
    summary = detail.groupby("driver_class").size().rename("n").reset_index()
    summary["pct"] = 100 * summary["n"] / len(detail)
    summary.to_csv(output / "driver_class_summary.csv", index=False)
    summarize_by_lead(detail).to_csv(
        output / "lead_sv_driver_summary.non_ratio_primary.csv", index=False
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
