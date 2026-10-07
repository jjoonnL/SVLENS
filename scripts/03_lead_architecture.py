"""Build significant SV associations, classify lead removal, and audit residual SVs."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from svarch.acat import SV_TYPES, acat_pvalue, run_acat_gene


# Association master and exact lead-SV recurrence
def annotation_layer(row):
    if bool(row["is_cds"]) or bool(row["is_utr"]) or bool(row["is_promoter"]):
        return "Functional"
    return "Intronic" if bool(row["is_intronic"]) else "Other"


def sv_representation_tag(sv_id):
    match = re.search(r"<([^:>]+):SVSIZE=([^:>]+):([^>]+)>", str(sv_id))
    return match.group(3) if match else None


def join_unique(values):
    return "; ".join(sorted(pd.Series(values).dropna().astype(str).unique()))


def read_significant_trait(trait_dir, category):
    acat = pd.read_parquet(trait_dir / "acat_gene.parquet")
    threshold = 0.05 / len(acat)
    significant = acat.loc[acat["p_acat_o"].lt(threshold)].copy()
    if significant.empty:
        return significant

    weighted = pd.read_parquet(trait_dir / "sv_weighted.parquet").copy()
    weighted["lead_sv_layer"] = weighted.apply(annotation_layer, axis=1)
    lead = (
        weighted.sort_values(["gene_id", "pval", "maf"])
        .drop_duplicates("gene_id", keep="first")
        .rename(columns={
            "sv_id": "lead_sv_id", "chrom": "lead_sv_chrom",
            "sv_start": "lead_sv_start", "sv_end": "lead_sv_end",
            "sv_type": "lead_sv_type", "svsize": "lead_sv_size",
            "Beta": "lead_sv_beta", "SE": "lead_sv_se",
            "pval": "lead_sv_p", "maf": "lead_sv_maf",
        })
    )
    lead_columns = [
        "gene_id", "lead_sv_id", "lead_sv_chrom", "lead_sv_start", "lead_sv_end",
        "lead_sv_type", "lead_sv_size", "lead_sv_beta", "lead_sv_se",
        "lead_sv_p", "lead_sv_maf", "is_cds", "is_utr", "is_promoter",
        "is_intronic", "lead_sv_layer",
    ]
    out = significant.merge(lead[lead_columns], on="gene_id", validate="one_to_one")
    out["trait"] = trait_dir.name
    out["trait_category"] = category
    out["sv_bonf_threshold"] = threshold
    out["sv_bonf"] = True
    out["lead_sv_pos"] = (
        (out["lead_sv_start"] + out["lead_sv_end"]) / 2
    ).round().astype("Int64")
    out["sv_repr_tag"] = out["lead_sv_id"].map(sv_representation_tag)
    return out


def build_master(results_dir, metadata):
    category = metadata.set_index("trait")["category"].to_dict()
    trait_dirs = sorted(
        path for path in results_dir.iterdir()
        if path.is_dir() and (path / "acat_gene.parquet").exists()
        and (path / "sv_weighted.parquet").exists()
    )
    unknown = {path.name for path in trait_dirs} - set(category)
    if unknown:
        raise ValueError(f"Traits missing from metadata: {sorted(unknown)}")
    missing = set(category) - {path.name for path in trait_dirs}
    if missing:
        raise ValueError(f"Missing per-trait ACAT results: {sorted(missing)}")
    rows = [read_significant_trait(path, category[path.name]) for path in trait_dirs]
    nonempty = [row for row in rows if not row.empty]
    if not nonempty:
        raise ValueError("No significant SV gene–trait associations")
    master = pd.concat(nonempty, ignore_index=True)
    genes = pd.read_parquet(results_dir / "gene_table.parquet")
    gene_columns = [
        "gene_id", "chrom", "gene_start", "gene_end",
        "promoter_start", "promoter_end",
    ]
    master = master.merge(genes[gene_columns], on="gene_id", validate="many_to_one")
    master["gene_length"] = master["gene_end"] - master["gene_start"]
    master["support_class_final"] = "not_assessed"
    return master


def exact_lead_recurrence(master):
    lead = (
        master.groupby("lead_sv_id", dropna=False)
        .agg(
            genes=("gene_name", join_unique), n_genes=("gene_name", "nunique"),
            n_traits_non_ratio=("trait", "nunique"),
            traits_non_ratio=("trait", join_unique),
            fine_categories_non_ratio=("trait_category", join_unique),
            broad_modules_non_ratio=("broad_module", join_unique),
            lead_sv_chrom=("lead_sv_chrom", "first"),
            lead_sv_pos=("lead_sv_pos", "median"),
            lead_sv_type=("lead_sv_type", "first"),
            lead_sv_layer=("lead_sv_layer", join_unique),
            lead_sv_size=("lead_sv_size", "first"),
            lead_sv_maf=("lead_sv_maf", "median"),
            min_lead_sv_p_non_ratio=("lead_sv_p", "min"),
            min_p_acat_o_non_ratio=("p_acat_o", "min"),
            support_classes_non_ratio=("support_class_final", join_unique),
        )
        .reset_index()
    )
    lead["is_recurrent_non_ratio"] = lead["n_traits_non_ratio"].ge(2)
    return lead.sort_values(
        ["is_recurrent_non_ratio", "n_traits_non_ratio", "min_p_acat_o_non_ratio"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def _stage_primary(args):
    metadata = pd.read_csv(args.traits, sep="\t")
    results_dir = args.project_root / "results"
    output = results_dir / "sv_pleiotropy"
    master = build_master(results_dir, metadata)
    metadata_indexed = metadata.set_index("trait")
    master["broad_module"] = master["trait"].map(metadata_indexed["broad_module"])
    master["is_ratio_trait"] = master["trait"].map(metadata_indexed["is_ratio"]).astype(bool)
    non_ratio = master.loc[~master["is_ratio_trait"]].copy()
    recurrence = exact_lead_recurrence(non_ratio)
    all_counts = master.groupby("lead_sv_id")["trait"].nunique()
    recurrence["n_traits_all"] = recurrence["lead_sv_id"].map(all_counts)
    recurrence["n_ratio_trait_hits_removed"] = (
        recurrence["n_traits_all"] - recurrence["n_traits_non_ratio"]
    )
    output.mkdir(parents=True, exist_ok=True)
    master.to_csv(output / "sv_pleiotropy_master.csv", index=False)
    non_ratio.to_csv(output / "sv_pleiotropy_master.non_ratio_primary.csv", index=False)
    recurrence.to_csv(output / "lead_sv_level_pleiotropy_summary.non_ratio_primary.csv", index=False)
    recurrence.loc[recurrence["is_recurrent_non_ratio"]].to_csv(
        output / "candidate_recurrent_lead_svs.non_ratio_primary.csv", index=False
    )
    print(
        f"{len(non_ratio)} associations, {non_ratio['gene_id'].nunique()} genes, "
        f"{non_ratio['lead_sv_id'].nunique()} exact lead SVs"
    )


# Lead-SV removal classification
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


def _stage_decomposition(args):
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


# Random-removal and conditional residual-SV audit
N_PERMUTATIONS = 100_000


RANDOM_SEED = 20260914


def gene_acat(table: pd.DataFrame) -> float:
    """Reproduce the two-level SV-type-stratified gene-level ACAT."""
    if table.empty:
        return np.nan
    stratum_p = []
    for sv_type in SV_TYPES:
        group = table.loc[table["sv_type"].eq(sv_type)]
        if group.empty:
            continue
        stratum_p.append(
            acat_pvalue(
                group["pval"].to_numpy(dtype=float),
                weights=group["w_final"].to_numpy(dtype=float),
            )
        )
    if not stratum_p:
        return np.nan
    return float(acat_pvalue(np.asarray(stratum_p), weights=None))


def interval_distance(
    start_a: int, end_a: int, start_b: int, end_b: int
) -> int:
    return max(0, start_a - end_b, start_b - end_a)


def load_trait_sv(trait: str) -> pd.DataFrame:
    path = RESULTS_DIR / trait / "sv_weighted.parquet"
    if not path.exists():
        raise FileNotFoundError(path)
    columns = [
        "sv_id", "sv_type", "sv_start", "sv_end", "svsize", "maf", "pval",
        "w_final", "gene_id", "gene_name",
    ]
    return pd.read_parquet(path, columns=columns)


def audit_associations(
    driver: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    removal_rows: list[dict[str, object]] = []
    association_rows: list[dict[str, object]] = []
    locality_rows: list[dict[str, object]] = []

    for trait, trait_hits in driver.groupby("trait", sort=True):
        trait_sv = load_trait_sv(trait)
        for hit in trait_hits.itertuples(index=False):
            gene_sv = trait_sv.loc[trait_sv["gene_id"].eq(hit.gene_id)].copy()
            gene_sv = (
                gene_sv.sort_values(["pval", "maf", "sv_id"])
                .drop_duplicates("sv_id", keep="first")
                .reset_index(drop=True)
            )
            if gene_sv.empty:
                raise ValueError(f"No SV rows for {trait}, {hit.gene_id}")
            if hit.lead_sv_id not in set(gene_sv["sv_id"]):
                raise ValueError(f"Lead SV absent for {trait}, {hit.gene_id}")

            full_p = gene_acat(gene_sv)
            threshold = float(hit.sv_bonf_threshold)
            lead = gene_sv.loc[gene_sv["sv_id"].eq(hit.lead_sv_id)].iloc[0]
            association_removals: list[dict[str, object]] = []

            for removed in gene_sv.itertuples(index=False):
                remaining = gene_sv.loc[gene_sv["sv_id"].ne(removed.sv_id)]
                p_without = gene_acat(remaining)
                loses_significance = bool(
                    np.isnan(p_without) or p_without >= threshold
                )
                is_lead = removed.sv_id == hit.lead_sv_id
                distance_to_lead = (
                    0 if is_lead else interval_distance(
                        int(lead["sv_start"]), int(lead["sv_end"]),
                        int(removed.sv_start), int(removed.sv_end),
                    )
                )
                row = {
                    "trait": trait,
                    "gene_id": hit.gene_id,
                    "gene_name": hit.gene_name,
                    "lead_sv_id": hit.lead_sv_id,
                    "removed_sv_id": removed.sv_id,
                    "removed_sv_type": removed.sv_type,
                    "removed_sv_p": float(removed.pval),
                    "is_lead": is_lead,
                    "p_without_removed_sv": p_without,
                    "loses_significance": loses_significance,
                    "distance_to_lead_bp": distance_to_lead,
                    "same_type_as_lead": removed.sv_type == lead["sv_type"],
                    "n_svs": len(gene_sv),
                    "full_gene_p_recomputed": full_p,
                    "bonferroni_threshold": threshold,
                }
                removal_rows.append(row)
                association_removals.append(row)

            removal = pd.DataFrame(association_removals)
            lead_removal = removal.loc[removal["is_lead"]].iloc[0]
            nonlead = removal.loc[~removal["is_lead"]]
            n_critical = int(removal["loses_significance"].sum())
            n_nonlead_critical = int(nonlead["loses_significance"].sum())
            recalculated_class = (
                "single_sv_only" if len(gene_sv) == 1
                else "lead_anchored" if lead_removal["loses_significance"]
                else "multi_sv_supported"
            )

            association_rows.append({
                "trait": trait,
                "trait_category": hit.trait_category,
                "gene_id": hit.gene_id,
                "gene_name": hit.gene_name,
                "lead_sv_id": hit.lead_sv_id,
                "lead_sv_type": lead["sv_type"],
                "n_svs": len(gene_sv),
                "full_gene_p_reported": float(hit.p_acat_o),
                "full_gene_p_recomputed": full_p,
                "bonferroni_threshold": threshold,
                "reported_driver_class": hit.driver_class,
                "recalculated_driver_class": recalculated_class,
                "lead_is_critical": bool(lead_removal["loses_significance"]),
                "n_critical_removals": n_critical,
                "n_nonlead_critical_removals": n_nonlead_critical,
                "critical_removal_fraction": n_critical / len(gene_sv),
                "nonlead_critical_fraction": (
                    n_nonlead_critical / len(nonlead) if len(nonlead) else np.nan
                ),
                "lead_is_sole_critical_sv": bool(
                    lead_removal["loses_significance"] and n_nonlead_critical == 0
                ),
                "lead_second_log10_margin": hit.lead_second_log10_margin,
                "delta_neglog10_after_lead_removal": (
                    hit.delta_neglog10_original_minus_without
                ),
            })

            if recalculated_class == "multi_sv_supported":
                residual_candidates = gene_sv.loc[
                    gene_sv["sv_id"].ne(hit.lead_sv_id)
                ].copy()
                residual_candidates["distance_to_lead_bp"] = residual_candidates.apply(
                    lambda row: interval_distance(
                        int(lead["sv_start"]), int(lead["sv_end"]),
                        int(row["sv_start"]), int(row["sv_end"]),
                    ),
                    axis=1,
                )
                residual_candidates["same_type_as_lead"] = residual_candidates[
                    "sv_type"
                ].eq(lead["sv_type"])
                strongest = residual_candidates.sort_values(
                    ["pval", "maf", "sv_id"]
                ).iloc[0]
                observed_distance = int(strongest["distance_to_lead_bp"])
                distances = residual_candidates["distance_to_lead_bp"].to_numpy()
                locality_rows.append({
                    "trait": trait,
                    "trait_category": hit.trait_category,
                    "gene_id": hit.gene_id,
                    "gene_name": hit.gene_name,
                    "gene_length_bp": int(hit.gene_length),
                    "lead_sv_id": hit.lead_sv_id,
                    "lead_sv_type": lead["sv_type"],
                    "residual_sv_id": strongest["sv_id"],
                    "residual_sv_type": strongest["sv_type"],
                    "residual_sv_p": float(strongest["pval"]),
                    "observed_distance_bp": observed_distance,
                    "observed_distance_over_gene_length": (
                        observed_distance / float(hit.gene_length)
                        if hit.gene_length > 0 else np.nan
                    ),
                    "observed_overlap": observed_distance == 0,
                    "observed_within_50kb": observed_distance <= 50_000,
                    "observed_within_250kb": observed_distance <= 250_000,
                    "observed_same_type": bool(strongest["same_type_as_lead"]),
                    "observed_same_type_within_50kb": bool(
                        strongest["same_type_as_lead"] and observed_distance <= 50_000
                    ),
                    "n_nonlead_candidates": len(residual_candidates),
                    "candidate_median_distance_bp": float(np.median(distances)),
                    "candidate_mean_distance_bp": float(np.mean(distances)),
                    "candidate_probability_overlap": float(np.mean(distances == 0)),
                    "candidate_probability_within_50kb": float(np.mean(distances <= 50_000)),
                    "candidate_probability_within_250kb": float(np.mean(distances <= 250_000)),
                    "candidate_probability_same_type": float(np.mean(
                        residual_candidates["same_type_as_lead"].to_numpy()
                    )),
                    "candidate_probability_same_type_within_50kb": float(np.mean(
                        residual_candidates["same_type_as_lead"].to_numpy()
                        & (distances <= 50_000)
                    )),
                    "distance_lower_tail_fraction": float(np.mean(distances <= observed_distance)),
                    "candidate_ids": "; ".join(residual_candidates["sv_id"].astype(str)),
                    "candidate_distances_bp": "; ".join(distances.astype(str)),
                    "candidate_same_type": "; ".join(
                        residual_candidates["same_type_as_lead"].astype(int).astype(str)
                    ),
                })

    return (
        pd.DataFrame(removal_rows),
        pd.DataFrame(association_rows),
        pd.DataFrame(locality_rows),
    )


def criticality_permutation(
    association: pd.DataFrame, rng: np.random.Generator
) -> tuple[pd.DataFrame, pd.DataFrame]:
    eligible = association.loc[association["n_svs"].gt(1)].copy()
    probabilities = eligible["critical_removal_fraction"].to_numpy(dtype=float)
    null_counts = np.zeros(N_PERMUTATIONS, dtype=np.int32)
    for probability in probabilities:
        null_counts += rng.random(N_PERMUTATIONS) < probability

    observed = int(eligible["lead_is_critical"].sum())
    null_mean = float(null_counts.mean())
    summary = pd.DataFrame([{
        "metric": "associations_losing_significance_after_one_removal",
        "n_associations": len(eligible),
        "observed_lead_removal_count": observed,
        "random_removal_null_mean": null_mean,
        "observed_to_null_ratio": observed / null_mean if null_mean else np.inf,
        "empirical_upper_tail_p": (1 + int((null_counts >= observed).sum()))
        / (N_PERMUTATIONS + 1),
        "n_lead_as_sole_critical_sv": int(eligible["lead_is_sole_critical_sv"].sum()),
        "n_with_any_nonlead_critical_sv": int(
            eligible["n_nonlead_critical_removals"].gt(0).sum()
        ),
    }])
    distribution = pd.DataFrame({"random_removal_loss_count": null_counts})
    return summary, distribution


def parse_semicolon_ints(value: str) -> np.ndarray:
    return np.asarray([int(item) for item in str(value).split("; ")], dtype=int)


def parse_semicolon_bools(value: str) -> np.ndarray:
    return np.asarray([item == "1" for item in str(value).split("; ")], dtype=bool)


def locality_permutation(
    locality: pd.DataFrame, rng: np.random.Generator
) -> tuple[pd.DataFrame, pd.DataFrame]:
    null_overlap = np.zeros(N_PERMUTATIONS, dtype=np.int16)
    null_within_50 = np.zeros(N_PERMUTATIONS, dtype=np.int16)
    null_within_250 = np.zeros(N_PERMUTATIONS, dtype=np.int16)
    null_same_type = np.zeros(N_PERMUTATIONS, dtype=np.int16)
    null_same_type_50 = np.zeros(N_PERMUTATIONS, dtype=np.int16)
    null_total_distance = np.zeros(N_PERMUTATIONS, dtype=np.float64)

    for row in locality.itertuples(index=False):
        distances = parse_semicolon_ints(row.candidate_distances_bp)
        same_type = parse_semicolon_bools(row.candidate_same_type)
        selected = rng.integers(0, len(distances), size=N_PERMUTATIONS)
        sampled_distance = distances[selected]
        sampled_same_type = same_type[selected]
        null_overlap += sampled_distance == 0
        null_within_50 += sampled_distance <= 50_000
        null_within_250 += sampled_distance <= 250_000
        null_same_type += sampled_same_type
        null_same_type_50 += sampled_same_type & (sampled_distance <= 50_000)
        null_total_distance += sampled_distance

    null_mean_distance = null_total_distance / len(locality)
    distribution = pd.DataFrame({
        "overlap_count": null_overlap,
        "within_50kb_count": null_within_50,
        "within_250kb_count": null_within_250,
        "same_type_count": null_same_type,
        "same_type_within_50kb_count": null_same_type_50,
        "mean_distance_bp": null_mean_distance,
    })

    metric_specs = [
        ("overlap", "observed_overlap", "overlap_count", "upper"),
        ("within_50kb", "observed_within_50kb", "within_50kb_count", "upper"),
        ("within_250kb", "observed_within_250kb", "within_250kb_count", "upper"),
        ("same_type", "observed_same_type", "same_type_count", "upper"),
        (
            "same_type_within_50kb", "observed_same_type_within_50kb",
            "same_type_within_50kb_count", "upper",
        ),
    ]
    rows: list[dict[str, object]] = []
    for metric, observed_column, null_column, tail in metric_specs:
        observed = int(locality[observed_column].sum())
        null_values = distribution[null_column].to_numpy()
        rows.append({
            "metric": metric,
            "n_associations": len(locality),
            "observed": observed,
            "null_mean": float(null_values.mean()),
            "observed_to_null_ratio": observed / null_values.mean()
            if null_values.mean() else np.inf,
            "empirical_p": (1 + int((null_values >= observed).sum()))
            / (N_PERMUTATIONS + 1),
            "tail": tail,
        })

    observed_mean_distance = float(locality["observed_distance_bp"].mean())
    rows.append({
        "metric": "mean_distance_bp",
        "n_associations": len(locality),
        "observed": observed_mean_distance,
        "null_mean": float(null_mean_distance.mean()),
        "observed_to_null_ratio": observed_mean_distance / null_mean_distance.mean()
        if null_mean_distance.mean() else np.nan,
        "empirical_p": (1 + int((null_mean_distance <= observed_mean_distance).sum()))
        / (N_PERMUTATIONS + 1),
        "tail": "lower",
    })
    return pd.DataFrame(rows), distribution


def _stage_criticality(args):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    driver = pd.read_csv(DRIVER_PATH)
    if len(driver) != 404:
        raise ValueError(f"Expected 404 primary associations, found {len(driver)}")

    removal, association, locality = audit_associations(driver)
    rng = np.random.default_rng(RANDOM_SEED)
    criticality_summary, criticality_null = criticality_permutation(association, rng)
    locality_summary, locality_null = locality_permutation(locality, rng)

    driver_counts = association["recalculated_driver_class"].value_counts()
    decision = pd.DataFrame([
        {
            "question": "Is lead removal more disruptive than random member removal?",
            "result": (
                f"{int(criticality_summary.iloc[0]['observed_lead_removal_count'])} "
                f"observed vs {criticality_summary.iloc[0]['random_removal_null_mean']:.2f} "
                "expected"
            ),
            "interpretation_limit": (
                "Supports statistical concentration on the minimum-P SV; does not "
                "separate biology from the behavior of ACAT under sparse signals."
            ),
        },
        {
            "question": "Are strongest residual SVs unusually local within their gene sets?",
            "result": (
                f"{int(locality['observed_within_250kb'].sum())}/"
                f"{len(locality)} within 250 kb; see conditional geometry null"
            ),
            "interpretation_limit": (
                "Conditions on observed within-gene candidate SVs; does not establish "
                "independent alleles, LD, haplotypes, or causality."
            ),
        },
    ])

    qc = pd.DataFrame([
        {"check": "404 primary associations", "value": len(association),
         "expected": 404, "passed": len(association) == 404},
        {"check": "all reported associations remain significant on recomputation",
         "value": int((association["full_gene_p_recomputed"] < association["bonferroni_threshold"]).sum()),
         "expected": 404,
         "passed": bool((association["full_gene_p_recomputed"] < association["bonferroni_threshold"]).all())},
        {"check": "driver classes reproduced", "value": int(
            association["reported_driver_class"].eq(association["recalculated_driver_class"]).sum()
        ), "expected": 404, "passed": bool(
            association["reported_driver_class"].eq(association["recalculated_driver_class"]).all()
        )},
        {"check": "single-SV-only count", "value": int(driver_counts.get("single_sv_only", 0)),
         "expected": 8, "passed": int(driver_counts.get("single_sv_only", 0)) == 8},
        {"check": "lead-anchored count", "value": int(driver_counts.get("lead_anchored", 0)),
         "expected": 339, "passed": int(driver_counts.get("lead_anchored", 0)) == 339},
        {"check": "multi-SV-supported count", "value": len(locality),
         "expected": 57, "passed": len(locality) == 57},
        {"check": "one removal row per qualifying SV", "value": len(removal),
         "expected": int(association["n_svs"].sum()),
         "passed": len(removal) == int(association["n_svs"].sum())},
        {"check": "100000 criticality permutations", "value": len(criticality_null),
         "expected": N_PERMUTATIONS, "passed": len(criticality_null) == N_PERMUTATIONS},
        {"check": "100000 locality permutations", "value": len(locality_null),
         "expected": N_PERMUTATIONS, "passed": len(locality_null) == N_PERMUTATIONS},
    ])

    outputs = {
        "step33_leave_one_sv_out_detail.csv.gz": removal,
        "step33_association_criticality.csv": association,
        "step33_criticality_summary.csv": criticality_summary,
        "step33_criticality_null_distribution.csv.gz": criticality_null,
        "step33_locality_association_detail.csv": locality,
        "step33_locality_summary.csv": locality_summary,
        "step33_locality_null_distribution.csv.gz": locality_null,
        "step33_decision_summary.csv": decision,
        "step33_qc.csv": qc,
    }
    for filename, table in outputs.items():
        table.to_csv(OUT_DIR / filename, index=False)

    print("Driver classes")
    print(driver_counts.to_string())
    print("\nLead criticality")
    print(criticality_summary.to_string(index=False))
    print("\nConditional locality")
    print(locality_summary.to_string(index=False))
    print("\nQC")
    print(qc.to_string(index=False))
    if not qc["passed"].all():
        failed = qc.loc[~qc["passed"], "check"].tolist()
        raise AssertionError(f"Step 33 QC failed: {failed}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument(
        "--traits", type=Path,
        default=Path(__file__).resolve().parents[1] / "metadata/traits.tsv",
    )
    args = parser.parse_args()
    _stage_primary(args)
    _stage_decomposition(args)

    global PROJECT_ROOT, RESULTS_DIR, PLEIOTROPY_DIR, DRIVER_PATH, OUT_DIR
    PROJECT_ROOT = args.project_root.expanduser().resolve()
    RESULTS_DIR = PROJECT_ROOT / "results"
    PLEIOTROPY_DIR = RESULTS_DIR / "sv_pleiotropy"
    DRIVER_PATH = (
        PLEIOTROPY_DIR
        / "driver_decomposition/lead_sv_driver_decomposition.non_ratio_primary.csv"
    )
    OUT_DIR = PLEIOTROPY_DIR / "lead_criticality_locality_audit"
    _stage_criticality(args)


if __name__ == "__main__":
    main()
