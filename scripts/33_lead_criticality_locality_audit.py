"""Audit lead-SV criticality and residual-SV locality in the primary results.

Part 1 performs exhaustive leave-one-SV-out gene-level ACAT recalculation for
all 404 significant gene-trait associations. It compares loss of significance
after removal of the predefined lead SV with the loss expected after removing
a uniformly selected member of the same SV set.

Part 2 is restricted to the 57 multi-SV-supported associations. It compares
the distance and SV-type concordance between the lead SV and the strongest
residual SV with a conditional null obtained by selecting another nonlead SV
from the same gene-trait set. This preserves the observed gene, lead position,
SV density, candidate intervals, and candidate SV types. Neither analysis
establishes independence or causality, which would require individual-level
conditional or haplotype data.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from svarch.acat import SV_TYPES, acat_pvalue


if not os.environ.get("SVARCH_PROJECT_ROOT"):
    raise SystemExit("Set SVARCH_PROJECT_ROOT to the private analysis directory")
PROJECT_ROOT = Path(os.environ["SVARCH_PROJECT_ROOT"]).expanduser().resolve()
RESULTS_DIR = PROJECT_ROOT / "results"
PLEIOTROPY_DIR = RESULTS_DIR / "sv_pleiotropy"
DRIVER_PATH = (
    PLEIOTROPY_DIR
    / "driver_decomposition/lead_sv_driver_decomposition.non_ratio_primary.csv"
)
OUT_DIR = PLEIOTROPY_DIR / "lead_criticality_locality_audit"

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


def main() -> None:
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


if __name__ == "__main__":
    main()
