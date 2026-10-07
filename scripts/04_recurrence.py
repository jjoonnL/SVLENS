"""Audit recurrent exact lead SVs and test phenotype structure."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


# Recurrent-lead-SV QC caution
SUFFIX = "_adjAgeSexYobPC_InvNorm.txt.gz"


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def scan_raw_leads(master, raw_dir):
    records = []
    for trait, hits in master.groupby("trait"):
        target_ids = set(hits["lead_sv_id"])
        path = raw_dir / f"{trait}{SUFFIX}"
        found = set()
        with gzip.open(path, "rt") as handle:
            columns = handle.readline().rstrip("\n").split("\t")
            index = {name: position for position, name in enumerate(columns)}
            required = {"Name", "effectAlleleFreq", "Beta", "pval", "SE", "N", "info"}
            if missing := required - set(index):
                raise ValueError(f"{path} is missing columns: {sorted(missing)}")
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                sv_id = fields[index["Name"]]
                if sv_id not in target_ids:
                    continue
                found.add(sv_id)
                eaf = to_float(fields[index["effectAlleleFreq"]])
                n = to_float(fields[index["N"]])
                maf = min(eaf, 1 - eaf) if pd.notna(eaf) else np.nan
                records.append({
                    "trait": trait, "lead_sv_id": sv_id,
                    "raw_beta": to_float(fields[index["Beta"]]),
                    "raw_pval": to_float(fields[index["pval"]]),
                    "raw_se": to_float(fields[index["SE"]]),
                    "raw_N": n,
                    "raw_info": to_float(fields[index["info"]]),
                    "raw_maf": maf,
                    "approx_mac": 2 * n * maf if pd.notna(n) and pd.notna(maf) else np.nan,
                })
        if missing := target_ids - found:
            raise ValueError(f"{trait}: {len(missing)} lead SVs missing from raw summary")
    return pd.DataFrame(records)


def assign_qc_caution(master, raw_records):
    detail = master.merge(
        raw_records, on=["trait", "lead_sv_id"], how="left", validate="many_to_one"
    )
    if detail[["raw_N", "raw_info", "approx_mac"]].isna().any(axis=None):
        raise ValueError("Missing raw QC measurements for at least one lead-SV trait")
    summary = (
        detail.groupby("lead_sv_id")
        .agg(
            genes=("gene_name", lambda x: "; ".join(sorted(set(x)))),
            n_traits=("trait", "nunique"),
            lead_sv_chrom=("lead_sv_chrom", "first"),
            lead_sv_pos=("lead_sv_pos", "first"),
            min_raw_pval=("raw_pval", "min"),
            min_approx_mac=("approx_mac", "min"),
            median_approx_mac=("approx_mac", "median"),
            min_info=("raw_info", "min"),
        )
        .reset_index()
    )
    summary["is_hla_region"] = (
        summary["lead_sv_chrom"].eq("chr6")
        & summary["lead_sv_pos"].between(28_000_000, 34_000_000)
    )
    summary["has_olfactory_gene"] = summary["genes"].map(
        lambda value: bool(re.search(r"(?:^|; )OR\d", str(value)))
    )
    summary["low_mac_flag"] = summary["min_approx_mac"].lt(20)
    summary["low_info_flag"] = summary["min_info"].lt(0.8)
    reasons = ["is_hla_region", "has_olfactory_gene", "low_mac_flag", "low_info_flag"]
    summary["qc_caution_reason"] = summary.apply(
        lambda row: ";".join(name for name in reasons if row[name]), axis=1
    )
    summary["case_tier"] = np.where(
        summary[reasons].any(axis=1), "QC_caution", "eligible"
    )
    return summary


def _stage_qc(args):
    root = args.project_root.expanduser().resolve()
    output = root / "results/sv_pleiotropy"
    master = pd.read_csv(output / "sv_pleiotropy_master.csv")
    recurrent_all = set(
        master.groupby("lead_sv_id")["trait"].nunique().loc[lambda x: x.ge(2)].index
    )
    hits = master.loc[master["lead_sv_id"].isin(recurrent_all)]
    raw = scan_raw_leads(hits, root / "SV_association")
    qc = assign_qc_caution(hits, raw)
    qc.to_csv(output / "recurrent_lead_sv_case_priority.csv", index=False)
    lead_path = output / "lead_sv_level_pleiotropy_summary.non_ratio_primary.csv"
    lead = pd.read_csv(lead_path)
    fields = ["lead_sv_id", "case_tier", "qc_caution_reason", "median_approx_mac", "min_approx_mac", "min_info"]
    lead = lead.merge(qc[fields], on="lead_sv_id", how="left", validate="one_to_one")
    lead.to_csv(lead_path, index=False)
    lead.loc[lead["is_recurrent_non_ratio"].astype(bool)].to_csv(
        output / "candidate_recurrent_lead_svs.non_ratio_primary.csv", index=False
    )
    print(f"{int(qc['case_tier'].eq('QC_caution').sum())} QC-caution lead SVs")


# Trait-category and phenotype-family enrichment
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


def _stage_enrichment(args):
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


# Target ambiguity and origin audit
IMMUNE_LOCI = [
    ("IGH", "14", 105_586_937, 106_879_844),
    ("IGK", "2", 88_857_361, 90_235_368),
    ("IGL", "22", 22_026_076, 22_922_913),
    ("TRA", "14", 21_621_904, 22_552_132),
    ("TRB", "7", 142_299_011, 142_813_287),
    ("TRG", "7", 38_240_024, 38_368_055),
    ("TRD", "14", 21_924_138, 22_469_614),
]


CALR_TYPE1 = ("19", 12_943_750, 12_943_802, "DEL", 52)


def target_assignment(master):
    target = (
        master.groupby("lead_sv_id", dropna=False)
        .agg(
            n_genes=("gene_name", "nunique"), genes=("gene_name", join_unique),
            n_traits=("trait", "nunique"), traits=("trait", join_unique),
            n_hits=("trait", "size"),
            categories=("trait_category", join_unique),
            broad_modules=("broad_module", join_unique),
            min_gene_p=("p_acat_o", "min"),
            min_lead_sv_p=("lead_sv_p", "min"),
            lead_sv_chrom=("lead_sv_chrom", "first"),
            lead_sv_pos=("lead_sv_pos", "first"),
            lead_sv_start=("lead_sv_start", "first"),
            lead_sv_end=("lead_sv_end", "first"),
            lead_sv_type=("lead_sv_type", "first"),
            lead_sv_layer=("lead_sv_layer", "first"),
            lead_sv_maf=("lead_sv_maf", "first"),
            lead_sv_size=("lead_sv_size", "first"),
            support_classes=("support_class_final", join_unique),
        )
        .reset_index()
    )
    target["is_recurrent"] = target["n_traits"].ge(2)
    target["target_ambiguity_class"] = np.where(
        target["n_genes"].eq(1), "single_gene_target", "ambiguous_multi_gene_target"
    )
    return target


def origin_audit(lead, category_assignment):
    parsed = lead["lead_sv_id"].str.extract(
        r"^chr(?P<chrom>[^:]+):(?P<start>\d+):<(?P<sv_type>[^:]+):SVSIZE=(?P<size>\d+)"
    )
    if parsed.isna().any(axis=None):
        raise ValueError("At least one exact lead-SV ID cannot be parsed")
    annotated = lead.copy()
    annotated["chrom_norm"] = parsed["chrom"].str.replace("chr", "", regex=False)
    annotated["sv_start_audit"] = parsed["start"].astype(int)
    annotated["sv_type_audit"] = parsed["sv_type"]
    annotated["sv_size_audit"] = parsed["size"].astype(int)
    annotated["sv_end_audit"] = annotated["sv_start_audit"] + np.where(
        annotated["sv_type_audit"].eq("INS"), 1, annotated["sv_size_audit"]
    )
    loci = pd.DataFrame(
        IMMUNE_LOCI, columns=["immune_locus", "chrom", "locus_start", "locus_end"]
    )

    def immune_overlap(row):
        same_chrom = loci.loc[loci["chrom"].eq(row["chrom_norm"])]
        hits = same_chrom.loc[
            same_chrom["locus_end"].ge(row["sv_start_audit"])
            & same_chrom["locus_start"].le(row["sv_end_audit"])
        ]
        if hits.empty:
            return pd.Series({
                "immune_loci": "", "max_immune_overlap_bp": 0,
                "within_immune_locus": False,
            })
        overlap = np.minimum(hits["locus_end"], row["sv_end_audit"]) - np.maximum(
            hits["locus_start"], row["sv_start_audit"]
        ) + 1
        within = (
            hits["locus_start"].le(row["sv_start_audit"])
            & hits["locus_end"].ge(row["sv_end_audit"])
        ).any()
        return pd.Series({
            "immune_loci": ";".join(hits["immune_locus"]),
            "max_immune_overlap_bp": int(overlap.max()),
            "within_immune_locus": bool(within),
        })

    annotated = pd.concat([annotated, annotated.apply(immune_overlap, axis=1)], axis=1)
    chrom, start, end, sv_type, size = CALR_TYPE1
    calr = (
        annotated["chrom_norm"].eq(chrom)
        & annotated["sv_start_audit"].eq(start)
        & annotated["sv_end_audit"].eq(end)
        & annotated["sv_type_audit"].eq(sv_type)
        & annotated["sv_size_audit"].eq(size)
    )
    annotated["known_somatic_label"] = np.where(calr, "CALR_type1_52bp_del", "")
    annotated["origin_class"] = np.select(
        [calr, annotated["immune_loci"].ne("")],
        ["known_somatic_driver", "immune_receptor_locus_overlap"],
        default="other",
    )
    annotated["blood_derived_candidate"] = annotated["origin_class"].ne("other")
    assignment = category_assignment[[
        "lead_sv_id", "n_significant_categories", "plot_group", "significant_categories",
    ]].copy()
    assignment["is_category_enriched"] = assignment["n_significant_categories"].gt(0)
    return annotated.merge(
        assignment.drop(columns="n_significant_categories"), on="lead_sv_id", how="left",
        validate="one_to_one",
    )


def exclusion_sensitivity(master, target, origin, category_pairs):
    excluded = set(origin.loc[origin["blood_derived_candidate"], "lead_sv_id"])
    filtered = category_pairs.loc[~category_pairs["lead_sv_id"].isin(excluded)].copy()
    p = filtered["hypergeom_p"].to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    adjusted = np.empty(len(p))
    adjusted[order] = np.clip(q, 0, 1)
    filtered["category_fdr_filtered"] = adjusted
    filtered["category_fdr_filtered_le_0_05"] = filtered["category_fdr_filtered"].le(0.05)
    remaining = set(origin["lead_sv_id"]) - excluded
    original_enriched = set(category_pairs.loc[
        category_pairs["category_fdr_lt_0_05"], "lead_sv_id"
    ])
    filtered_enriched = set(filtered.loc[
        filtered["category_fdr_filtered_le_0_05"], "lead_sv_id"
    ])

    def metrics(label, lead_ids, enriched_ids):
        lead = origin.loc[origin["lead_sv_id"].isin(lead_ids)]
        hits = master.loc[master["lead_sv_id"].isin(lead_ids)]
        recurrent = set(lead.loc[lead["is_recurrent_non_ratio"], "lead_sv_id"])
        recurrent_target = target.loc[target["lead_sv_id"].isin(recurrent)]
        return {
            "analysis_set": label,
            "n_gene_trait_associations": len(hits),
            "n_traits": hits["trait"].nunique(),
            "n_genes": hits["gene_name"].nunique(),
            "n_lead_svs": len(lead),
            "n_recurrent_lead_svs": len(recurrent),
            "n_category_enriched_recurrent_lead_svs": len(enriched_ids & recurrent),
            "n_single_gene_recurrent_lead_svs": int(
                recurrent_target["target_ambiguity_class"].eq("single_gene_target").sum()
            ),
        }

    summary = pd.DataFrame([
        metrics("Original", set(origin["lead_sv_id"]), original_enriched),
        metrics("Exclude known somatic and IG/TR overlap", remaining, filtered_enriched),
    ])
    return filtered, summary


def _stage_locus_audit(args):
    base = args.project_root / "results/sv_pleiotropy"
    master = pd.read_csv(base / "sv_pleiotropy_master.non_ratio_primary.csv")
    lead = pd.read_csv(base / "lead_sv_level_pleiotropy_summary.non_ratio_primary.csv")
    assignment = pd.read_csv(
        base / "category_enrichment/recurrent_lead_sv_category_assignment.csv"
    )
    category_pairs = pd.read_csv(
        base / "category_enrichment/recurrent_lead_sv_category_enrichment.csv"
    )
    target = target_assignment(master)
    origin = origin_audit(lead, assignment)
    filtered_pairs, sensitivity = exclusion_sensitivity(master, target, origin, category_pairs)
    target_out = base / "target_ambiguity"
    origin_out = base / "somatic_immune_audit"
    target_out.mkdir(parents=True, exist_ok=True)
    origin_out.mkdir(parents=True, exist_ok=True)
    target.to_csv(target_out / "lead_sv_target_ambiguity.all_non_ratio.csv", index=False)
    target.loc[target["is_recurrent"]].to_csv(
        target_out / "recurrent_lead_sv_target_ambiguity.csv", index=False
    )
    origin.to_csv(origin_out / "lead_sv_origin_audit.csv", index=False)
    filtered_pairs.to_csv(
        origin_out / "category_enrichment_excluding_blood_derived_candidates.csv", index=False
    )
    sensitivity.to_csv(origin_out / "core_architecture_exclusion_sensitivity.csv", index=False)
    print(f"{len(target)} target rows; {int(origin['blood_derived_candidate'].sum())} flagged lead SVs")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument(
        "--traits", type=Path,
        default=Path(__file__).resolve().parents[1] / "metadata/traits.tsv",
    )
    args = parser.parse_args()
    _stage_qc(args)
    _stage_enrichment(args)
    _stage_locus_audit(args)


if __name__ == "__main__":
    main()
