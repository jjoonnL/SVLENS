"""Audit significant-gene assignment and blood-derived lead-SV candidates."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# IMGT GRCh38 receptor-locus coordinates used by the manuscript analysis.
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


def join_unique(values):
    return "; ".join(sorted(pd.Series(values).dropna().astype(str).unique()))


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()
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


if __name__ == "__main__":
    main()
