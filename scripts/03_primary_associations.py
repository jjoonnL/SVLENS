"""Build the non-ratio SV association master and exact lead-SV recurrence."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument(
        "--traits", type=Path,
        default=Path(__file__).resolve().parents[1] / "metadata" / "traits.tsv",
    )
    args = parser.parse_args()
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


if __name__ == "__main__":
    main()
