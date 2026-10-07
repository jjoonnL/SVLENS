"""Assign the prespecified QC-caution flag to recurrent exact lead SVs."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.expanduser().resolve()
    output = root / "results/sv_pleiotropy"
    master = pd.read_csv(output / "sv_pleiotropy_master.csv")
    recurrent_all = set(
        master.groupby("lead_sv_id")["trait"].nunique().loc[lambda x: x.ge(2)].index
    )
    hits = master.loc[master["lead_sv_id"].isin(recurrent_all)]
    raw = scan_raw_leads(hits, root / "SV_association/NFE/qtbig.A")
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


if __name__ == "__main__":
    main()
