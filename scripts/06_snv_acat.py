#!/usr/bin/env python3
"""Uniform rare-SNV gene-level ACAT for the full 77-trait SV panel.

The exact NFE trait-to-GCST mapping was audited against the first sheet of
``GWAS catalog GCST list.xlsx``. Traits without an exact phenotype-matched
summary remain unavailable; related or derived phenotypes are never silently
substituted. Ready traits are recomputed uniformly using MAF <1%, MAC >=10,
Beta(MAF; 1,25) weights, GENCODE v49 gene-body membership, and separate
overall/functional/intronic ACAT outputs. The primary input is restricted to
biallelic SNVs. A sensitivity analysis adds sequence-resolved small indels with
an allele-length difference <50 bp; variants meeting the conventional >=50-bp
SV boundary are excluded from both analyses.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd


TRAIT_TO_GCST: dict[str, str | None] = {
    "Alanine_aminotransferase": "GCST90474304",
    "Albumin": "GCST90474294",
    "Alkaline_phosphatase": "GCST90474299",
    "Apolipoprotein_A": "GCST90474309",
    "Apolipoprotein_B": "GCST90474314",
    "Aspartate_aminotransferase": "GCST90474319",
    "Basophill_count": "GCST90474541",
    "Basophill_percentage": "GCST90474566",
    "Body_mass_index": "GCST90474606",
    "C_reactive_protein": "GCST90474349",
    "Calcium": "GCST90474334",
    "Cholesterol_total": "GCST90474339",
    "Creatinine": "GCST90474344",
    "Cystatin_C": "GCST90474354",
    "Diastolic_blood_pressure": "GCST90474631",
    "Direct_bilirubin": "GCST90474324",
    "Eosinophill_count": "GCST90474536",
    "Eosinophill_percentage": "GCST90474561",
    "Gamma_glutamyltransferase": "GCST90474359",
    "Glucose": "GCST90474364",
    "Glycated_haemoglobin_HbA1c": "GCST90474369",
    "Haematocrit_percentage": "GCST90474476",
    "Haemoglobin_concentration": "GCST90474471",
    "High_density_lipoprotein": "GCST90474374",
    "High_light_scatter_reticulocyte_count": "GCST90474601",
    "High_light_scatter_reticulocyte_percentage": "GCST90474596",
    "Hip_circumference": "GCST90474616",
    "IGF_1": "GCST90474379",
    "Immature_reticulocyte_fraction": "GCST90474591",
    "Low_density_lipoprotein_calculated": None,
    "Low_density_lipoprotein_measured": "GCST90474384",
    "Lymphocyte_count": "GCST90474521",
    "Lymphocyte_percentage": "GCST90474546",
    "Mean_arterial_pressure": None,
    "Mean_corpuscular_haemoglobin": "GCST90474486",
    "Mean_corpuscular_haemoglobin_concentration": "GCST90474491",
    "Mean_corpuscular_volume": "GCST90474481",
    "Mean_platelet_thrombocyte_volume": "GCST90474511",
    "Mean_reticulocyte_volume": "GCST90474581",
    "Mean_sphered_cell_volume": "GCST90474586",
    "Monocyte_count": "GCST90474526",
    "Monocyte_percentage": "GCST90474551",
    "Monocyte_to_lymphocyte_ratio": None,
    "Neutrophil_to_lymphocyte_ratio": None,
    "Neutrophill_count": "GCST90474531",
    "Neutrophill_percentage": "GCST90474556",
    "Non_high_density_lipoprotein": None,
    "Non_high_density_lipoprotein_Apolipoprotein_B_ratio": None,
    "Oestradiol": "GCST90474397",
    "Phosphate": "GCST90474402",
    "Platelet_count": "GCST90474501",
    "Platelet_crit": "GCST90474506",
    "Platelet_distribution_width": "GCST90474516",
    "Platelet_to_lymphocyte_ratio": None,
    "Pulse_pressure": None,
    "Pulse_rate": None,
    "Red_blood_cell_erythrocyte_count": "GCST90474466",
    "Red_blood_cell_erythrocyte_distribution_width": "GCST90474496",
    "Remnant_cholesterol": None,
    "Reticulocyte_count": "GCST90474576",
    "Reticulocyte_percentage": "GCST90474571",
    "Rheumatoid_factor": "GCST90474404",
    "Seated_height": None,
    "Sex_hormon_binding_globulin_SHBG": "GCST90474408",
    "Sitting_height": None,
    "Systolic_blood_pressure": "GCST90474636",
    "Testosterone": "GCST90474418",
    "Total_bilirubin": "GCST90474413",
    "Total_protein": "GCST90474423",
    "Triglycerides": "GCST90474428",
    "Urate": "GCST90474433",
    "Urea": "GCST90474329",
    "Vitamin_D": "GCST90474456",
    "Waist_circumference": "GCST90474611",
    "Waist_hip_ratio": None,
    "Weight": None,
    "White_blood_cell_leukocyte_count": "GCST90474461",
}

NO_EXACT_REASON = {
    "Low_density_lipoprotein_calculated": "only direct/measured LDL GWAS identified",
    "Mean_arterial_pressure": "no exact harmonized GWAS identified",
    "Monocyte_to_lymphocyte_ratio": "no exact harmonized GWAS identified",
    "Neutrophil_to_lymphocyte_ratio": "no exact harmonized GWAS identified",
    "Non_high_density_lipoprotein": "no exact harmonized GWAS identified",
    "Non_high_density_lipoprotein_Apolipoprotein_B_ratio": "no exact harmonized GWAS identified",
    "Platelet_to_lymphocyte_ratio": "no exact harmonized GWAS identified",
    "Pulse_pressure": "no exact harmonized GWAS identified",
    "Pulse_rate": "no exact harmonized GWAS identified",
    "Remnant_cholesterol": "no exact harmonized GWAS identified",
    "Seated_height": "only standing-height GWAS identified",
    "Sitting_height": "only standing-height GWAS identified",
    "Waist_hip_ratio": "no exact harmonized GWAS identified",
    "Weight": "no exact harmonized GWAS identified",
}

MAF_MAX = 0.01
MAC_MIN = 10
BONF_ALPHA = 0.05
PRIMARY_GENE_P = 2.5e-6
P_SMALL = 1e-15
VARIANT_SCHEMA = "biallelic_snv_primary_plus_lt50bp_indel_gene_body_v2"
VALID_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX"}
GWAS_COLUMNS = [
    "chromosome",
    "base_pair_location",
    "effect_allele",
    "other_allele",
    "p_value",
    "effect_allele_frequency",
    "n",
]

GENE_BY_CHROM = None
GENE_ID_TO_NAME = None
CDS_BY_GENE = None
UTR_BY_GENE = None
GENE_TABLE = None
CONFIG = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gwas-dir", type=Path, required=True)
    parser.add_argument("--gene-table", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--scratch-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--chunk-size", type=int, default=1_000_000)
    parser.add_argument("--traits", help="Optional comma-separated subset")
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-scratch-copy", action="store_true")
    return parser.parse_args()


def acat_pvalue(pvals: np.ndarray, weights: np.ndarray) -> float:
    pvals = np.clip(np.asarray(pvals, dtype=float), 1e-300, 1 - 1e-15)
    weights = np.asarray(weights, dtype=float)
    cauchy = np.where(
        pvals < P_SMALL,
        1.0 / (pvals * np.pi),
        np.tan((0.5 - pvals) * np.pi),
    )
    statistic = np.sum(weights * cauchy) / np.sum(weights)
    return float(np.clip(0.5 - np.arctan(statistic) / np.pi, 0, 1))


def maf_weights(mafs: np.ndarray) -> np.ndarray:
    return 25.0 * np.power(1.0 - np.asarray(mafs, dtype=float), 24)


def parse_gene_id(attributes: str) -> str | None:
    match = re.search(r'gene_id "([^"]+)"', str(attributes))
    return match.group(1) if match else None


def load_reference(gene_table_path: Path, gtf_path: Path) -> tuple:
    genes = pd.read_parquet(gene_table_path)
    required = {
        "chrom", "gene_id", "gene_name", "gene_start", "gene_end",
        "promoter_start", "promoter_end",
    }
    missing = required - set(genes.columns)
    if missing:
        raise ValueError(f"Gene table missing columns: {sorted(missing)}")
    genes = genes[genes["chrom"].isin(VALID_CHROMS)].copy().reset_index(drop=True)
    if genes["gene_id"].duplicated().any():
        raise ValueError("Gene table contains duplicate gene_id values")

    by_chrom = {}
    for chrom, group in genes.groupby("chrom", sort=False):
        group = group.sort_values("gene_start").reset_index(drop=True)
        by_chrom[chrom] = {
            "gene_ids": group["gene_id"].astype(str).to_numpy(),
            "starts": group["gene_start"].to_numpy(dtype=np.int64),
            "ends": group["gene_end"].to_numpy(dtype=np.int64),
            "promoter_starts": group["promoter_start"].to_numpy(dtype=np.int64),
            "promoter_ends": group["promoter_end"].to_numpy(dtype=np.int64),
        }

    names = [
        "chrom", "source", "feature", "start", "end",
        "score", "strand", "frame", "attributes",
    ]
    valid_gene_ids = set(genes["gene_id"].astype(str))
    feature_frames = []
    with pd.read_csv(
        gtf_path, sep="\t", comment="#", header=None, names=names,
        chunksize=100_000, low_memory=False,
    ) as reader:
        for chunk in reader:
            keep = chunk[chunk["feature"].isin({"CDS", "UTR"})].copy()
            if not keep.empty:
                feature_frames.append(keep)
    features = pd.concat(feature_frames, ignore_index=True)
    features["gene_id"] = features["attributes"].map(parse_gene_id)
    features = features[features["gene_id"].isin(valid_gene_ids)].copy()
    features["feature_start"] = features["start"].astype(np.int64) - 1
    features["feature_end"] = features["end"].astype(np.int64)

    intervals = {}
    for feature in ("CDS", "UTR"):
        target = {}
        for gene_id, group in features[features["feature"].eq(feature)].groupby(
            "gene_id", sort=False
        ):
            target[str(gene_id)] = list(
                zip(group["feature_start"].to_numpy(), group["feature_end"].to_numpy())
            )
        intervals[feature] = target

    gene_to_name = genes.set_index("gene_id")["gene_name"].astype(str).to_dict()
    return (
        by_chrom, gene_to_name, intervals["CDS"], intervals["UTR"], genes,
    )


def init_worker(reference: tuple, config: dict) -> None:
    global GENE_BY_CHROM, GENE_ID_TO_NAME, CDS_BY_GENE, UTR_BY_GENE
    global GENE_TABLE, CONFIG
    (
        GENE_BY_CHROM, GENE_ID_TO_NAME, CDS_BY_GENE, UTR_BY_GENE,
        GENE_TABLE,
    ) = reference
    CONFIG = config


def resolve_raw_file(gwas_dir: Path, gcst: str) -> Path:
    matches = sorted(gwas_dir.glob(f"*/{gcst}.tsv.gz"))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one raw file for {gcst}; found {len(matches)}: {matches}"
        )
    return matches[0]


def test_gzip(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as handle:
            handle.read(1)
        return True
    except (OSError, EOFError):
        return False


def normalize_chromosomes(values: pd.Series) -> np.ndarray:
    chrom = values.astype(str).str.replace(r"\.0$", "", regex=True)
    chrom = chrom.str.replace(r"^chr", "", regex=True, case=False)
    chrom = chrom.replace({"23": "X", "x": "X"})
    return ("chr" + chrom).to_numpy()


def overlaps_any(position: int, intervals: list[tuple[int, int]] | None) -> bool:
    return bool(intervals) and any(start <= position < end for start, end in intervals)


def atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_parquet(temporary, index=False)
    temporary.replace(path)


def dicts_to_cache(pvalues: dict, mafs: dict, path: Path) -> None:
    records = [
        (gene_id, pvalue, maf)
        for gene_id, values in pvalues.items()
        for pvalue, maf in zip(values, mafs[gene_id])
    ]
    atomic_parquet(pd.DataFrame(records, columns=["gene_id", "pval", "maf"]), path)


def cache_to_dicts(path: Path) -> tuple[defaultdict, defaultdict]:
    frame = pd.read_parquet(path)
    if frame.empty:
        return defaultdict(list), defaultdict(list)
    pvalues = defaultdict(list, frame.groupby("gene_id")["pval"].apply(list).to_dict())
    mafs = defaultdict(list, frame.groupby("gene_id")["maf"].apply(list).to_dict())
    return pvalues, mafs


def layer_acat(
    pvalues: dict,
    mafs: dict,
    extra_pvalues: dict | None = None,
    extra_mafs: dict | None = None,
    count_column: str = "n_variant",
) -> pd.DataFrame:
    extra_pvalues = extra_pvalues or {}
    extra_mafs = extra_mafs or {}
    records = []
    for gene_id in set(pvalues) | set(extra_pvalues):
        values = pvalues.get(gene_id, []) + extra_pvalues.get(gene_id, [])
        if not values or gene_id not in GENE_ID_TO_NAME:
            continue
        pvalue_array = np.asarray(values, dtype=float)
        maf_array = np.asarray(
            mafs.get(gene_id, []) + extra_mafs.get(gene_id, []), dtype=float
        )
        records.append({
            "gene_id": gene_id,
            "gene_name": GENE_ID_TO_NAME[gene_id],
            count_column: len(values),
            "p_acat_snp": acat_pvalue(pvalue_array, maf_weights(maf_array)),
        })
    frame = pd.DataFrame(
        records, columns=["gene_id", "gene_name", count_column, "p_acat_snp"]
    )
    threshold = BONF_ALPHA / len(frame) if len(frame) else PRIMARY_GENE_P
    frame["sig"] = frame["p_acat_snp"] < threshold
    return frame.sort_values("p_acat_snp").reset_index(drop=True)


def overall_acat(
    functional_p,
    functional_m,
    intronic_p,
    intronic_m,
    extra_functional_p=None,
    extra_functional_m=None,
    extra_intronic_p=None,
    extra_intronic_m=None,
    count_column="n_rare_variant",
) -> pd.DataFrame:
    extra_functional_p = extra_functional_p or {}
    extra_functional_m = extra_functional_m or {}
    extra_intronic_p = extra_intronic_p or {}
    extra_intronic_m = extra_intronic_m or {}
    tested_gene_ids = (
        set(functional_p) | set(intronic_p)
        | set(extra_functional_p) | set(extra_intronic_p)
    )
    records = []
    for gene_id in tested_gene_ids:
        if gene_id not in GENE_ID_TO_NAME:
            continue
        pvalues = np.asarray(
            functional_p.get(gene_id, [])
            + intronic_p.get(gene_id, [])
            + extra_functional_p.get(gene_id, [])
            + extra_intronic_p.get(gene_id, []),
            dtype=float,
        )
        mafs = np.asarray(
            functional_m.get(gene_id, [])
            + intronic_m.get(gene_id, [])
            + extra_functional_m.get(gene_id, [])
            + extra_intronic_m.get(gene_id, []),
            dtype=float,
        )
        records.append({
            "gene_id": gene_id,
            "gene_name": GENE_ID_TO_NAME[gene_id],
            count_column: len(pvalues),
            "p_acat_snp": acat_pvalue(pvalues, maf_weights(mafs)),
        })
    tested = pd.DataFrame(
        records, columns=["gene_id", "gene_name", count_column, "p_acat_snp"]
    )
    null_gene_ids = set(GENE_TABLE["gene_id"].astype(str)) - tested_gene_ids
    null = pd.DataFrame({
        "gene_id": list(null_gene_ids),
        "gene_name": [GENE_ID_TO_NAME[gene_id] for gene_id in null_gene_ids],
        count_column: 0,
        "p_acat_snp": 1.0,
    })
    frame = pd.concat([tested, null], ignore_index=True)
    frame["sig"] = frame["p_acat_snp"] < PRIMARY_GENE_P
    return frame.sort_values("p_acat_snp").reset_index(drop=True)


def output_complete(output_dir: Path, gcst: str) -> bool:
    required = [
        output_dir / "snv_acat_gene.parquet",
        output_dir / "snv_acat_functional.parquet",
        output_dir / "snv_acat_intronic.parquet",
        output_dir / "small_variant_acat_gene.parquet",
        output_dir / "small_variant_acat_functional.parquet",
        output_dir / "small_variant_acat_intronic.parquet",
        output_dir / "run_metadata.json",
    ]
    if not all(path.exists() for path in required):
        return False
    metadata = json.loads(required[-1].read_text())
    return (
        metadata.get("status") == "complete"
        and metadata.get("gcst") == gcst
        and metadata.get("variant_schema") == VARIANT_SCHEMA
        and metadata.get("gene_membership")
        == "GENCODE_v49_gene_body_0_based_half_open"
    )


def assign_variants_to_genes(
    chunk: pd.DataFrame,
    valid: np.ndarray,
    positions_all: np.ndarray,
    pvalues_all: np.ndarray,
    mafs_all: np.ndarray,
    functional_p: defaultdict,
    functional_m: defaultdict,
    intronic_p: defaultdict,
    intronic_m: defaultdict,
) -> int:
    if not valid.any():
        return 0
    chromosomes = normalize_chromosomes(chunk.loc[valid, "chromosome"])
    # GWAS positions are 1-based; reference intervals are 0-based and half-open.
    positions = positions_all[valid].astype(np.int64) - 1
    pvalues = pvalues_all[valid].astype(float)
    mafs = mafs_all[valid].astype(float)
    assignments = 0

    for chromosome in np.unique(chromosomes):
        if chromosome not in GENE_BY_CHROM:
            continue
        genes = GENE_BY_CHROM[chromosome]
        mask = chromosomes == chromosome
        chrom_positions = positions[mask]
        chrom_pvalues = pvalues[mask]
        chrom_mafs = mafs[mask]
        for start in range(0, len(chrom_positions), 5_000):
            batch_positions = chrom_positions[start : start + 5_000]
            overlaps = (
                (batch_positions[:, None] >= genes["starts"])
                & (batch_positions[:, None] < genes["ends"])
            )
            variant_indices, gene_indices = np.where(overlaps)
            assignments += len(variant_indices)
            for variant_index, gene_index in zip(variant_indices, gene_indices):
                position = int(batch_positions[variant_index])
                gene_id = genes["gene_ids"][gene_index]
                pvalue = chrom_pvalues[start + variant_index]
                maf = chrom_mafs[start + variant_index]
                functional = (
                    genes["promoter_starts"][gene_index]
                    <= position < genes["promoter_ends"][gene_index]
                )
                if not functional:
                    functional = overlaps_any(position, CDS_BY_GENE.get(gene_id))
                if not functional:
                    functional = overlaps_any(position, UTR_BY_GENE.get(gene_id))
                if functional:
                    functional_p[gene_id].append(pvalue)
                    functional_m[gene_id].append(maf)
                else:
                    intronic_p[gene_id].append(pvalue)
                    intronic_m[gene_id].append(maf)
    return assignments


def process_trait(task: tuple[str, str, str]) -> dict:
    trait, gcst, raw_file_string = task
    raw_file = Path(raw_file_string)
    output_dir = Path(CONFIG["out_dir"]) / trait
    output_dir.mkdir(parents=True, exist_ok=True)
    if output_complete(output_dir, gcst) and not CONFIG["overwrite"]:
        metadata = json.loads((output_dir / "run_metadata.json").read_text())
        metadata["status"] = "skipped_complete"
        return metadata

    started = time.time()
    pid = os.getpid()
    print(f"[PID {pid}] START {trait} ({gcst})", flush=True)
    snv_cache_functional = output_dir / "snv_cache_func.parquet"
    snv_cache_intronic = output_dir / "snv_cache_intr.parquet"
    indel_cache_functional = output_dir / "small_indel_cache_func.parquet"
    indel_cache_intronic = output_dir / "small_indel_cache_intr.parquet"
    cache_metadata = output_dir / "cache_metadata.json"
    scan_stats_path = output_dir / "scan_stats.json"
    expected_cache = {
        "gcst": gcst,
        "maf_max": MAF_MAX,
        "mac_min": MAC_MIN,
        "variant_schema": VARIANT_SCHEMA,
        "gene_membership": "GENCODE_v49_gene_body_0_based_half_open",
        "small_indel_length_difference_max_exclusive": 50,
    }
    cache_valid = (
        snv_cache_functional.exists()
        and snv_cache_intronic.exists()
        and indel_cache_functional.exists()
        and indel_cache_intronic.exists()
        and cache_metadata.exists()
        and scan_stats_path.exists()
        and json.loads(cache_metadata.read_text()) == expected_cache
    )

    stats = {
        "n_rows": 0,
        "n_frequency_eligible": 0,
        "n_rare_snv": 0,
        "n_rare_small_indel": 0,
        "n_rare_small_variant": 0,
        "n_excluded_mnv": 0,
        "n_excluded_indel_ge50bp": 0,
        "n_excluded_other_variant_class": 0,
        "n_snv_gene_assignments": 0,
        "n_small_indel_gene_assignments": 0,
    }
    if cache_valid and not CONFIG["overwrite"]:
        snv_functional_p, snv_functional_m = cache_to_dicts(snv_cache_functional)
        snv_intronic_p, snv_intronic_m = cache_to_dicts(snv_cache_intronic)
        indel_functional_p, indel_functional_m = cache_to_dicts(
            indel_cache_functional
        )
        indel_intronic_p, indel_intronic_m = cache_to_dicts(indel_cache_intronic)
        stats.update(json.loads(scan_stats_path.read_text()))
    else:
        snv_functional_p, snv_functional_m = defaultdict(list), defaultdict(list)
        snv_intronic_p, snv_intronic_m = defaultdict(list), defaultdict(list)
        indel_functional_p, indel_functional_m = defaultdict(list), defaultdict(list)
        indel_intronic_p, indel_intronic_m = defaultdict(list), defaultdict(list)
        scratch_file = Path(CONFIG["scratch_dir"]) / f"{gcst}.tsv.gz"
        scan_file = raw_file
        copied = False
        try:
            if CONFIG["copy_to_scratch"]:
                shutil.copy2(raw_file, scratch_file)
                scan_file = scratch_file
                copied = True
            reader = pd.read_csv(
                scan_file, sep="\t", usecols=GWAS_COLUMNS,
                chunksize=CONFIG["chunk_size"], low_memory=False,
            )
            for chunk_index, chunk in enumerate(reader, start=1):
                stats["n_rows"] += len(chunk)
                eaf = pd.to_numeric(chunk["effect_allele_frequency"], errors="coerce").to_numpy()
                n = pd.to_numeric(chunk["n"], errors="coerce").to_numpy()
                pvalues = pd.to_numeric(chunk["p_value"], errors="coerce").to_numpy()
                positions = pd.to_numeric(chunk["base_pair_location"], errors="coerce").to_numpy()
                mafs = np.minimum(eaf, 1 - eaf)
                macs = np.rint(2 * n * mafs)
                valid_frequency = (
                    np.isfinite(mafs) & np.isfinite(macs) & np.isfinite(pvalues)
                    & np.isfinite(positions) & (pvalues >= 0) & (pvalues <= 1)
                    & (mafs >= 0) & (mafs < MAF_MAX) & (macs >= MAC_MIN)
                )
                effect = chunk["effect_allele"].astype("string").str.strip().str.upper()
                other = chunk["other_allele"].astype("string").str.strip().str.upper()
                sequence_resolved = (
                    effect.str.fullmatch(r"[ACGT]+", na=False)
                    & other.str.fullmatch(r"[ACGT]+", na=False)
                ).to_numpy(dtype=bool)
                different = effect.ne(other).fillna(False).to_numpy(dtype=bool)
                effect_length = effect.str.len().fillna(0).to_numpy(dtype=np.int64)
                other_length = other.str.len().fillna(0).to_numpy(dtype=np.int64)
                is_snv = (
                    sequence_resolved & different
                    & (effect_length == 1) & (other_length == 1)
                )
                is_small_indel = (
                    sequence_resolved & different
                    & (effect_length != other_length)
                    & (np.abs(effect_length - other_length) < 50)
                )
                is_mnv = (
                    sequence_resolved & different
                    & (effect_length == other_length) & (effect_length > 1)
                )
                is_indel_ge50 = (
                    sequence_resolved & different
                    & (effect_length != other_length)
                    & (np.abs(effect_length - other_length) >= 50)
                )
                valid_snv = valid_frequency & is_snv
                valid_small_indel = valid_frequency & is_small_indel
                valid_mnv = valid_frequency & is_mnv
                valid_indel_ge50 = valid_frequency & is_indel_ge50
                classified = is_snv | is_small_indel | is_mnv | is_indel_ge50
                valid_other = valid_frequency & ~classified
                stats["n_frequency_eligible"] += int(valid_frequency.sum())
                stats["n_excluded_mnv"] += int(valid_mnv.sum())
                stats["n_excluded_indel_ge50bp"] += int(valid_indel_ge50.sum())
                stats["n_excluded_other_variant_class"] += int(valid_other.sum())
                if not valid_snv.any() and not valid_small_indel.any():
                    continue
                stats["n_rare_snv"] += int(valid_snv.sum())
                stats["n_rare_small_indel"] += int(valid_small_indel.sum())
                stats["n_snv_gene_assignments"] += assign_variants_to_genes(
                    chunk, valid_snv, positions, pvalues, mafs,
                    snv_functional_p, snv_functional_m,
                    snv_intronic_p, snv_intronic_m,
                )
                stats["n_small_indel_gene_assignments"] += assign_variants_to_genes(
                    chunk, valid_small_indel, positions, pvalues, mafs,
                    indel_functional_p, indel_functional_m,
                    indel_intronic_p, indel_intronic_m,
                )
                if chunk_index % 10 == 0:
                    print(
                        f"[PID {pid}] {trait}: {stats['n_rows']/1e6:.1f}M rows; "
                        f"{time.time()-started:.0f}s", flush=True,
                    )
        finally:
            if copied and scratch_file.exists():
                scratch_file.unlink()
        stats["n_rare_small_variant"] = (
            stats["n_rare_snv"] + stats["n_rare_small_indel"]
        )
        classified_total = (
            stats["n_rare_small_variant"]
            + stats["n_excluded_mnv"]
            + stats["n_excluded_indel_ge50bp"]
            + stats["n_excluded_other_variant_class"]
        )
        if classified_total != stats["n_frequency_eligible"]:
            raise AssertionError(
                "Variant-class counts do not sum to n_frequency_eligible"
            )
        dicts_to_cache(snv_functional_p, snv_functional_m, snv_cache_functional)
        dicts_to_cache(snv_intronic_p, snv_intronic_m, snv_cache_intronic)
        dicts_to_cache(
            indel_functional_p, indel_functional_m, indel_cache_functional
        )
        dicts_to_cache(indel_intronic_p, indel_intronic_m, indel_cache_intronic)
        cache_metadata.write_text(json.dumps(expected_cache, indent=2) + "\n")
        scan_stats_path.write_text(json.dumps(stats, indent=2) + "\n")

    snv_overall = overall_acat(
        snv_functional_p, snv_functional_m, snv_intronic_p, snv_intronic_m,
        count_column="n_rare_snv",
    )
    snv_functional = layer_acat(
        snv_functional_p, snv_functional_m, count_column="n_snv"
    )
    snv_intronic = layer_acat(
        snv_intronic_p, snv_intronic_m, count_column="n_snv"
    )
    small_overall = overall_acat(
        snv_functional_p, snv_functional_m, snv_intronic_p, snv_intronic_m,
        indel_functional_p, indel_functional_m,
        indel_intronic_p, indel_intronic_m,
        count_column="n_rare_small_variant",
    )
    small_functional = layer_acat(
        snv_functional_p, snv_functional_m,
        indel_functional_p, indel_functional_m,
        count_column="n_small_variant",
    )
    small_intronic = layer_acat(
        snv_intronic_p, snv_intronic_m,
        indel_intronic_p, indel_intronic_m,
        count_column="n_small_variant",
    )
    atomic_parquet(snv_overall, output_dir / "snv_acat_gene.parquet")
    atomic_parquet(snv_functional, output_dir / "snv_acat_functional.parquet")
    atomic_parquet(snv_intronic, output_dir / "snv_acat_intronic.parquet")
    atomic_parquet(small_overall, output_dir / "small_variant_acat_gene.parquet")
    atomic_parquet(
        small_functional, output_dir / "small_variant_acat_functional.parquet"
    )
    atomic_parquet(
        small_intronic, output_dir / "small_variant_acat_intronic.parquet"
    )
    metadata = {
        "trait": trait,
        "gcst": gcst,
        "raw_file": str(raw_file),
        "status": "complete",
        "variant_schema": VARIANT_SCHEMA,
        "gene_membership": "GENCODE_v49_gene_body_0_based_half_open",
        **stats,
        "n_genes_tested": int((snv_overall["n_rare_snv"] > 0).sum()),
        "n_overall_hits": int(snv_overall["sig"].sum()),
        "n_functional_hits": int(snv_functional["sig"].sum()),
        "n_intronic_hits": int(snv_intronic["sig"].sum()),
        "n_small_variant_genes_tested": int(
            (small_overall["n_rare_small_variant"] > 0).sum()
        ),
        "n_small_variant_overall_hits": int(small_overall["sig"].sum()),
        "n_small_variant_functional_hits": int(small_functional["sig"].sum()),
        "n_small_variant_intronic_hits": int(small_intronic["sig"].sum()),
        "maf_max_exclusive": MAF_MAX,
        "mac_min_inclusive": MAC_MIN,
        "primary_gene_p_threshold": PRIMARY_GENE_P,
        "snv_functional_bonferroni_threshold": (
            BONF_ALPHA / len(snv_functional) if len(snv_functional) else PRIMARY_GENE_P
        ),
        "snv_intronic_bonferroni_threshold": (
            BONF_ALPHA / len(snv_intronic) if len(snv_intronic) else PRIMARY_GENE_P
        ),
        "small_variant_functional_bonferroni_threshold": (
            BONF_ALPHA / len(small_functional)
            if len(small_functional) else PRIMARY_GENE_P
        ),
        "small_variant_intronic_bonferroni_threshold": (
            BONF_ALPHA / len(small_intronic)
            if len(small_intronic) else PRIMARY_GENE_P
        ),
        "elapsed_seconds": round(time.time() - started, 1),
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(
        f"[PID {pid}] DONE {trait}: overall hits={metadata['n_overall_hits']}; "
        f"{metadata['elapsed_seconds']:.0f}s", flush=True,
    )
    return metadata


def main() -> None:
    args = parse_args()
    if args.workers < 1 or args.chunk_size < 1:
        raise ValueError("--workers and --chunk-size must be positive")
    if len(TRAIT_TO_GCST) != 77:
        raise AssertionError(f"Expected 77 traits, found {len(TRAIT_TO_GCST)}")
    if sum(gcst is not None for gcst in TRAIT_TO_GCST.values()) != 63:
        raise AssertionError("Expected 63 exact phenotype-matched traits")

    selected = TRAIT_TO_GCST.copy()
    if args.traits:
        requested = {item.strip() for item in args.traits.split(",") if item.strip()}
        unknown = requested - set(selected)
        if unknown:
            raise ValueError(f"Unknown traits: {sorted(unknown)}")
        selected = {trait: gcst for trait, gcst in selected.items() if trait in requested}

    gwas_dir = Path(args.gwas_dir)
    gene_table_path = Path(args.gene_table)
    gtf_path = Path(args.gtf)
    out_dir = Path(args.out_dir)
    scratch_dir = Path(args.scratch_dir)
    if not gwas_dir.exists():
        raise FileNotFoundError(gwas_dir)
    for required in (gene_table_path, gtf_path):
        if not required.exists():
            raise FileNotFoundError(required)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_files = {}
    preflight_rows = []
    for trait, gcst in selected.items():
        if gcst is None:
            preflight_rows.append({
                "trait": trait,
                "gcst": None,
                "raw_file": None,
                "compressed_size_gb": np.nan,
                "gzip_ok": np.nan,
                "status": "no_exact_harmonized_gwas",
                "note": NO_EXACT_REASON[trait],
            })
            continue
        try:
            raw_file = resolve_raw_file(gwas_dir, gcst)
            gzip_ok = test_gzip(raw_file)
            status = "ready" if gzip_ok else "invalid_gzip"
            note = ""
        except RuntimeError as error:
            raw_file = None
            gzip_ok = np.nan
            status = "raw_file_resolution_error"
            note = str(error)
        if raw_file is not None:
            raw_files[trait] = raw_file
        preflight_rows.append({
            "trait": trait,
            "gcst": gcst,
            "raw_file": str(raw_file) if raw_file is not None else None,
            "compressed_size_gb": (
                round(raw_file.stat().st_size / 1e9, 3) if raw_file is not None else np.nan
            ),
            "gzip_ok": gzip_ok,
            "status": status,
            "note": note,
        })
    preflight = pd.DataFrame(preflight_rows)
    preflight.to_csv(out_dir / "step27_full_rare_snv_preflight.csv", index=False)
    print(preflight.to_string(index=False), flush=True)
    failed = preflight[~preflight["status"].isin({"ready", "no_exact_harmonized_gwas"})]
    if not failed.empty:
        raise RuntimeError("At least one mapped raw GWAS failed preflight")
    if args.preflight_only:
        return

    scratch_dir.mkdir(parents=True, exist_ok=True)
    print("Loading shared gene and CDS/UTR references...", flush=True)
    reference = load_reference(gene_table_path, gtf_path)
    config = {
        "out_dir": str(out_dir),
        "scratch_dir": str(scratch_dir),
        "chunk_size": args.chunk_size,
        "copy_to_scratch": not args.no_scratch_copy,
        "overwrite": args.overwrite,
    }
    tasks = [
        (trait, gcst, str(raw_files[trait]))
        for trait, gcst in selected.items() if gcst is not None
    ]
    with Pool(
        processes=args.workers,
        initializer=init_worker,
        initargs=(reference, config),
    ) as pool:
        results = pool.map(process_trait, tasks)

    summary = pd.DataFrame(results).sort_values("trait").reset_index(drop=True)
    summary.to_csv(out_dir / "step27_full_rare_snv_run_summary.csv", index=False)
    audit = preflight.merge(
        summary.drop(columns=["gcst", "raw_file"], errors="ignore").rename(
            columns={"status": "run_status"}
        ),
        on="trait", how="left", validate="one_to_one",
    )
    audit["final_status"] = np.where(
        audit["status"].eq("no_exact_harmonized_gwas"),
        "no_exact_harmonized_gwas", audit["run_status"],
    )
    audit.to_csv(out_dir / "step27_full_rare_snv_audit.csv", index=False)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
