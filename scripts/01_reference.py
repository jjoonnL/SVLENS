"""Build the GENCODE v49 basic protein-coding gene-body reference."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


GTF_COLUMNS = [
    "chrom", "source", "feature", "start", "end", "score", "strand",
    "frame", "attributes",
]
GENE_COLUMNS = [
    "chrom", "gene_start", "gene_end", "strand", "gene_id", "gene_name",
    "tss", "promoter_start", "promoter_end",
]
VALID_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def attribute(text: str, name: str) -> str | None:
    match = re.search(rf'{name} "([^"]+)"', text)
    return match.group(1) if match else None


def build_gene_table(gtf_path: Path) -> pd.DataFrame:
    gene_chunks = []
    for chunk in pd.read_csv(
        gtf_path, sep="\t", comment="#", header=None, names=GTF_COLUMNS,
        chunksize=100_000, low_memory=False,
    ):
        genes = chunk.loc[chunk["feature"].eq("gene")].copy()
        if not genes.empty:
            gene_chunks.append(genes)
    if not gene_chunks:
        raise ValueError(f"No gene features found in {gtf_path}")

    genes = pd.concat(gene_chunks, ignore_index=True)
    for name in ("gene_id", "gene_name", "gene_type"):
        genes[name] = genes["attributes"].map(lambda value: attribute(value, name))
    genes = genes.loc[
        genes["gene_type"].eq("protein_coding") & genes["chrom"].isin(VALID_CHROMS)
    ].copy()
    genes["tss"] = np.where(genes["strand"].eq("+"), genes["start"] - 1, genes["end"])
    genes["promoter_start"] = (genes["tss"] - 2_000).clip(lower=0)
    genes["promoter_end"] = genes["tss"] + 2_000
    genes["gene_start"] = genes["start"] - 1
    genes["gene_end"] = genes["end"]
    return genes[GENE_COLUMNS].reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True, help="GENCODE v49 basic GTF")
    parser.add_argument("--out", type=Path, required=True, help="Output gene_table.parquet")
    args = parser.parse_args()
    table = build_gene_table(args.gtf)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(args.out, index=False)
    print(f"Wrote {len(table)} genes to {args.out}")


if __name__ == "__main__":
    main()
