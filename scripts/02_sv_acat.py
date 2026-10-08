"""
NFE Quantitative Traits Pipeline
Run Steps 2 (sv_annotation) and 3 (weight_acat) across traits

Usage:
    python scripts/02_sv_acat.py --project-root /path/to/private-data-root
    python scripts/02_sv_acat.py --project-root /path/to/private-data-root --step3-only

Options:
    --skip-done   : Skip traits with existing outputs (default: True)
    --step2-only  : Run only Step 2
    --step3-only  : Run only Step 3 (requires sv_annotated.parquet)
"""

import os, sys, re, time, argparse
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import beta as beta_dist
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from acat import run_acat_gene, SV_TYPES

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = ASSOC_DIR = TRAITS_DIR = GENE_TABLE = GTF_PATH = None
MAF_THRESHOLD = 0.05
PROMOTER_WIN  = 2000
SUMSTATS_SUFFIX = "_adjAgeSexYobPC_InvNorm"

# ── Traits ───────────────────────────────────────────────────────────────────
def set_project_root(root):
    global BASE, ASSOC_DIR, TRAITS_DIR, GENE_TABLE, GTF_PATH
    BASE = Path(root).expanduser().resolve()
    ASSOC_DIR = BASE / "SV_association"
    TRAITS_DIR = BASE / "results" / "traits"
    GENE_TABLE = BASE / "results" / "gene_table.parquet"
    GTF_PATH = BASE / "Gene_ref/gencode.v49.basic.annotation.gtf.gz"


# ════════════════════════════════════════════════════════════════════════════════
# Step 2: SV annotation
# ════════════════════════════════════════════════════════════════════════════════

def load_gene_table():
    gt = pd.read_parquet(GENE_TABLE)
    return gt


def load_feat_track():
    """Load CDS/UTR coordinates from the GTF once for all traits."""
    GTF_COLS     = ["chrom","source","feature","start","end","score","strand","frame","attributes"]
    TARGET_FEATS = {"CDS", "UTR"}
    chunks = []
    with pd.read_csv(GTF_PATH, sep="\t", comment="#", header=None,
                     names=GTF_COLS, chunksize=100_000, low_memory=False) as reader:
        for chunk in reader:
            sub = chunk[chunk["feature"].isin(TARGET_FEATS)]
            if len(sub) > 0:
                chunks.append(sub[["chrom","feature","start","end","attributes"]])
    feat_track = pd.concat(chunks, ignore_index=True)

    def parse_gene_id(attr):
        m = re.search(r'gene_id "([^"]+)"', attr)
        return m.group(1) if m else None

    feat_track["gene_id"]    = feat_track["attributes"].apply(parse_gene_id)
    feat_track["feat_start"] = feat_track["start"] - 1   # 0-based
    feat_track["feat_end"]   = feat_track["end"]
    feat_track["feat_type"]  = feat_track["feature"]
    return feat_track[["chrom","gene_id","feat_type","feat_start","feat_end"]].dropna()


def parse_sv_allele(allele):
    m = re.match(r"<(DEL|DUP|INS):SVSIZE=(\d+)\+?:", str(allele))
    if m:
        return m.group(1), int(m.group(2))
    return None, None


def find_overlaps_chrom(sv_chr, gene_chr):
    sv_s  = sv_chr["sv_start"].values
    sv_e  = sv_chr["sv_end"].values
    g_s   = gene_chr["gene_start"].values
    g_e   = gene_chr["gene_end"].values
    sv_i, g_i = np.where(
        (sv_s[:, None] < g_e[None, :]) &
        (sv_e[:, None] > g_s[None, :])
    )
    return sv_chr.index[sv_i], gene_chr.index[g_i]


def compute_feat_overlap(sv_gene_df, feat_track, flag_col):
    left = sv_gene_df[["sv_id","gene_id","sv_start","sv_end"]].copy()
    left["pair_idx"] = np.arange(len(left))
    merged = left.merge(feat_track, on=["gene_id"], how="left")
    overlap = merged[
        (merged["sv_start"] < merged["feat_end"]) &
        (merged["sv_end"]   > merged["feat_start"])
    ]["pair_idx"].unique()
    flag = np.zeros(len(sv_gene_df), dtype=bool)
    flag[overlap] = True
    sv_gene_df[flag_col] = flag
    return sv_gene_df


def run_step2(trait_name, gene_table, feat_track):
    sumstats_path = f"{ASSOC_DIR}/{trait_name}{SUMSTATS_SUFFIX}.txt.gz"
    output_dir    = f"{TRAITS_DIR}/{trait_name}"
    out_path      = f"{output_dir}/sv_annotated.parquet"
    os.makedirs(output_dir, exist_ok=True)

    # Load summary statistics
    raw = pd.read_csv(sumstats_path, sep="\t", low_memory=False)

    # Parse SV type and size
    parsed = raw["effectAllele"].apply(parse_sv_allele)
    raw["sv_type"] = [p[0] for p in parsed]
    raw["svsize"]  = [p[1] for p in parsed]
    raw = raw.dropna(subset=["sv_type"]).copy()

    # Standardize columns
    sv = raw.rename(columns={"Chrom": "chrom", "Pos": "pos"}).copy()
    sv["chrom"] = sv["chrom"].astype(str).apply(lambda x: x if x.startswith("chr") else "chr" + x)
    sv["maf"]   = sv["effectAlleleFreq"].clip(0, 1)
    sv["maf"]   = sv["maf"].apply(lambda x: min(x, 1 - x))
    sv["Beta"]  = sv["beta"] if "beta" in sv.columns else sv.get("oddsRatio", np.nan)
    sv["SE"]    = sv.get("se", np.nan)
    sv["N"]     = sv.get("n", np.nan)
    sv["info"]  = sv.get("info", np.nan)
    sv["pval"]  = sv["pval"].clip(lower=0)

    # SV coordinates (0-based half-open)
    sv["sv_start"] = sv["pos"].astype(int) - 1
    sv["sv_end"]   = sv["sv_start"] + sv["svsize"]
    sv["sv_id"]    = sv["Name"].astype(str)

    sv = sv[["sv_id","chrom","sv_start","sv_end","sv_type","svsize",
             "Beta","SE","pval","maf","N","info"]].drop_duplicates("sv_id").copy()

    # gene body overlap
    pairs = []
    for chrom in sv["chrom"].unique():
        sv_c   = sv[sv["chrom"] == chrom]
        gene_c = gene_table[gene_table["chrom"] == chrom]
        if len(sv_c) == 0 or len(gene_c) == 0:
            continue
        si, gi = find_overlaps_chrom(sv_c, gene_c)
        if len(si) == 0:
            continue
        pair = sv_c.loc[si].reset_index(drop=True).join(
            gene_c.loc[gi][["gene_id","gene_name","tss","promoter_start","promoter_end",
                            "gene_start","gene_end"]].reset_index(drop=True)
        )
        pairs.append(pair)

    if not pairs:
        print(f"  [WARN] {trait_name}: no SV-gene pairs found")
        return False

    sv_gene = pd.concat(pairs, ignore_index=True)

    # Promoter overlap
    sv_gene["is_promoter"] = (
        (sv_gene["sv_start"] < sv_gene["promoter_end"]) &
        (sv_gene["sv_end"]   > sv_gene["promoter_start"])
    )

    # CDS/UTR overlap
    feat_cds = feat_track[feat_track["feat_type"] == "CDS"]
    feat_utr = feat_track[feat_track["feat_type"] == "UTR"]
    sv_gene  = compute_feat_overlap(sv_gene, feat_cds, "is_cds")
    sv_gene  = compute_feat_overlap(sv_gene, feat_utr, "is_utr")

    # intronic
    sv_gene["is_intronic"] = ~(sv_gene["is_cds"] | sv_gene["is_utr"])
    sv_gene["n_layers"]    = (sv_gene["is_cds"].astype(int) +
                               sv_gene["is_utr"].astype(int) +
                               sv_gene["is_promoter"].astype(int) +
                               sv_gene["is_intronic"].astype(int))

    out_cols = ["sv_id","chrom","sv_start","sv_end","sv_type","svsize",
                "Beta","SE","pval","maf","N","info",
                "gene_id","gene_name","tss","promoter_start","promoter_end",
                "is_cds","is_utr","is_promoter","is_intronic","n_layers"]
    sv_gene[out_cols].to_parquet(out_path, index=False)
    return True


# ════════════════════════════════════════════════════════════════════════════════
# Step 3: Weight + ACAT
# ════════════════════════════════════════════════════════════════════════════════

def run_step3(trait_name):
    output_dir = f"{TRAITS_DIR}/{trait_name}"
    in_path    = f"{output_dir}/sv_annotated.parquet"

    if not os.path.exists(in_path):
        print(f"  [SKIP] {trait_name}: sv_annotated.parquet not found")
        return False

    sv = pd.read_parquet(in_path)
    sv = sv[sv["maf"] < MAF_THRESHOLD].copy()
    sv["w_maf"]   = beta_dist.pdf(sv["maf"], a=1, b=25)
    sv["w_final"] = sv["w_maf"]

    sv[["sv_id","chrom","sv_start","sv_end","sv_type","svsize",
        "Beta","SE","pval","maf","gene_id","gene_name",
        "is_cds","is_utr","is_promoter","is_intronic","n_layers",
        "w_maf","w_final"]].to_parquet(f"{output_dir}/sv_weighted.parquet", index=False)

    # Full ACAT
    acat_stratum, acat_gene = run_acat_gene(sv, min_sv=1)
    bonf = 0.05 / len(acat_gene)
    _, fdr, _, _ = multipletests(acat_gene["p_acat_o"].values, method="fdr_bh")
    acat_gene["fdr"] = fdr
    acat_stratum.to_parquet(f"{output_dir}/acat_stratum.parquet", index=False)
    acat_gene.to_parquet(f"{output_dir}/acat_gene.parquet", index=False)

    # Functional mask
    func_mask = sv["is_cds"] | sv["is_utr"] | sv["is_promoter"]
    intr_mask = sv["is_intronic"] & ~sv["is_cds"] & ~sv["is_utr"] & ~sv["is_promoter"]
    _, acat_gene_func     = run_acat_gene(sv[func_mask],     min_sv=1)
    _, acat_gene_intronic = run_acat_gene(sv[intr_mask],     min_sv=1)
    acat_gene_func.to_parquet(f"{output_dir}/acat_gene_functional.parquet", index=False)
    acat_gene_intronic.to_parquet(f"{output_dir}/acat_gene_intronic.parquet", index=False)

    n_hits = (acat_gene["p_acat_o"] < bonf).sum()
    return n_hits


# ════════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True,
                        help="Directory containing SV_association/, Gene_ref/, and results/")
    parser.add_argument("--skip-done",   action="store_true", default=True)
    parser.add_argument("--step2-only",  action="store_true")
    parser.add_argument("--step3-only",  action="store_true")
    parser.add_argument("--traits",      nargs="+", help="Run only the specified traits")
    args = parser.parse_args()

    set_project_root(args.project_root)
    if not ASSOC_DIR.is_dir():
        raise FileNotFoundError(ASSOC_DIR)
    all_traits = sorted(
        f.name.replace(f"{SUMSTATS_SUFFIX}.txt.gz", "")
        for f in ASSOC_DIR.glob(f"*{SUMSTATS_SUFFIX}.txt.gz")
    )
    print(f"Total quantitative traits: {len(all_traits)}")

    traits = args.traits if args.traits else all_traits

    # Load shared resources when Step 2 is needed
    gene_table = feat_track = None
    if not args.step3_only:
        print("Loading gene table...")
        gene_table = load_gene_table()
        print("Loading GTF feature track (CDS/UTR)... this may take several minutes")
        feat_track = load_feat_track()
        print(f"GTF features loaded: {len(feat_track):,}")

    t0 = time.time()
    results = []

    for i, trait in enumerate(traits, 1):
        output_dir     = f"{TRAITS_DIR}/{trait}"
        annotated_path = f"{output_dir}/sv_annotated.parquet"
        weighted_path  = f"{output_dir}/sv_weighted.parquet"

        print(f"\n[{i:3d}/{len(traits)}] {trait}")

        # Step 2
        if not args.step3_only:
            if args.skip_done and os.path.exists(annotated_path):
                print("  Step 2: SKIP (already done)")
            else:
                t1 = time.time()
                ok = run_step2(trait, gene_table, feat_track)
                print(f"  Step 2: {'OK' if ok else 'FAIL'}  ({time.time()-t1:.0f}s)")

        # Step 3
        if not args.step2_only:
            if args.skip_done and os.path.exists(weighted_path):
                print("  Step 3: SKIP (already done)")
            else:
                t1 = time.time()
                n_hits = run_step3(trait)
                print(f"  Step 3: OK  hits={n_hits}  ({time.time()-t1:.0f}s)")
                results.append({"trait": trait, "n_hits": n_hits})

    print(f"\n{'='*60}")
    print(f"Done. Total elapsed: {(time.time()-t0)/60:.1f} min")
    if results:
        df = pd.DataFrame(results).sort_values("n_hits", ascending=False)
        print(df.to_string(index=False))
