"""ACAT functions used by the manuscript's SV gene-level analysis."""

import numpy as np
import pandas as pd


SV_TYPES = ["DEL", "DUP", "INS"]
P_SMALL = 1e-15


# ACAT P-value combination: Liu et al., Am J Hum Genet 2019, 104:410-421.
# doi:10.1016/j.ajhg.2019.01.002
def acat_pvalue(pvals, weights=None):
    """Combine P values by ACAT, using the small-P Cauchy approximation."""
    pvals = np.asarray(pvals, dtype=float)
    pvals = np.clip(pvals, 1e-300, 1 - 1e-15)
    cauchy = np.where(
        pvals < P_SMALL,
        1.0 / (pvals * np.pi),
        np.tan((0.5 - pvals) * np.pi),
    )
    if weights is None:
        statistic = cauchy.mean()
    else:
        weights = np.asarray(weights, dtype=float)
        statistic = np.sum(weights * cauchy) / np.sum(weights)
    return float(np.clip(0.5 - np.arctan(statistic) / np.pi, 0, 1))


# Beta(1,25) density evaluated at MAF: Wu et al., Am J Hum Genet 2011, 89:82-93.
# doi:10.1016/j.ajhg.2011.05.029 (weighting precedent, this test uses ACAT).
def run_acat_gene(sv, min_sv=1):
    """Run MAF-weighted ACAT within SV types, then equal-weight ACAT across types.

    ``w_final`` is the Beta(1,25) density evaluated at MAF in the primary analysis.
    The historical ``p_acat_o`` output name is retained for downstream schema
    compatibility. This procedure is not the burden/SKAT/ACAT-V omnibus ACAT-O.
    """
    sv_valid = sv[sv["w_final"].notna()].copy()
    stratum_records = []
    for (gene_id, gene_name, sv_type), group in sv_valid.groupby(
        ["gene_id", "gene_name", "sv_type"], sort=False
    ):
        if len(group) < min_sv:
            continue
        stratum_records.append({
            "gene_id": gene_id,
            "gene_name": gene_name,
            "sv_type": sv_type,
            "n_sv": len(group),
            "p_acat": acat_pvalue(group["pval"].values, group["w_final"].values),
        })
    acat_stratum = pd.DataFrame(stratum_records)
    if acat_stratum.empty:
        return acat_stratum, pd.DataFrame()

    gene_records = []
    for (gene_id, gene_name), group in acat_stratum.groupby(
        ["gene_id", "gene_name"], sort=False
    ):
        gene_records.append({
            "gene_id": gene_id,
            "gene_name": gene_name,
            "n_strata": len(group),
            "n_sv_total": group["n_sv"].sum(),
            "p_acat_o": acat_pvalue(group["p_acat"].values),
        })
    acat_gene = pd.DataFrame(gene_records).sort_values("p_acat_o")
    return acat_stratum, acat_gene
