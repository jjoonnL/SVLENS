"""Gene-body-matched rare-SNV convergence and recurrent-gene-preserving null."""

from __future__ import annotations

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()

    # Inputs, fixed analysis parameters and output location
    import json
    import numpy as np
    import pandas as pd

    PROJECT_ROOT = args.project_root.expanduser().resolve()
    SV_MASTER_PATH = PROJECT_ROOT / 'results/sv_pleiotropy/sv_pleiotropy_master.csv'
    GENE_TABLE_PATH = PROJECT_ROOT / 'results/gene_table.parquet'
    SNV_ROOT = PROJECT_ROOT / 'results_server/step27_full_rare_snv_gene_body_all_traits'
    AUDIT_PATH = SNV_ROOT / 'step27_full_rare_snv_audit.csv'
    OUTPUT_ROOT = PROJECT_ROOT / 'results/sv_pleiotropy/full_sv_rare_snv_convergence_gene_body'
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    POOL_SIZES = [50, 100, 200, 500]
    PRIMARY_POOL_SIZE = 200
    N_PERMUTATIONS = 100_000
    RANDOM_SEED = 20260720
    SV_PRIMARY_P = 2.5e-6
    SNV_PRIMARY_P = 2.5e-6
    OUTCOMES = ['overall_sig', 'functional_sig', 'intronic_sig', 'any_layer_sig', 'any_support_sig']
    VARIANT_SETS = ['snv', 'small_variant']
    ALL_OUTCOME_COLUMNS = OUTCOMES + [f'small_variant_{x}' for x in OUTCOMES]

    print('Project:', PROJECT_ROOT)
    print('Corrected SNV results:', SNV_ROOT)
    print(f'Permutations: {N_PERMUTATIONS:,}')

    # Analysis block 3
    audit = pd.read_csv(AUDIT_PATH)
    assert len(audit) == 77, f'Expected 77 audit rows, found {len(audit)}'
    accepted = {'complete', 'skipped_complete'}
    complete = audit[audit['final_status'].isin(accepted)].copy()
    no_exact = audit[audit['final_status'].eq('no_exact_harmonized_gwas')].copy()
    assert len(complete) == 63, f'Expected 63 exact completed traits, found {len(complete)}'
    assert len(no_exact) == 14, f'Expected 14 no-exact traits, found {len(no_exact)}'
    assert set(audit['final_status']).issubset(accepted | {'no_exact_harmonized_gwas'})

    required_trait_files = [
        'snv_acat_gene.parquet', 'snv_acat_functional.parquet', 'snv_acat_intronic.parquet',
        'small_variant_acat_gene.parquet', 'small_variant_acat_functional.parquet',
        'small_variant_acat_intronic.parquet', 'run_metadata.json', 'scan_stats.json',
    ]
    for row in complete.itertuples(index=False):
        trait_root = SNV_ROOT / row.trait
        for filename in required_trait_files:
            assert (trait_root / filename).exists(), f'Missing {row.trait}/{filename}'
        metadata = json.loads((trait_root / 'run_metadata.json').read_text())
        stats = json.loads((trait_root / 'scan_stats.json').read_text())
        assert metadata['gcst'] == row.gcst
        assert metadata['variant_schema'] == 'biallelic_snv_primary_plus_lt50bp_indel_gene_body_v2'
        assert metadata['gene_membership'] == 'GENCODE_v49_gene_body_0_based_half_open'
        assert metadata['maf_max_exclusive'] == 0.01
        assert metadata['mac_min_inclusive'] == 10
        assert metadata['primary_gene_p_threshold'] == SNV_PRIMARY_P
        for key in [
            'snv_functional_bonferroni_threshold', 'snv_intronic_bonferroni_threshold',
            'small_variant_functional_bonferroni_threshold',
            'small_variant_intronic_bonferroni_threshold',
        ]:
            assert 0 < float(metadata[key]) < 0.05, (row.trait, key, metadata[key])
        assert stats['n_rare_small_variant'] == stats['n_rare_snv'] + stats['n_rare_small_indel']
        assert stats['n_frequency_eligible'] == (
            stats['n_rare_small_variant'] + stats['n_excluded_mnv']
            + stats['n_excluded_indel_ge50bp'] + stats['n_excluded_other_variant_class']
        )

    sv_master = pd.read_csv(SV_MASTER_PATH)
    sv_master = sv_master.loc[sv_master['p_acat_o'] < SV_PRIMARY_P].copy()
    assert len(sv_master) == 427, f'Expected 427 fixed-threshold SV hits, found {len(sv_master)}'
    assert not sv_master.duplicated(['trait', 'gene_id']).any()
    evaluable_traits = sorted(set(complete['trait']) & set(sv_master['trait']))
    hit_master = sv_master[sv_master['trait'].isin(evaluable_traits)].copy()
    assert not hit_master.empty

    universe = pd.DataFrame([{
        'sv_traits_total': sv_master['trait'].nunique(),
        'exact_snv_traits': len(complete),
        'no_exact_snv_traits': len(no_exact),
        'evaluable_sv_hit_traits': hit_master['trait'].nunique(),
        'evaluable_sv_gene_trait_hits': len(hit_master),
        'evaluable_sv_hit_genes': hit_master['gene_id'].nunique(),
        'evaluable_exact_lead_svs': hit_master['lead_sv_id'].nunique(),
    }])
    universe.to_csv(OUTPUT_ROOT / 'step27_analysis_universe.csv', index=False)
    no_exact[['trait', 'note']].to_csv(OUTPUT_ROOT / 'step27_no_exact_gwas_traits.csv', index=False)

    # Analysis block 5
    gene_table = pd.read_parquet(GENE_TABLE_PATH)[['gene_id', 'gene_name', 'gene_start', 'gene_end']].copy()
    gene_table['membership_length_bp'] = gene_table['gene_end'] - gene_table['gene_start']
    gene_lengths = gene_table[['gene_id', 'membership_length_bp']]

    def read_variant_layer(trait, filename, prefix, count_column):
        frame = pd.read_parquet(SNV_ROOT / trait / filename)
        frame = frame[['gene_id', count_column, 'p_acat_snp', 'sig']].copy()
        return frame.rename(columns={
            count_column: f'{prefix}_n_variants',
            'p_acat_snp': f'{prefix}_p',
            'sig': f'{prefix}_sig',
        })

    background_frames = []
    for trait in evaluable_traits:
        sv = pd.read_parquet(PROJECT_ROOT / 'results' / trait / 'acat_gene.parquet')
        sv = sv[['gene_id', 'n_strata', 'n_sv_total', 'p_acat_o']].copy()
        overall = read_variant_layer(trait, 'snv_acat_gene.parquet', 'overall', 'n_rare_snv')
        functional = read_variant_layer(trait, 'snv_acat_functional.parquet', 'functional', 'n_snv')
        intronic = read_variant_layer(trait, 'snv_acat_intronic.parquet', 'intronic', 'n_snv')
        small_overall = read_variant_layer(
            trait, 'small_variant_acat_gene.parquet', 'small_variant_overall',
            'n_rare_small_variant',
        )
        small_functional = read_variant_layer(
            trait, 'small_variant_acat_functional.parquet', 'small_variant_functional',
            'n_small_variant',
        )
        small_intronic = read_variant_layer(
            trait, 'small_variant_acat_intronic.parquet', 'small_variant_intronic',
            'n_small_variant',
        )
        joint = sv.merge(overall, on='gene_id', how='inner', validate='one_to_one')
        joint = joint.merge(functional, on='gene_id', how='left', validate='one_to_one')
        joint = joint.merge(intronic, on='gene_id', how='left', validate='one_to_one')
        joint = joint.merge(small_overall, on='gene_id', how='inner', validate='one_to_one')
        joint = joint.merge(small_functional, on='gene_id', how='left', validate='one_to_one')
        joint = joint.merge(small_intronic, on='gene_id', how='left', validate='one_to_one')
        for layer in ['functional', 'intronic', 'small_variant_functional', 'small_variant_intronic']:
            joint[f'{layer}_n_variants'] = joint[f'{layer}_n_variants'].fillna(0).astype(int)
            joint[f'{layer}_p'] = joint[f'{layer}_p'].fillna(1.0)
            joint[f'{layer}_sig'] = joint[f'{layer}_sig'].eq(True)
        joint['overall_sig'] = joint['overall_p'] < SNV_PRIMARY_P
        joint['any_layer_sig'] = joint['functional_sig'] | joint['intronic_sig']
        joint['any_support_sig'] = joint['overall_sig'] | joint['any_layer_sig']
        joint['small_variant_overall_sig'] = joint['small_variant_overall_p'] < SNV_PRIMARY_P
        joint['small_variant_any_layer_sig'] = (
            joint['small_variant_functional_sig'] | joint['small_variant_intronic_sig']
        )
        joint['small_variant_any_support_sig'] = (
            joint['small_variant_overall_sig'] | joint['small_variant_any_layer_sig']
        )
        joint['trait'] = trait
        joint['overall_percentile'] = joint['overall_p'].rank(method='average', pct=True)
        joint['small_variant_overall_percentile'] = (
            joint['small_variant_overall_p'].rank(method='average', pct=True)
        )
        background_frames.append(joint)

    background = pd.concat(background_frames, ignore_index=True)
    background = background.merge(gene_lengths, on='gene_id', how='left', validate='many_to_one')
    assert background['membership_length_bp'].notna().all()
    assert not background.duplicated(['trait', 'gene_id']).any()
    background['log_membership_length'] = np.log1p(background['membership_length_bp'])
    background['log_overall_n_snv'] = np.log1p(background['overall_n_variants'])
    background['log_functional_n_snv'] = np.log1p(background['functional_n_variants'])
    background['log_intronic_n_snv'] = np.log1p(background['intronic_n_variants'])
    background['log_n_sv_total'] = np.log1p(background['n_sv_total'])

    master_columns = [
        'trait', 'trait_category', 'gene_id', 'gene_name', 'p_acat_o',
        'lead_sv_id', 'lead_sv_type', 'lead_sv_size', 'lead_sv_maf',
        'lead_sv_p', 'lead_sv_layer', 'n_strata', 'n_sv_total',
    ]
    candidate = hit_master[master_columns].merge(
        background.drop(columns=['p_acat_o', 'n_strata', 'n_sv_total']),
        on=['trait', 'gene_id'], how='left', validate='one_to_one',
    )
    missing_candidates = candidate['overall_p'].isna()
    assert not missing_candidates.any(), candidate.loc[missing_candidates, ['trait', 'gene_id']]
    candidate['candidate_row_index'] = np.arange(len(candidate))
    candidate.to_csv(OUTPUT_ROOT / 'step27_candidate_gene_trait_evidence.csv', index=False)

    concordance_rows = []
    for variant_set in VARIANT_SETS:
        prefix = '' if variant_set == 'snv' else 'small_variant_'
        concordance_rows.append({
            'variant_set': variant_set, 'n_gene_trait_hits': len(candidate),
            **{f'n_{outcome}': int(candidate[f'{prefix}{outcome}'].sum()) for outcome in OUTCOMES},
            **{f'rate_{outcome}': float(candidate[f'{prefix}{outcome}'].mean()) for outcome in OUTCOMES},
        })
    concordance = pd.DataFrame(concordance_rows)
    concordance.to_csv(OUTPUT_ROOT / 'step27_observed_concordance.csv', index=False)

    # Analysis block 7
    MATCH_FEATURES = [
        'log_membership_length', 'log_overall_n_snv', 'log_functional_n_snv',
        'log_intronic_n_snv', 'log_n_sv_total', 'n_strata',
    ]
    candidate_gene_ids = set(candidate['gene_id'])
    slot = (candidate[['gene_id', 'gene_name']].drop_duplicates()
            .sort_values('gene_id').reset_index(drop=True))
    slot['slot_index'] = np.arange(len(slot))
    slot_index_by_gene = slot.set_index('gene_id')['slot_index'].to_dict()
    candidate['slot_index'] = candidate['gene_id'].map(slot_index_by_gene).astype(int)

    gene_ids = sorted(background['gene_id'].unique())
    gene_to_index = {gene_id: index for index, gene_id in enumerate(gene_ids)}
    slot_pools = []
    diagnostic_rows = []
    max_pool_size = max(POOL_SIZES)

    for slot_row in slot.itertuples(index=False):
        relevant_traits = sorted(candidate.loc[candidate['gene_id'].eq(slot_row.gene_id), 'trait'].unique())
        relevant = background[background['trait'].isin(relevant_traits)].copy()
        target = relevant[relevant['gene_id'].eq(slot_row.gene_id)][MATCH_FEATURES].mean().to_numpy(float)
        controls = relevant[~relevant['gene_id'].isin(candidate_gene_ids)].copy()
        grouped = controls.groupby('gene_id', sort=False)
        control_features = grouped[MATCH_FEATURES].mean()
        n_traits_tested = grouped['trait'].nunique()
        control_features = control_features.loc[n_traits_tested.eq(len(relevant_traits))].copy()
        assert len(control_features) >= max_pool_size, (slot_row.gene_id, len(control_features))
        matrix = control_features.to_numpy(float)
        center = np.median(matrix, axis=0)
        q25, q75 = np.quantile(matrix, [0.25, 0.75], axis=0)
        scale = q75 - q25
        fallback = np.std(matrix, axis=0)
        scale = np.where(scale > 0, scale, np.where(fallback > 0, fallback, 1.0))
        distance = np.sqrt(np.mean(np.square((matrix - target) / scale), axis=1))
        nearest = np.argsort(distance)[:max_pool_size]
        pool_gene_ids = control_features.index.to_numpy()[nearest]
        pool_indices = np.array([gene_to_index[x] for x in pool_gene_ids], dtype=np.int32)
        slot_pools.append(pool_indices)
        diagnostic_rows.append({
            'slot_index': slot_row.slot_index, 'gene_id': slot_row.gene_id,
            'gene_name': slot_row.gene_name, 'n_relevant_traits': len(relevant_traits),
            'n_eligible_controls': len(control_features),
            'min_distance': float(distance[nearest[0]]),
            'median_distance_pool_200': float(np.median(distance[nearest[:PRIMARY_POOL_SIZE]])),
            'max_distance_pool_200': float(distance[nearest[PRIMARY_POOL_SIZE - 1]]),
        })

    match_diagnostics = pd.DataFrame(diagnostic_rows)
    assert len(slot_pools) == len(slot)
    match_diagnostics.to_csv(OUTPUT_ROOT / 'step27_match_diagnostics.csv', index=False)

    # Analysis block 9
    traits = sorted(candidate['trait'].unique())
    trait_to_index = {trait: index for index, trait in enumerate(traits)}
    n_traits, n_genes = len(traits), len(gene_ids)
    sig_arrays = {outcome: np.zeros((n_traits, n_genes), dtype=bool) for outcome in ALL_OUTCOME_COLUMNS}
    score_arrays = {
        'snv': np.full((n_traits, n_genes), np.nan, dtype=np.float32),
        'small_variant': np.full((n_traits, n_genes), np.nan, dtype=np.float32),
    }
    for row in background[background['trait'].isin(traits)].itertuples(index=False):
        ti, gi = trait_to_index[row.trait], gene_to_index[row.gene_id]
        for outcome in ALL_OUTCOME_COLUMNS:
            sig_arrays[outcome][ti, gi] = bool(getattr(row, outcome))
        score_arrays['snv'][ti, gi] = -np.log10(max(float(row.overall_percentile), 1e-12))
        score_arrays['small_variant'][ti, gi] = -np.log10(
            max(float(row.small_variant_overall_percentile), 1e-12)
        )

    candidate_trait_indices = candidate['trait'].map(trait_to_index).to_numpy(int)
    candidate_slot_indices = candidate['slot_index'].to_numpy(int)

    pair = (candidate.groupby(['lead_sv_id', 'trait'], as_index=False)
            .agg(trait_category=('trait_category', 'first'),
                 n_candidate_genes=('gene_id', 'size'),
                 **{outcome: (outcome, 'max') for outcome in ALL_OUTCOME_COLUMNS}))
    pair['pair_index'] = np.arange(len(pair))
    pair_key_to_index = pair.set_index(['lead_sv_id', 'trait'])['pair_index'].to_dict()
    candidate['pair_index'] = [pair_key_to_index[(x, y)] for x, y in zip(candidate['lead_sv_id'], candidate['trait'])]
    candidate_rows_by_pair = [candidate.index[candidate['pair_index'].eq(i)].to_numpy(int) for i in range(len(pair))]
    candidate_rows_by_slot = [candidate.index[candidate['slot_index'].eq(i)].to_numpy(int) for i in range(len(slot))]

    def summarize_null(pool_size, variant_set, unit, outcome, observed, denominator, null_counts):
        expected = float(np.mean(null_counts))
        return {
            'pool_size': pool_size, 'variant_set': variant_set, 'unit': unit, 'outcome': outcome,
            'observed_hits': int(observed), 'denominator': int(denominator),
            'observed_rate': observed / denominator,
            'null_mean_hits': expected,
            'null_2.5pct_hits': float(np.quantile(null_counts, 0.025)),
            'null_97.5pct_hits': float(np.quantile(null_counts, 0.975)),
            'fold_enrichment': observed / expected if expected > 0 else np.inf,
            'empirical_upper_p': (1 + int((null_counts >= observed).sum())) / (N_PERMUTATIONS + 1),
        }

    summary_rows, null_frames = [], []
    primary_gene_matrices = {}
    primary_sampled_slots = None

    for pool_size in POOL_SIZES:
        rng = np.random.default_rng(RANDOM_SEED + pool_size)
        sampled_slots = [
            pool[:pool_size][rng.integers(0, pool_size, size=N_PERMUTATIONS)].astype(np.int32)
            for pool in slot_pools
        ]
        null_frame = pd.DataFrame({'pool_size': pool_size, 'permutation': np.arange(N_PERMUTATIONS)})
        for variant_set in VARIANT_SETS:
            prefix = '' if variant_set == 'snv' else 'small_variant_'
            for outcome in OUTCOMES:
                outcome_column = f'{prefix}{outcome}'
                gene_matrix = np.vstack([
                    sig_arrays[outcome_column][trait_index, sampled_slots[slot_index]]
                    for trait_index, slot_index in zip(candidate_trait_indices, candidate_slot_indices)
                ])
                pair_matrix = np.vstack([gene_matrix[rows].any(axis=0) for rows in candidate_rows_by_pair])
                cluster_matrix = np.vstack([gene_matrix[rows].any(axis=0) for rows in candidate_rows_by_slot])
                values = {
                    'gene_trait': (int(candidate[outcome_column].sum()), len(candidate), gene_matrix.sum(axis=0)),
                    'candidate_gene': (
                        int(candidate.groupby('gene_id')[outcome_column].max().sum()),
                        len(slot), cluster_matrix.sum(axis=0),
                    ),
                    'lead_sv_trait_pair': (int(pair[outcome_column].sum()), len(pair), pair_matrix.sum(axis=0)),
                }
                for unit, (observed, denominator, null_counts) in values.items():
                    summary_rows.append(summarize_null(
                        pool_size, variant_set, unit, outcome, observed, denominator, null_counts
                    ))
                    if outcome in ['overall_sig', 'any_support_sig']:
                        null_frame[f'{variant_set}_{unit}_{outcome}_hits'] = null_counts
                if pool_size == PRIMARY_POOL_SIZE and outcome in ['overall_sig', 'any_support_sig']:
                    primary_gene_matrices[(variant_set, outcome)] = gene_matrix.copy()
        if pool_size == PRIMARY_POOL_SIZE:
            primary_sampled_slots = [x.copy() for x in sampled_slots]
        null_frames.append(null_frame)

    enrichment = pd.DataFrame(summary_rows)
    null_distribution = pd.concat(null_frames, ignore_index=True)
    enrichment.to_csv(OUTPUT_ROOT / 'step27_cluster_enrichment_summary.csv', index=False)
    null_distribution.to_csv(OUTPUT_ROOT / 'step27_cluster_null_distribution.csv.gz', index=False, compression='gzip')

    # Analysis block 11
    continuous_rows = []
    for variant_set in VARIANT_SETS:
        percentile_column = 'overall_percentile' if variant_set == 'snv' else 'small_variant_overall_percentile'
        observed_scores = -np.log10(np.clip(candidate[percentile_column].to_numpy(float), 1e-12, 1))
        null_score_matrix = np.vstack([
            score_arrays[variant_set][trait_index, primary_sampled_slots[slot_index]]
            for trait_index, slot_index in zip(candidate_trait_indices, candidate_slot_indices)
        ])
        assert np.isfinite(null_score_matrix).all()
        null_mean_scores = null_score_matrix.mean(axis=0)
        observed_mean_score = float(observed_scores.mean())
        continuous_rows.append({
            'variant_set': variant_set,
            'endpoint': 'mean_negative_log10_overall_empirical_percentile',
            'observed': observed_mean_score, 'null_mean': float(null_mean_scores.mean()),
            'null_2.5pct': float(np.quantile(null_mean_scores, 0.025)),
            'null_97.5pct': float(np.quantile(null_mean_scores, 0.975)),
            'empirical_upper_p': (1 + int((null_mean_scores >= observed_mean_score).sum())) / (N_PERMUTATIONS + 1),
        })
    continuous = pd.DataFrame(continuous_rows)
    continuous.to_csv(OUTPUT_ROOT / 'step27_continuous_rank_summary.csv', index=False)

    def subset_summary(variant_set, exclusion_type, excluded_label, mask, outcome):
        prefix = '' if variant_set == 'snv' else 'small_variant_'
        outcome_column = f'{prefix}{outcome}'
        observed = int(candidate.loc[mask, outcome_column].sum())
        null_counts = primary_gene_matrices[(variant_set, outcome)][mask].sum(axis=0)
        expected = float(null_counts.mean())
        return {
            'variant_set': variant_set, 'exclusion_type': exclusion_type, 'excluded_label': excluded_label,
            'outcome': outcome, 'n_gene_trait_rows_retained': int(mask.sum()),
            'observed_hits': observed, 'null_mean_hits': expected,
            'fold_enrichment': observed / expected if expected > 0 else np.inf,
            'empirical_upper_p': (1 + int((null_counts >= observed).sum())) / (N_PERMUTATIONS + 1),
        }

    robustness_rows = []
    for variant_set in VARIANT_SETS:
        for category in sorted(candidate['trait_category'].unique()):
            mask = ~candidate['trait_category'].eq(category)
            for outcome in ['overall_sig', 'any_support_sig']:
                robustness_rows.append(subset_summary(variant_set, 'trait_category', category, mask, outcome))
        for trait in sorted(candidate['trait'].unique()):
            mask = ~candidate['trait'].eq(trait)
            for outcome in ['overall_sig', 'any_support_sig']:
                robustness_rows.append(subset_summary(variant_set, 'trait', trait, mask, outcome))
        for gene_id in sorted(candidate['gene_id'].unique()):
            mask = ~candidate['gene_id'].eq(gene_id)
            for outcome in ['overall_sig', 'any_support_sig']:
                robustness_rows.append(subset_summary(variant_set, 'candidate_gene', gene_id, mask, outcome))
    distance_cutoff = match_diagnostics['median_distance_pool_200'].quantile(0.90)
    well_matched_slots = set(match_diagnostics.loc[match_diagnostics['median_distance_pool_200'] <= distance_cutoff, 'slot_index'])
    mask = candidate['slot_index'].isin(well_matched_slots)
    for variant_set in VARIANT_SETS:
        for outcome in ['overall_sig', 'any_support_sig']:
            robustness_rows.append(subset_summary(
                variant_set, 'worst_10pct_match', f'median_distance>{distance_cutoff:.6g}', mask, outcome
            ))
    robustness = pd.DataFrame(robustness_rows)
    robustness.to_csv(OUTPUT_ROOT / 'step27_leave_out_robustness.csv', index=False)

    # Analysis block 13
    def bh_adjust(pvalues):
        pvalues = np.asarray(pvalues, float)
        order = np.argsort(pvalues)
        ranked = pvalues[order]
        adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
        result = np.empty_like(adjusted)
        result[order] = np.minimum(adjusted, 1.0)
        return result

    category_rows = []
    for variant_set in VARIANT_SETS:
      prefix = '' if variant_set == 'snv' else 'small_variant_'
      for outcome in ['overall_sig', 'any_support_sig']:
        outcome_column = f'{prefix}{outcome}'
        matrix = primary_gene_matrices[(variant_set, outcome)]
        for category_name in sorted(candidate['trait_category'].unique()):
            mask = candidate['trait_category'].eq(category_name).to_numpy()
            observed = int(candidate.loc[mask, outcome_column].sum())
            null_counts = matrix[mask].sum(axis=0)
            expected = float(null_counts.mean())
            category_rows.append({
                'variant_set': variant_set, 'trait_category': category_name, 'outcome': outcome,
                'n_gene_trait_rows': int(mask.sum()), 'observed_hits': observed,
                'observed_rate': observed / mask.sum(), 'null_mean_hits': expected,
                'fold_enrichment': observed / expected if expected > 0 else np.inf,
                'empirical_upper_p': (1 + int((null_counts >= observed).sum())) / (N_PERMUTATIONS + 1),
            })
    category = pd.DataFrame(category_rows)
    category['fdr_bh_within_outcome'] = category.groupby(
        ['variant_set', 'outcome']
    )[ 'empirical_upper_p'].transform(bh_adjust)
    category.to_csv(OUTPUT_ROOT / 'step27_category_enrichment.csv', index=False)

    heterogeneity_rows = []
    for variant_set in VARIANT_SETS:
      prefix = '' if variant_set == 'snv' else 'small_variant_'
      for outcome in ['overall_sig', 'any_support_sig']:
        outcome_column = f'{prefix}{outcome}'
        matrix = primary_gene_matrices[(variant_set, outcome)]
        categories = sorted(candidate['trait_category'].unique())
        masks = [candidate['trait_category'].eq(x).to_numpy() for x in categories]
        sizes = np.array([mask.sum() for mask in masks], float)
        observed_rates = np.array([candidate.loc[mask, outcome_column].mean() for mask in masks])
        observed_global = candidate[outcome_column].mean()
        observed_stat = float(np.sum(sizes * np.square(observed_rates - observed_global)))
        null_category_rates = np.vstack([matrix[mask].mean(axis=0) for mask in masks])
        null_global_rates = matrix.mean(axis=0)
        null_stats = np.sum(sizes[:, None] * np.square(null_category_rates - null_global_rates), axis=0)
        heterogeneity_rows.append({
            'variant_set': variant_set, 'outcome': outcome, 'observed_weighted_rate_variance': observed_stat,
            'null_mean': float(null_stats.mean()),
            'empirical_upper_p': (1 + int((null_stats >= observed_stat).sum())) / (N_PERMUTATIONS + 1),
        })
    heterogeneity = pd.DataFrame(heterogeneity_rows)
    heterogeneity.to_csv(OUTPUT_ROOT / 'step27_category_heterogeneity.csv', index=False)

    # Analysis block 15
    primary = enrichment[(enrichment['pool_size'] == PRIMARY_POOL_SIZE)
                         & (enrichment['unit'] == 'gene_trait')
                         & enrichment['outcome'].isin(['overall_sig', 'any_support_sig'])].copy()

    lead_significance_rows = []
    lead_groups = {
        'all': np.ones(len(candidate), dtype=bool),
        'lead_sv_p_lt_5e-8': candidate['lead_sv_p'].lt(5e-8).to_numpy(),
        'lead_sv_p_ge_5e-8': candidate['lead_sv_p'].ge(5e-8).to_numpy(),
    }
    for variant_set in VARIANT_SETS:
        for group_name, mask in lead_groups.items():
            row = subset_summary(variant_set, 'lead_significance', group_name, mask, 'overall_sig')
            row['n_lead_svs'] = candidate.loc[mask, 'lead_sv_id'].nunique()
            lead_significance_rows.append(row)
    lead_significance = pd.DataFrame(lead_significance_rows)
    lead_significance.to_csv(OUTPUT_ROOT / 'step27_lead_sv_significance_sensitivity.csv', index=False)

    size_thresholds = [('all', 0), ('ge_100bp', 100), ('ge_500bp', 500), ('ge_1kb', 1_000),
                       ('ge_10kb', 10_000), ('ge_50kb', 50_000)]
    lead_size_rows = []
    for variant_set in VARIANT_SETS:
        for label, threshold in size_thresholds:
            mask = candidate['lead_sv_size'].abs().ge(threshold).to_numpy()
            row = subset_summary(variant_set, 'lead_sv_size', label, mask, 'overall_sig')
            row['minimum_lead_sv_size_bp'] = threshold
            row['n_lead_svs'] = candidate.loc[mask, 'lead_sv_id'].nunique()
            lead_size_rows.append(row)
    lead_size = pd.DataFrame(lead_size_rows)
    lead_size.to_csv(OUTPUT_ROOT / 'step27_lead_sv_size_sensitivity.csv', index=False)
    qc = pd.DataFrame([
        ('audit rows', len(audit), 77),
        ('exact completed traits', len(complete), 63),
        ('no-exact traits', len(no_exact), 14),
        ('candidate rows represented in joint background', candidate['overall_p'].notna().sum(), len(candidate)),
        ('matching slots with full 500-gene pools', len(slot_pools), len(slot)),
        ('pool-size analyses', enrichment['pool_size'].nunique(), len(POOL_SIZES)),
        ('variant-set analyses', enrichment['variant_set'].nunique(), len(VARIANT_SETS)),
        ('permutations per pool', N_PERMUTATIONS, N_PERMUTATIONS),
    ], columns=['check', 'observed', 'expected'])
    qc['passed'] = qc['observed'] == qc['expected']
    qc.to_csv(OUTPUT_ROOT / 'step27_qc.csv', index=False)
    assert qc['passed'].all(), 'Step 27 QC failed'

    decision = primary[[
        'variant_set', 'outcome', 'observed_hits', 'denominator', 'observed_rate', 'null_mean_hits',
        'fold_enrichment', 'empirical_upper_p',
    ]].copy()
    decision = decision.merge(
        continuous[['variant_set', 'empirical_upper_p']].rename(
            columns={'empirical_upper_p': 'continuous_rank_empirical_p'}
        ),
        on='variant_set', how='left', validate='many_to_one',
    )
    decision['interpretation_boundary'] = 'Gene-level convergence only; not colocalization or shared causality'
    decision.to_csv(OUTPUT_ROOT / 'step27_decision_summary.csv', index=False)


if __name__ == "__main__":
    main()
