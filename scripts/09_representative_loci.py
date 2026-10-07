"""Select representative complete, partial and absent rare-SNV-support loci."""

from __future__ import annotations

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()

    # Inputs and output location
    import numpy as np
    import pandas as pd

    PROJECT_ROOT = args.project_root.expanduser().resolve()
    RESULT_ROOT = PROJECT_ROOT / 'results' / 'sv_pleiotropy'
    OUTPUT_DIR = RESULT_ROOT / 'candidate_locus_dossiers_gene_body'
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    PATHS = {
        'master': RESULT_ROOT / 'sv_pleiotropy_master.non_ratio_primary.csv',
        'recurrent': RESULT_ROOT / 'candidate_recurrent_lead_svs.non_ratio_primary.csv',
        'origin_audit': RESULT_ROOT / 'somatic_immune_audit/lead_sv_origin_audit.csv',
        'target_all': RESULT_ROOT / 'target_ambiguity/recurrent_lead_sv_target_ambiguity.csv',
        'driver_all': RESULT_ROOT / 'driver_decomposition/lead_sv_driver_summary.non_ratio_primary.csv',
        'redundancy_all': RESULT_ROOT / 'lead_sv_redundancy_sensitivity.csv',
        'gene_table': PROJECT_ROOT / 'results/gene_table.parquet',
        'step27_evidence': RESULT_ROOT / 'full_sv_rare_snv_convergence_gene_body/step27_candidate_gene_trait_evidence.csv',
        'step27_decision': RESULT_ROOT / 'full_sv_rare_snv_convergence_gene_body/step27_decision_summary.csv',
    }
    missing_required = [str(path) for path in PATHS.values() if not path.exists()]
    if missing_required:
        raise FileNotFoundError('Missing required inputs:\n' + '\n'.join(missing_required))

    # Analysis block 3
    master = pd.read_csv(PATHS['master'])
    recurrent = pd.read_csv(PATHS['recurrent'])
    origin_audit = pd.read_csv(PATHS['origin_audit'])
    target_all = pd.read_csv(PATHS['target_all'])
    driver_all = pd.read_csv(PATHS['driver_all'])
    redundancy_all = pd.read_csv(PATHS['redundancy_all'])
    gene_table = pd.read_parquet(PATHS['gene_table'])
    step27_evidence = pd.read_csv(PATHS['step27_evidence'])
    step27_decision = pd.read_csv(PATHS['step27_decision'])

    required_step27_columns = {
        'gene_id', 'gene_name', 'trait', 'overall_sig', 'small_variant_overall_sig',
        'overall_p', 'small_variant_overall_p', 'p_acat_o', 'lead_sv_id',
    }
    if not required_step27_columns.issubset(step27_evidence.columns):
        raise ValueError(f'Missing corrected Step 27 columns: {sorted(required_step27_columns - set(step27_evidence.columns))}')
    assert not step27_evidence.duplicated(['gene_id', 'trait']).any()

    if recurrent['lead_sv_id'].nunique() != 64:
        raise AssertionError(f"Expected 64 non-ratio recurrent lead SVs, found {recurrent['lead_sv_id'].nunique()}")

    spectrum = step27_evidence.groupby(['lead_sv_id', 'gene_id', 'gene_name'], as_index=False).agg(
        n_evaluable_traits=('trait', 'nunique'),
        n_snv_supported=('overall_sig', 'sum'),
        n_small_variant_supported=('small_variant_overall_sig', 'sum'),
        min_sv_gene_p=('p_acat_o', 'min'),
    )
    spectrum['snv_support_fraction'] = spectrum['n_snv_supported'] / spectrum['n_evaluable_traits']
    spectrum['small_variant_support_fraction'] = (
        spectrum['n_small_variant_supported'] / spectrum['n_evaluable_traits']
    )
    spectrum['snv_support_class'] = np.select(
        [spectrum['n_snv_supported'].eq(spectrum['n_evaluable_traits']), spectrum['n_snv_supported'].eq(0)],
        ['complete', 'absent'], default='partial',
    )

    family = redundancy_all[redundancy_all['scenario'].eq('family_no_ratio')][[
        'lead_sv_id', 'n_units', 'broad_fdr_lt_0_05', 'fine_fdr_lt_0_05',
        'any_module_fdr_lt_0_05', 'broad_purity_fdr', 'fine_purity_fdr',
    ]].rename(columns={
        'n_units': 'n_families_non_ratio',
        'broad_fdr_lt_0_05': 'family_broad_fdr_lt_0_05',
        'fine_fdr_lt_0_05': 'family_fine_fdr_lt_0_05',
        'any_module_fdr_lt_0_05': 'family_any_module_fdr_lt_0_05',
        'broad_purity_fdr': 'family_broad_purity_fdr',
        'fine_purity_fdr': 'family_fine_purity_fdr',
    })

    selection = recurrent.merge(
        spectrum, on='lead_sv_id', how='left', validate='one_to_many'
    ).merge(
        origin_audit[['lead_sv_id', 'blood_derived_candidate', 'origin_class']],
        on='lead_sv_id', how='left', validate='many_to_one'
    ).merge(
        target_all[['lead_sv_id', 'target_ambiguity_class']],
        on='lead_sv_id', how='left', validate='many_to_one'
    ).merge(
        driver_all[['lead_sv_id', 'gene_name', 'frac_lead_anchored_or_single', 'any_without_lead_bonf']],
        left_on=['lead_sv_id', 'gene_name'], right_on=['lead_sv_id', 'gene_name'],
        how='left', validate='many_to_one'
    ).merge(family, on='lead_sv_id', how='left', validate='many_to_one')

    selection['origin_ok'] = selection['blood_derived_candidate'].eq(False)
    selection['qc_ok'] = selection['case_tier'].ne('QC_caution')
    selection['single_gene_target'] = selection['target_ambiguity_class'].eq('single_gene_target')
    selection['lead_dependency_ok'] = (
        selection['frac_lead_anchored_or_single'].eq(1.0)
        & selection['any_without_lead_bonf'].eq(False)
    )
    selection['step27_evaluable_ok'] = selection['n_evaluable_traits'].ge(2)
    selection['eligible_for_dossier'] = selection[[
        'origin_ok', 'qc_ok', 'single_gene_target', 'lead_dependency_ok', 'step27_evaluable_ok',
    ]].all(axis=1)

    def exclusion_reason(row):
        reasons = []
        if not row['origin_ok']:
            reasons.append('origin_audit_flag')
        if not row['qc_ok']:
            reasons.append('qc_caution')
        if not row['single_gene_target']:
            reasons.append('multi_gene_target')
        if not row['lead_dependency_ok']:
            reasons.append('not_fully_lead_anchored')
        if not row['step27_evaluable_ok']:
            reasons.append('fewer_than_two_evaluable_traits')
        return ';'.join(reasons) if reasons else 'eligible'

    selection['selection_status'] = selection.apply(exclusion_reason, axis=1)
    selection['selected_for_dossier'] = False

    # Uniform ranking within each support class: family-collapsed module robustness,
    # non-ratio module enrichment, evaluable trait count, then SV association strength.
    rank_columns = [
        'family_any_module_fdr_lt_0_05', 'any_module_fdr_lt_0_05',
        'n_evaluable_traits', 'min_sv_gene_p', 'gene_name',
    ]
    rank_ascending = [False, False, False, True, True]
    selected_rows = []
    for support_class in ['complete', 'partial', 'absent']:
        ranked = selection[
            selection['eligible_for_dossier'] & selection['snv_support_class'].eq(support_class)
        ].sort_values(rank_columns, ascending=rank_ascending)
        if ranked.empty:
            raise AssertionError(f'No eligible {support_class} rare-SNV-support locus')
        selected_rows.append(ranked.iloc[0])

    selected_table = pd.DataFrame(selected_rows)
    selected_keys = set(zip(selected_table['lead_sv_id'], selected_table['gene_name']))
    selection['selected_for_dossier'] = [
        (lead_sv_id, gene_name) in selected_keys
        for lead_sv_id, gene_name in zip(selection['lead_sv_id'], selection['gene_name'])
    ]
    CANDIDATES = selected_table['gene_name'].tolist()
    if len(CANDIDATES) != 3 or len(set(CANDIDATES)) != 3:
        raise AssertionError(f'Expected three distinct spectrum anchors, found {CANDIDATES}')
    role_by_class = {
        'complete': 'Complete rare-SNV support',
        'partial': 'Partial rare-SNV support',
        'absent': 'No significant rare-SNV support',
    }
    candidate_classes = selected_table.set_index('gene_name')['snv_support_class'].to_dict()
    CANDIDATE_ROLES = {gene: role_by_class[candidate_classes[gene]] for gene in CANDIDATES}
    selection.to_csv(OUTPUT_DIR / 'candidate_selection_universe.csv', index=False)
    selection[selection['n_evaluable_traits'].notna()].to_csv(
        OUTPUT_DIR / 'candidate_selection_spectrum.csv', index=False
    )

    candidate_overview = selected_table.copy()
    candidate_overview['dossier_role'] = candidate_overview['gene_name'].map(CANDIDATE_ROLES)
    candidate_overview['approx_mac_median'] = candidate_overview['median_approx_mac']
    categories_by_unit = master.groupby(['lead_sv_id', 'gene_name'])['trait_category'].agg(
        lambda values: '; '.join(sorted(set(values)))
    )
    candidate_overview['categories_non_ratio'] = [
        categories_by_unit.loc[(lead_sv_id, gene_name)]
        for lead_sv_id, gene_name in zip(
            candidate_overview['lead_sv_id'], candidate_overview['gene_name']
        )
    ]
    candidate_overview['genes'] = candidate_overview['gene_name']
    candidate_overview = candidate_overview[[
        'genes', 'dossier_role', 'lead_sv_id', 'lead_sv_type', 'lead_sv_layer',
        'lead_sv_size', 'lead_sv_maf', 'approx_mac_median', 'n_traits_non_ratio',
        'n_families_non_ratio', 'categories_non_ratio', 'support_classes_non_ratio',
    ]].reset_index(drop=True)
    candidate_overview.to_csv(OUTPUT_DIR / 'candidate_overview.csv', index=False)

    # Analysis block 4
    # Export the 32 eligible units in the Supplementary Table S8 column schema.
    eligible = selection.loc[selection['eligible_for_dossier']].copy()
    if len(eligible) != 32:
        raise AssertionError(f'Expected 32 eligible lead-SV–gene units, found {len(eligible)}')
    eligible['class_order'] = eligible['snv_support_class'].map({'complete': 0, 'partial': 1, 'absent': 2})
    eligible = eligible.sort_values(
        ['class_order', *rank_columns], ascending=[True, *rank_ascending]
    ).reset_index(drop=True)
    eligible['selection_rank_within_support_class'] = eligible.groupby('snv_support_class').cumcount() + 1
    eligible['phenotype_categories'] = [
        categories_by_unit.loc[(lead_sv_id, gene_name)]
        for lead_sv_id, gene_name in zip(eligible['lead_sv_id'], eligible['gene_name'])
    ]
    table_s8 = eligible.rename(columns={
        'snv_support_class': 'rare_snv_support_class',
        'lead_sv_size': 'lead_sv_size_bp',
        'n_families_non_ratio': 'n_trait_families',
        'family_any_module_fdr_lt_0_05': 'family_collapsed_module_support',
        'any_module_fdr_lt_0_05': 'individual_trait_module_support',
        'n_snv_supported': 'n_rare_snv_supported',
        'snv_support_fraction': 'rare_snv_support_fraction',
        'min_sv_gene_p': 'min_sv_gene_level_p',
        'origin_ok': 'eligibility_origin_ok',
        'qc_ok': 'eligibility_qc_ok',
        'single_gene_target': 'eligibility_single_gene_target',
        'lead_dependency_ok': 'eligibility_lead_dependency_ok',
        'step27_evaluable_ok': 'eligibility_evaluable_traits_ge_2',
    })[[
        'rare_snv_support_class', 'selection_rank_within_support_class', 'selected_for_dossier',
        'gene_name', 'gene_id', 'lead_sv_id', 'lead_sv_type', 'lead_sv_size_bp',
        'lead_sv_maf', 'phenotype_categories', 'n_trait_families',
        'family_collapsed_module_support', 'individual_trait_module_support',
        'n_evaluable_traits', 'n_rare_snv_supported', 'rare_snv_support_fraction',
        'min_sv_gene_level_p', 'eligibility_origin_ok', 'eligibility_qc_ok',
        'eligibility_single_gene_target', 'eligibility_lead_dependency_ok',
        'eligibility_evaluable_traits_ge_2',
    ]]
    table_s8.to_csv(OUTPUT_DIR / 'supplementary_table_rare_snv_locus_selection.csv', index=False)


if __name__ == "__main__":
    main()
