"""Small data-free checks for the published analysis primitives."""

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from acat import acat_pvalue, run_acat_gene  # noqa: E402


def load_reference_module():
    path = ROOT / "scripts" / "01_reference.py"
    spec = importlib.util.spec_from_file_location("prepare_reference", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_script(name):
    path = ROOT / "scripts" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CoreMethodTests(unittest.TestCase):
    def test_output_layout(self):
        module = load_script("02_sv_acat.py")
        data_root = ROOT / "synthetic_data_root"
        module.set_project_root(data_root)
        self.assertEqual(module.TRAITS_DIR, data_root / "results" / "traits")
        self.assertEqual(module.GENE_TABLE, data_root / "results" / "gene_table.parquet")

    def test_residual_sv_summary_retains_null_interval(self):
        module = load_script("03_lead_architecture.py")
        module.N_PERMUTATIONS = 1_000
        locality = pd.DataFrame([{
            "candidate_distances_bp": "0; 100000",
            "candidate_same_type": "1; 0",
            "observed_distance_bp": 0,
            "observed_overlap": 1,
            "observed_within_50kb": 1,
            "observed_within_250kb": 1,
            "observed_same_type": 1,
            "observed_same_type_within_50kb": 1,
        }])
        summary, _ = module.locality_permutation(locality, np.random.default_rng(1))
        self.assertTrue((summary["null_95pct_low"] <= summary["null_mean"]).all())
        self.assertTrue((summary["null_mean"] <= summary["null_95pct_high"]).all())

    def test_sv_type_stratified_acat(self):
        variants = pd.DataFrame({
            "gene_id": ["G1", "G1", "G1"],
            "gene_name": ["GENE", "GENE", "GENE"],
            "sv_type": ["DEL", "DEL", "DUP"],
            "pval": [0.001, 0.1, 0.02],
            "w_final": [2.0, 1.0, 3.0],
        })
        strata, genes = run_acat_gene(variants)
        del_p = acat_pvalue([0.001, 0.1], [2.0, 1.0])
        expected = acat_pvalue([del_p, 0.02])
        self.assertEqual(strata["sv_type"].tolist(), ["DEL", "DUP"])
        self.assertEqual(int(genes.iloc[0]["n_sv_total"]), 3)
        self.assertAlmostEqual(float(genes.iloc[0]["p_acat_o"]), expected)

    def test_gene_body_coordinates_on_both_strands(self):
        module = load_reference_module()
        attributes = lambda gene: (
            f'gene_id "{gene}"; gene_name "{gene}"; gene_type "protein_coding";'
        )
        records = pd.DataFrame([
            ["chr1", "GENCODE", "gene", 101, 200, ".", "+", ".", attributes("PLUS")],
            ["chr2", "GENCODE", "gene", 301, 450, ".", "-", ".", attributes("MINUS")],
        ], columns=module.GTF_COLUMNS)
        with patch.object(module.pd, "read_csv", return_value=[records]):
            table = module.build_gene_table(Path("synthetic.gtf"))
        self.assertEqual(table["gene_start"].tolist(), [100, 300])
        self.assertEqual(table["gene_end"].tolist(), [200, 450])
        self.assertEqual(table["tss"].tolist(), [100, 450])
        self.assertTrue(np.all(table["promoter_start"].to_numpy() >= 0))

    def test_insertion_anchor_and_effect_columns(self):
        module = load_script("02_sv_acat.py")
        raw = pd.DataFrame({
            "Name": ["ins_left", "ins_right", "del_span"],
            "Chrom": ["1", "1", "1"],
            "Pos": [101, 102, 100],
            "effectAllele": [
                "<INS:SVSIZE=50:TEST>", "<INS:SVSIZE=50:TEST>",
                "<DEL:SVSIZE=3:TEST>",
            ],
            "effectAlleleFreq": [0.01, 0.02, 0.03],
            "Beta": [0.1, 0.2, 0.3],
            "SE": [0.01, 0.02, 0.03],
            "N": [100, 200, 300],
            "pval": [0.1, 0.2, 0.3],
        })
        genes = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "gene_id": ["left", "right"],
            "gene_name": ["LEFT", "RIGHT"],
            "tss": [100, 101],
            "promoter_start": [100, 101],
            "promoter_end": [101, 102],
            "gene_start": [100, 101],
            "gene_end": [101, 110],
        })
        features = pd.DataFrame(columns=[
            "chrom", "gene_id", "feat_type", "feat_start", "feat_end",
        ])
        saved = []
        with tempfile.TemporaryDirectory() as tmp:
            module.set_project_root(tmp)
            with patch.object(module.pd, "read_csv", return_value=raw):
                with patch.object(pd.DataFrame, "to_parquet", lambda frame, *args, **kwargs: saved.append(frame.copy())):
                    self.assertTrue(module.run_step2("synthetic", genes, features))
        annotated = saved[0].set_index(["sv_id", "gene_id"])
        self.assertEqual(set(annotated.index), {
            ("ins_left", "left"), ("ins_right", "right"),
            ("del_span", "left"), ("del_span", "right"),
        })
        self.assertEqual(int(annotated.loc[("ins_left", "left"), "sv_end"]), 101)
        self.assertEqual(int(annotated.loc[("del_span", "right"), "sv_end"]), 102)
        self.assertEqual(annotated.loc[("ins_left", "left"), ["Beta", "SE", "N"]].tolist(), [0.1, 0.01, 100])

    def test_merged_stage_order(self):
        cases = [
            ("03_lead_architecture.py", [
                "_stage_primary", "_stage_decomposition", "_stage_criticality",
            ]),
            ("04_recurrence.py", [
                "_stage_qc", "_stage_enrichment", "_stage_locus_audit",
            ]),
        ]
        for filename, names in cases:
            with self.subTest(script=filename):
                module = load_script(filename)
                called = []
                with patch("sys.argv", [filename, "--project-root", "/private/example"]):
                    with patch.multiple(module, **{
                        name: lambda args, stage=name: called.append(stage)
                        for name in names
                    }):
                        module.main()
                self.assertEqual(called, names)


if __name__ == "__main__":
    unittest.main()
