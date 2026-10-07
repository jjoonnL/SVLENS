"""Small data-free checks for the published analysis primitives."""

import importlib.util
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from svarch.acat import acat_pvalue, run_acat_gene  # noqa: E402


def load_reference_module():
    path = ROOT / "scripts" / "01_prepare_reference.py"
    spec = importlib.util.spec_from_file_location("prepare_reference", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CoreMethodTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
