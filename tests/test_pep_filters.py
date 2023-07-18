import unittest
from pathlib import Path

import pandas as pd
from rdkit import Chem

from chemFilters import PeptideFilters


class TestPeptideFilters(unittest.TestCase):
    def setUp(self) -> None:
        self.testroot = Path(__file__).parent
        self.filterFunc = PeptideFilters(filter_type="all", n_jobs=1, from_smi=True)
        self.test_smiles = (
            (self.testroot / "resources/testSmiles.smi").read_text().splitlines()
        )
        self.test_mols = [Chem.MolFromSmiles(smi) for smi in self.test_smiles]
        self.expected_filtered = pd.read_csv(
            self.testroot / "resources/testSmiles_pepfiltered.csv"
        )

    def test_available_filters(self):
        filters = self.filterFunc.available_filters
        self.assertIsInstance(filters, dict)
        # test initializing an invalid filter
        with self.assertRaises(ValueError):
            PeptideFilters(filter_type=6)
        with self.assertRaises(ValueError):
            PeptideFilters(filter_type="invalid")

    def test_filter_mols(self):
        filtered_mols = self.filterFunc.filter_mols(self.test_smiles)
        self.assertIsNotNone(filtered_mols)

    def test_get_flagging_df(self):
        result_df = self.filterFunc.get_flagging_df(self.test_smiles)
        pd.testing.assert_frame_equal(result_df, self.expected_filtered)
        filterFunc_from_mols = PeptideFilters(  # Filter true true
            from_smi=False, filter_type="all", n_jobs=1
        )
        from_mol_results_df = filterFunc_from_mols.get_flagging_df(self.test_mols)
        pd.testing.assert_frame_equal(result_df, from_mol_results_df)

    def tearDown(self) -> None:
        return super().tearDown()


if __name__ == "__main__":
    unittest.main()
