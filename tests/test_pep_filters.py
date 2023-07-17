import unittest
from pathlib import Path

import pandas as pd

from chemFilters import PeptideFilters


class TestPeptideFilters(unittest.TestCase):
    def setUp(self) -> None:
        self.testroot = Path(__file__).parent
        self.filterFunc = PeptideFilters(filterType="all", njobs=1, from_smi=True)
        self.test_smiles = (
            (self.testroot / "resources/testSmiles.smi").read_text().splitlines()
        )
        self.expected_filtered = pd.read_csv(
            self.testroot / "resources/testSmiles_pepfiltered.csv"
        )

    def test_available_filters(self):
        filters = self.filterFunc.available_filters
        self.assertIsInstance(filters, dict)  # Expected type may vary

    def test_filter_mols(self):
        filtered_mols = self.filterFunc.filter_mols(self.test_smiles)
        self.assertIsNotNone(
            filtered_mols
        )  # Depending on the filter_mols functionality, you may need to adjust the assertion

    def test_get_flagging_df(self):
        result_df = self.filterFunc.get_flagging_df(self.test_smiles)
        pd.testing.assert_frame_equal(result_df, self.expected_filtered)

    def tearDown(self) -> None:
        return super().tearDown()


if __name__ == "__main__":
    unittest.main()
