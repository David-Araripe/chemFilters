import unittest
from pathlib import Path

import pandas as pd

from chemFilters import MolbloomFilters


class TestMolbloomFilters(unittest.TestCase):
    def setUp(self) -> None:
        self.testroot = Path(__file__).parent
        self.filterFunc = MolbloomFilters(from_smi=True, standardize=True, std_method="chembl", n_jobs=1)
        self.test_smiles = (
            (self.testroot / "resources/testSmiles.smi").read_text().splitlines()
        )
        self.expected_filtered = (
            pd.read_csv(self.testroot / "resources/testSmiles_bloomfiltered.csv")
        )

    def test_get_catalogs(self):
        catalogs = self.filterFunc.get_catalogs()
        self.assertIsInstance(catalogs, dict)

    def test_buy_mol(self):
        bought_mol = self.filterFunc.buy_mol(self.test_smiles[0])
        self.assertIsNotNone(bought_mol)
        can_buy = self.filterFunc.buy_mol("CCCO", catalog="zinc-instock")
        self.assertTrue(can_buy)

    def test_get_flagging_df(self):
        result_df = self.filterFunc.get_flagging_df(self.test_smiles)
        result_df.to_csv(self.testroot / "resources/testSmiles_bloomfiltered.csv", index=False)
        pd.testing.assert_frame_equal(result_df, self.expected_filtered)

    def tearDown(self) -> None:
        return super().tearDown()


if __name__ == '__main__':
    unittest.main()
