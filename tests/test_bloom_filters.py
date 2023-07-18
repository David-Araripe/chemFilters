import unittest
from pathlib import Path

import pandas as pd
from rdkit import Chem

from chemFilters import MolbloomFilters


class TestMolbloomFilters(unittest.TestCase):
    def setUp(self) -> None:
        self.testroot = Path(__file__).parent
        self.filterFunc = MolbloomFilters(
            from_smi=True, standardize=True, std_method="chembl", n_jobs=1
        )
        self.test_smiles = (
            (self.testroot / "resources/testSmiles.smi").read_text().splitlines()
        )
        self.test_mols = [Chem.MolFromSmiles(smi) for smi in self.test_smiles]
        self.expected_filtered_df = pd.read_csv(
            self.testroot / "resources/testSmiles_bloomfiltered.csv"
        )

    def test_get_catalogs(self):
        catalogs = self.filterFunc.get_catalogs()
        self.assertIsInstance(catalogs, dict)

    def test_buy_smi(self):
        bought_mol = self.filterFunc.buy_smi(self.test_smiles[0])
        self.assertIsNotNone(bought_mol)
        can_buy = self.filterFunc.buy_smi("CCCO", catalog="zinc-instock")
        self.assertTrue(can_buy)

    def test_get_flagging_df(self):
        result_df = self.filterFunc.get_flagging_df(self.test_smiles)
        pd.testing.assert_frame_equal(result_df, self.expected_filtered_df)

    def test_get_flagging_df_variations(self):
        filterTT = MolbloomFilters(  # Filter true true
            from_smi=True, standardize=True, std_method="chembl", n_jobs=1
        ).get_flagging_df(self.test_smiles)
        filterTF = MolbloomFilters(
            from_smi=True, standardize=False, std_method="chembl", n_jobs=1
        ).get_flagging_df(self.test_smiles)
        filterFT = MolbloomFilters(
            from_smi=False, standardize=True, std_method="chembl", n_jobs=1
        ).get_flagging_df(self.test_mols)
        filterFF = MolbloomFilters(
            from_smi=False, standardize=False, std_method="chembl", n_jobs=1
        ).get_flagging_df(self.test_mols)
        pd.testing.assert_frame_equal(filterTT, filterTF)
        pd.testing.assert_frame_equal(filterTT, filterFT)
        pd.testing.assert_frame_equal(filterTT, filterFF)

    def test_get_flagging_df_no_standardizer(self):
        bloom_no_std = MolbloomFilters(from_smi=True, standardize=False)
        result_df = bloom_no_std.get_flagging_df(self.test_smiles)
        pd.testing.assert_frame_equal(result_df, self.expected_filtered_df)

    def tearDown(self) -> None:
        return super().tearDown()


if __name__ == "__main__":
    unittest.main()
