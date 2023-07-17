"""Unit tests for the chemFilters.rdkit_filters module."""
import unittest
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalogParams

from chemFilters import RdkitFilters


class TestRdkitFilters(unittest.TestCase):
    def setUp(self) -> None:
        self.testroot = Path(__file__).parent
        self.filterFunc_from_smi = RdkitFilters(filter_type="ALL", from_smi=True)
        self.filterFunc_from_mol = RdkitFilters(filter_type="ALL", from_smi=False)
        self.availFilters = self.filterFunc_from_mol.available_filters
        self.test_smiles = (
            (self.testroot / "resources/testSmiles.smi").read_text().splitlines()
        )
        self.expected_filtered_df = pd.read_csv(
            self.testroot / "resources/testSmiles_rdfiltered.csv"
        )

    def test_get_filter(self):
        for _filter in self.availFilters:
            try:
                newFilter = RdkitFilters(_filter)
                result = self.filterFunc_from_mol._get_filter()
            except ValueError:
                self.fail(f"Filter type {_filter} not available.")
            self.assertIsInstance(newFilter, RdkitFilters)
            self.assertIsInstance(result, FilterCatalogParams.FilterCatalogs)

    def test_filter_mols(self):
        mols = [Chem.MolFromSmiles(smi) for smi in self.test_smiles]
        filter_names, descriptions, substructs = self.filterFunc_from_mol.filter_mols(
            mols
        )
        self.assertIsInstance(filter_names, tuple)
        self.assertIsInstance(descriptions, tuple)
        self.assertIsInstance(substructs, tuple)
        self.assertIsInstance(filter_names[0], list)
        self.assertIsInstance(descriptions[0], list)
        self.assertIsInstance(substructs[0], list)
        self.assertIsInstance(filter_names[0][0], str)
        self.assertIsInstance(descriptions[0][0], str)
        self.assertIsInstance(substructs[0][0], Chem.Mol)

    def test_get_flagging_df_from_smi(self):
        result_df = self.filterFunc_from_smi.get_flagging_df(self.test_smiles)
        pd.testing.assert_frame_equal(result_df, self.expected_filtered_df)

    def test_get_flagging_df_from_mol(self):
        mols = [Chem.MolFromSmiles(smi) for smi in self.test_smiles]
        result_df = self.filterFunc_from_mol.get_flagging_df(mols)
        pd.testing.assert_frame_equal(result_df, self.expected_filtered_df)

    def test_unvalid_smiles(self):
        mock_smiles = ["hihihi", "hohoho"]
        result_df = self.filterFunc_from_smi.get_flagging_df(mock_smiles)
        self.assertTrue(result_df["SMILES"].isnull().all())

    def tearDown(self) -> None:
        return super().tearDown()


if __name__ == "__main__":
    unittest.main()
