import unittest
from pathlib import Path

import pandas as pd
from rdkit import Chem

from chemFilters import SillyMolSpotterFilter


class TestSillyMolSpotterFilter(unittest.TestCase):
    def setUp(self) -> None:
        self.testroot = Path(__file__).parent
        self.filterFunc = SillyMolSpotterFilter(
            from_smi=True, standardize=True, std_method="chembl", n_jobs=1
        )
        self.test_smiles = (
            (self.testroot / "resources/testSmiles.smi").read_text().splitlines()
        )
        self.test_mols = [Chem.MolFromSmiles(smi) for smi in self.test_smiles]
        self.expected_scoring_df = pd.read_csv(
            self.testroot / "resources/testSmiles_sillyfiltered.csv"
        )

    def test_get_spotters(self):
        spotters = self.filterFunc._get_spotters()
        self.assertIsInstance(spotters, dict)

    def test_score_smi(self):
        try:
            score = self.filterFunc.score_smi(self.test_smiles[0])
            self.assertIsNotNone(score)
        except Exception as e:
            self.fail(f"test_score_smi raised an exception: {e}")

    def test_get_scoring_df(self):
        result_df = self.filterFunc.get_scoring_df(self.test_smiles)
        pd.testing.assert_frame_equal(result_df, self.expected_scoring_df)

    def test_get_scoring_df_variations(self):
        filterTT = SillyMolSpotterFilter(  # Filter true true
            from_smi=True, standardize=True, std_method="chembl", n_jobs=1
        ).get_scoring_df(self.test_smiles)
        filterTF = SillyMolSpotterFilter(
            from_smi=True, standardize=False, std_method="chembl", n_jobs=1
        ).get_scoring_df(self.test_smiles)
        filterFT = SillyMolSpotterFilter(
            from_smi=False, standardize=True, std_method="chembl", n_jobs=1
        ).get_scoring_df(self.test_mols)
        filterFF = SillyMolSpotterFilter(
            from_smi=False, standardize=False, std_method="chembl", n_jobs=1
        ).get_scoring_df(self.test_mols)
        pd.testing.assert_frame_equal(filterTT, filterTF)
        pd.testing.assert_frame_equal(filterTT, filterFT)
        pd.testing.assert_frame_equal(filterTT, filterFF)

    def tearDown(self) -> None:
        return super().tearDown()


if __name__ == "__main__":
    unittest.main()
