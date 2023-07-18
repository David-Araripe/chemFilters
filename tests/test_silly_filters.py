import unittest
from pathlib import Path

import pandas as pd

from chemFilters.silly_filters import SillyMolSpotterFilter


class TestSillyMolSpotterFilter(unittest.TestCase):
    def setUp(self) -> None:
        self.testroot = Path(__file__).parent
        self.filterFunc = SillyMolSpotterFilter(
            from_smi=True, standardize=True, std_method="chembl", n_jobs=1
        )
        self.test_smiles = (
            (self.testroot / "resources/testSmiles.smi").read_text().splitlines()
        )
        self.expected_scoring_df = pd.read_csv(
            self.testroot / "resources/testSmiles_sillyfiltered.csv"
        )

    def test_get_spotters(self):
        spotters = self.filterFunc._get_spotters()
        self.assertIsInstance(spotters, dict)

    def test_score_compound(self):
        try:
            score = self.filterFunc.score_compound(self.test_smiles[0])
            self.assertIsNotNone(score)
        except Exception as e:
            self.fail(f"test_score_compound raised an exception: {e}")

    def test_get_scoring_df(self):
        result_df = self.filterFunc.get_scoring_df(self.test_smiles)
        pd.testing.assert_frame_equal(result_df, self.expected_scoring_df)

    def tearDown(self) -> None:
        return super().tearDown()


if __name__ == "__main__":
    unittest.main()
