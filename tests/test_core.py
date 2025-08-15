"""Unit tests for the chemFilters.rdkit_filters module."""

import unittest
from pathlib import Path

import pandas as pd
from rdkit import Chem

from chemFilters.core import CoreFilter

RDCOLS = [
    "NIH",
    "PAINS_A",
    "PAINS_B",
    "PAINS_C",
    "PAINS_any",
    "ZINC",
    "Brenk",
    "ChEMBL23_Dundee",
    "ChEMBL23_BMS",
    "ChEMBL23_MLSMR",
    "ChEMBL23_Inpharmatica",
    "ChEMBL23_SureChEMBL",
    "ChEMBL23_LINT",
    "ChEMBL23_Glaxo",
    "ChEMBL23_any",
]
SIFT_COLS = [
    "SIFT_NaturalLAminoAcids",
    "SIFT_NaturalLDAminoAcids",
    "SIFT_NaturalAminoAcidDerivatives",
    "SIFT_NonNaturalAminoAcidDerivatives",
    "SIFT_AllAmineAndAcid",
    "SIFT_any",
]

BLOOM_COLS = [
    "BLOOM_zinc20",
    "BLOOM_zinc-instock",
    "BLOOM_zinc-instock-mini",
    "BLOOM_surechembl",
    "BLOOM_any",
]
SILLY_COLS = ["SILLY_chembl", "SILLY_excape", "SILLY_papyrus"]


class TestCoreFilter(unittest.TestCase):
    def setUp(self) -> None:
        self.testroot = Path(__file__).parent
        self.standard_params = {
            "rdkit_filter": True,
            "pep_filter": True,
            "silly_filter": True,
            "bloom_filter": True,
            "rdfilter_subset": "ALL",
            "rdfilter_output": "string",
            "std_mols": True,
            "std_method": "chembl",
            "n_jobs": 1,
        }
        self.coreFilter = CoreFilter(**self.standard_params)
        self.test_smiles = (
            (self.testroot / "resources/testSmiles.smi").read_text().splitlines()
        )
        self.test_mols = [Chem.MolFromSmiles(smi) for smi in self.test_smiles]
        self.expected_filtered_df = pd.read_csv(
            self.testroot / "resources/testSmiles_allFiltered.csv"
        )

    def test_filter_smiles(self):
        filtered_df = self.coreFilter.filter_smiles(self.test_smiles, chunksize=1)
        pd.testing.assert_frame_equal(filtered_df, self.expected_filtered_df)

    def test_rdfilterMols(self):
        filtered_df = self.coreFilter._rdfilterMols(self.test_mols)
        # Won't test for the smiles to be the same - not necessarily standardized
        absent_cols = ["SMILES"] + SILLY_COLS + BLOOM_COLS + SIFT_COLS
        expected = self.expected_filtered_df.drop(columns=absent_cols)
        pd.testing.assert_frame_equal(filtered_df.drop(columns="SMILES"), expected)

    def test_pepfilterMols(self):
        filtered_df = self.coreFilter._pepfilterMols(self.test_mols)
        # Won't test for the smiles to be the same - not necessarily standardized
        absent_cols = ["SMILES"] + RDCOLS + SILLY_COLS + BLOOM_COLS
        expected = self.expected_filtered_df.drop(columns=absent_cols)
        pd.testing.assert_frame_equal(filtered_df.drop(columns="SMILES"), expected)

    def test_bloomfilterMols(self):
        filtered_df = self.coreFilter._bloomfilterMols(self.test_mols)
        # Won't test for the smiles to be the same - not necessarily standardized
        absent_cols = ["SMILES"] + RDCOLS + SILLY_COLS + SIFT_COLS
        expected = self.expected_filtered_df.drop(columns=absent_cols)
        pd.testing.assert_frame_equal(filtered_df.drop(columns="SMILES"), expected)

    def test_sillyfilterMols(self):
        filtered_df = self.coreFilter._sillyfilterMols(self.test_mols)
        # Won't test for the smiles to be the same - not necessarily standardized
        absent_cols = ["SMILES"] + RDCOLS + BLOOM_COLS + SIFT_COLS
        expected = self.expected_filtered_df.drop(columns=absent_cols)
        pd.testing.assert_frame_equal(filtered_df.drop(columns="SMILES"), expected)

    def tearDown(self) -> None:
        return super().tearDown()


if __name__ == "__main__":
    unittest.main()
