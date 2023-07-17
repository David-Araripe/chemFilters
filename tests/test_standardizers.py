import unittest
from pathlib import Path

from rdkit import Chem

from chemFilters.chem.standardizers import ChemStandardizer


class TestChemStandardizer(unittest.TestCase):
    def setUp(self):
        self.testroot = Path(__file__).parent
        self.aspirin = "CC(=O)Oc1ccccc1C(=O)O"
        self.chembl_smi_standardizer = ChemStandardizer(method="chembl", from_smi=True)
        self.chembl_mol_standardizer = ChemStandardizer(method="chembl", from_smi=False)
        self.papyrus_smi_standardizer = ChemStandardizer(
            method="papyrus", from_smi=True, small_molecule_min_mw=150
        )
        self.test_smiles = (
            (self.testroot / "resources/testSmiles.smi").read_text().splitlines()
        )

    def test_standardizer_invalid_smiles(self):
        result = self.chembl_smi_standardizer(["Got", "some", "invalid", "smiles"])
        self.assertEqual(result, [None, None, None, None])

    def test_standardizer_none_smiles(self):
        result = self.chembl_smi_standardizer([self.aspirin, None])
        self.assertEqual(result, [self.aspirin, None])

    def test_standardizer_invalid_input(self):
        with self.assertRaises(ValueError):
            self.chembl_smi_standardizer([True, "that"])

    def test_standardizer_invalid_method(self):
        with self.assertRaises(ValueError):
            ChemStandardizer(method="invalid_method", from_smi=True)

    def test_chembl_standardizer(self):
        mols = [Chem.MolFromSmiles(smi) for smi in self.test_smiles]
        std_smiles = self.chembl_smi_standardizer(self.test_smiles)
        for smi in std_smiles:
            self.assertIsInstance(smi, str)
        std_smiles = self.chembl_mol_standardizer(mols)
        for smi in std_smiles:
            self.assertIsInstance(smi, str)
        with self.assertRaises(ValueError):
            self.chembl_mol_standardizer(self.test_smiles)

    def test_papyrus_standardizer(self):
        result = self.papyrus_smi_standardizer([self.aspirin])
        self.assertEqual(result, [self.aspirin])

    def tearDown(self):
        pass


if __name__ == "__main__":
    unittest.main()
