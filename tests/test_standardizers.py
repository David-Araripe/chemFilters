import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.MolStandardize.rdMolStandardize import ChargeParent

from chemFilters.chem.standardizers import ChemStandardizer


def pickable_std_smi_func(smi):
    mol = Chem.MolFromSmiles(smi)
    ChargeParent(mol)
    return mol


def pickable_std_mol_func(mol):
    ChargeParent(mol)
    return mol


class TestChemStandardizer(unittest.TestCase):
    def setUp(self):
        self.testroot = Path(__file__).parent
        self.aspirin = "CC(=O)Oc1ccccc1C(=O)O"
        self.chembl_smi_standardizer = ChemStandardizer(method="chembl", from_smi=True)
        self.chembl_mol_standardizer = ChemStandardizer(method="chembl", from_smi=False)
        self.molvs_smi_standardizer = ChemStandardizer(method="molvs", from_smi=True)
        self.custom_smi_standardizer = ChemStandardizer(
            method=pickable_std_smi_func, from_smi=True
        )
        self.custom_mol_standardizer = ChemStandardizer(
            method=pickable_std_mol_func, from_smi=False
        )
        self.papyrus_smi_standardizer = ChemStandardizer(
            method="papyrus", from_smi=True, small_molecule_min_mw=150
        )
        self.test_smiles = (
            (self.testroot / "resources/testSmiles.smi").read_text().splitlines()
        )

    def test_custom_standardizer(self):
        result = self.custom_smi_standardizer([self.aspirin] * 5)
        self.assertEqual(result, ["CC(=O)Oc1ccccc1C(=O)O"] * 5)
        result = self.custom_mol_standardizer([Chem.MolFromSmiles(self.aspirin)] * 5)
        self.assertEqual(result, ["CC(=O)Oc1ccccc1C(=O)O"] * 5)
        nopickle_custom_from_smi = ChemStandardizer(
            method=lambda smi: Chem.MolFromSmiles(smi), from_smi=True
        )
        nopicle_custom_from_mol = ChemStandardizer(
            method=lambda mol: ChargeParent(mol), from_smi=False
        )
        result = nopickle_custom_from_smi([self.aspirin] * 5)
        self.assertEqual(result, ["CC(=O)Oc1ccccc1C(=O)O"] * 5)
        result = nopicle_custom_from_mol([Chem.MolFromSmiles(self.aspirin)] * 5)
        self.assertEqual(result, ["CC(=O)Oc1ccccc1C(=O)O"] * 5)

    def test_canon_smiles(self):
        stdzer = ChemStandardizer(method="canon", from_smi=True)
        result = stdzer([self.aspirin] * 5)
        self.assertEqual(result, ["CC(=O)Oc1ccccc1C(=O)O"] * 5)
        stdzer = ChemStandardizer(method="canon", from_smi=False)
        result = stdzer([Chem.MolFromSmiles(self.aspirin)] * 5)
        self.assertEqual(result, ["CC(=O)Oc1ccccc1C(=O)O"] * 5)

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

    def test_molvs_standardizer(self):
        mols = [Chem.MolFromSmiles(smi) for smi in self.test_smiles]
        std_smiles = self.molvs_smi_standardizer(self.test_smiles)
        for smi in std_smiles:
            self.assertIsInstance(smi, str)
        std_smiles = self.chembl_mol_standardizer(mols)
        for smi in std_smiles:
            self.assertIsInstance(smi, str)
        with self.assertRaises(ValueError):
            self.chembl_mol_standardizer(self.test_smiles)

    def tearDown(self):
        pass


if __name__ == "__main__":
    unittest.main()
