import unittest

from rdkit import Chem

from chemFilters.chem.interface import MoleculeHandler


class TestMoleculeHandler(unittest.TestCase):
    def setUp(self):
        self.mol_handler_from_smi = MoleculeHandler(from_smi=True)
        self.mol_handler_from_mol = MoleculeHandler(from_smi=False)
        self.valid_smi = "CC(=O)Oc1ccccc1C(=O)O"
        self.invalid_smi = "invalid 😤"
        self.valid_mol = Chem.MolFromSmiles(self.valid_smi)

    def test_output_mol(self):
        # Test Molecule -> Molecule
        self.assertEqual(
            self.mol_handler_from_mol._output_mol(self.valid_mol), self.valid_mol
        )
        # Test None -> None
        self.assertIsNone(self.mol_handler_from_mol._output_mol(None))

        # Test Valid SMILES -> Molecule
        self.assertTrue(
            isinstance(self.mol_handler_from_smi._output_mol(self.valid_smi), Chem.Mol)
        )
        # Test Invalid SMILES -> None
        self.assertIsNone(self.mol_handler_from_smi._output_mol(self.invalid_smi))
        # Test wrong input -> error
        with self.assertRaises(ValueError):
            self.mol_handler_from_smi._output_mol(self.valid_mol)
        with self.assertRaises(ValueError):
            self.mol_handler_from_mol._output_mol(self.valid_smi)
        with self.assertRaises(ValueError):
            self.mol_handler_from_mol._output_mol(True)

    def test_output_smi(self):
        # Test Molecule -> SMILES
        self.assertEqual(
            self.mol_handler_from_mol._output_smi(self.valid_mol), self.valid_smi
        )
        # Test None -> None
        self.assertIsNone(self.mol_handler_from_mol._output_smi(None))
        # Test Valid SMILES -> SMILES
        self.assertEqual(
            self.mol_handler_from_smi._output_smi(self.valid_smi), self.valid_smi
        )
        # Test Invalid SMILES -> None
        self.assertIsNone(self.mol_handler_from_smi._output_smi(self.invalid_smi))

        # Test wrong input -> error
        with self.assertRaises(ValueError):
            self.mol_handler_from_smi._output_smi(self.valid_mol)
        with self.assertRaises(ValueError):
            self.mol_handler_from_mol._output_smi(self.valid_smi)


if __name__ == "__main__":
    unittest.main()
