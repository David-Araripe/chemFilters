import unittest

import PIL
from rdkit import Chem

from chemFilters.img_render import FontManager, MolPlotter


class TestFontManager(unittest.TestCase):
    def setUp(self) -> None:
        self.fm = FontManager()
        self.wsl_fm = FontManager(wsl="auto")

    def test_operating_system(self):
        op_sys = self.fm.operating_system
        # Check that the returned operating system is a string
        self.assertIsInstance(op_sys, str)
        # Check that the returned operating system is one of the expected values
        self.assertIn(op_sys, ["windows", "wsl", "linux", "macos", "Unknown"])

    def test_available_fonts(self):
        fonts = self.fm.available_fonts
        # Check that the returned fonts is a dictionary
        self.assertIsInstance(fonts, dict)
        # Optionally, you can check that the dictionary is not empty
        self.assertNotEqual(len(fonts), 0)

    def tearDown(self) -> None:
        return super().tearDown()


class TestMolPlotter(unittest.TestCase):
    def setUp(self):
        """Define the setup method."""
        self.plotter = MolPlotter(from_smi=False)
        self.molecule = Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O")

    def test_available_fonts(self):
        self.assertIsInstance(self.plotter.available_fonts, dict)
        self.assertTrue(len(self.plotter.available_fonts) > 0)

    def test_process_mols(self):
        mol_list = [
            Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O"),
            Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O"),
        ]
        result = self.plotter._process_mols(mol_list)
        self.assertEqual(len(result), 2)
        self.assertIsInstance(result[0], Chem.rdchem.Mol)

    def test_substructure_palette(self):
        result = self.plotter._substructure_palette(3, "rainbow", 1.0)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 3)

    def test_render_mol(self):
        result = self.plotter.render_mol(self.molecule, label="test_label")
        self.assertIsInstance(result, PIL.Image.Image)

    def test_render_with_matches(self):
        substructs = [Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O")]
        result = self.plotter.render_with_matches(
            self.molecule, substructs, label="test_label"
        )
        self.assertIsInstance(result, PIL.Image.Image)

    def test_render_with_colored_matches(self):
        descriptions = ["test_description"]
        substructs = [Chem.MolFromSmiles("C(=O)O")]
        result = self.plotter.render_with_colored_matches(
            self.molecule, descriptions, substructs, label="test_label"
        )
        self.assertIsInstance(result, PIL.Image.Image)

    def tearDown(self):
        """Define the teardown method."""
        del self.plotter
        del self.molecule


if __name__ == "__main__":
    unittest.main()
