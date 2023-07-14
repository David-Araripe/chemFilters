import warnings
from typing import Union

from rdkit import Chem


class MoleculeHandler:
    """Interface clas for handling molecules. implemented so I can use this `from_smi`
    functionalitiy on other classes used in the pipeline.
    """

    def __init__(self, from_smi) -> None:
        self.from_smi = from_smi

    def _mol_handler(self, stdin: Union[str, Chem.Mol]):
        """Treat `stdin` as [str|Mol] depending on self.from_smiles and return a Mol
        or a None in case the SMILES are invalid."""
        if self.from_smi:
            try:
                mol = Chem.MolFromSmiles(stdin)
            except TypeError:
                warnings.warn(f"Error converting SMILES to Mol: {stdin}")
                mol = None
        else:
            if isinstance(stdin, str):
                raise ValueError(
                    "If from_smi is False, inputs must be rdkit.Chem.Mol objects."
                )
            mol = stdin
        return mol
