import logging
from typing import Union

from rdkit import Chem


class MoleculeHandler:
    """Interface class for handling molecules. implemented so I can use this `from_smi`
    functionalitiy on other classes in the package.
    """

    def __init__(self, from_smi=False) -> None:
        self._from_smi = from_smi

    def _mol_handler(self, stdin: Union[str, Chem.Mol]):
        """Treat `stdin` as [str|Mol] depending on self.from_smiles and return a Mol
        or a None in case the SMILES are invalid."""
        if self._from_smi:
            if stdin is not None and not isinstance(stdin, str):
                raise ValueError(
                    "If from_smi is True, inputs must be None or SMILES strings."
                )
            try:
                mol = Chem.MolFromSmiles(stdin)
            except TypeError:
                logging.warn(f"Error converting SMILES to Mol: {stdin}")
                mol = None
        else:
            if isinstance(stdin, str):
                raise ValueError(
                    "If from_smi is False, inputs must be rdkit.Chem.Mol objects."
                )
            mol = stdin
        return mol
