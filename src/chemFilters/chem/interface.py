from typing import Union
from multiprocessing import Pool
from typing import Iterable, Callable
from loguru import logger
from rdkit import Chem
from tqdm import tqdm


def mol_from_smi(smi: str):
    """Convert a SMILES string to a rdkit.Chem.Mol object."""
    if smi is None:
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        logger.warning(f"Could not convert SMILES {smi} to rdkit.Chem.Mol")
    return mol


def mol_to_smi(mol, **kwargs):
    """Convert a rdkit.Chem.Mol object to a SMILES string."""
    try:
        smi = Chem.MolToSmiles(mol, **kwargs)
    except TypeError as e:
        logger.warning(f"Exception!! {e}\n" f"Could not convert {mol} to SMILES.")
        smi = None
    return smi


class MoleculeHandler:
    """Interface class for handling molecules. implemented so I can use this `from_smi`
    functionalitiy on other classes in the package.
    """

    def __init__(self, from_smi=False, isomeric=True) -> None:
        self._from_smi = from_smi
        self._return_isomeric = isomeric

    def _output_mol(self, stdin: Union[str, Chem.Mol]):
        """Treat `stdin` as [str|Mol] depending on self.from_smiles and return a Mol
        or a None in case the SMILES are invalid."""
        if self._from_smi:
            if stdin is not None and not isinstance(stdin, str):
                raise ValueError(
                    "If from_smi is True, inputs must be None or SMILES strings."
                )
            else:
                mol = mol_from_smi(stdin)
        else:
            if isinstance(stdin, str):
                raise ValueError(
                    "If from_smi is False, inputs must be rdkit.Chem.Mol objects."
                )
            elif stdin is None or isinstance(stdin, Chem.Mol):
                mol = stdin
            else:
                raise ValueError(
                    f"stdin of type {type(stdin)} not recognized. "
                    "Must be None, or rdkit.Chem.Mol object."
                )
        return mol

    def _output_smi(self, stdin: Union[str, Chem.Mol]):
        """Treat `stdin` as [str|Mol] depending on self.from_smiles and return a SMILES
        or a None in case the Mol is invalid."""
        if self._from_smi:
            if isinstance(stdin, Chem.Mol):
                raise ValueError(
                    "If from_smi is True, inputs must be rdkit.Chem.Mol objects."
                )
            elif any([stdin is None, isinstance(stdin, str)]):
                if mol_from_smi(stdin) is None:
                    smi = None
                else:
                    smi = stdin
            else:
                raise ValueError(
                    f"stdin of type {type(stdin)} not recognized. "
                    "Must be None, or SMILES string."
                )
        else:
            if stdin is None:
                smi = None
            elif isinstance(stdin, Chem.Mol):
                smi = mol_to_smi(
                    stdin,
                    kekuleSmiles=False,
                    canonical=True,
                    isomericSmiles=self._return_isomeric,
                )
            else:
                raise ValueError(
                    "If from_smi is False, inputs must be None or SMILES strings."
                )
        return smi

    def pmap(
        self,
        n_jobs: int,
        progress: bool,
        stdin: Iterable,
        func: Callable,
        pickable: bool = True,
    ):
        """Helper function to map a function to an iterable using a multiprocessing pool.

        Args:
            n_jobs: number of jobs for the pool.
            progress: display progress bar with tqdm.
            stdin: iterable to map the function to.
            func: function to be mapped to the variables.
            pickable: bool indicating whether it should parallel process or not.
                Defaults to True.

        Returns:
            A list of the results of the function mapped to the iterable.
        """  # TODO: convert this to use joblib instead
        if pickable:
            with Pool(n_jobs) as p:
                if progress:
                    vals = list(tqdm(p.imap(func, stdin), total=len(stdin)))
                else:
                    vals = p.map(func, stdin)
        else:
            if progress:
                vals = [func(_in) for _in in stdin]
                list(tqdm(map(func, stdin), total=len(stdin)))
            else:
                vals = list(map(func, stdin))
        return vals
