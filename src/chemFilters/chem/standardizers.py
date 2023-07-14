# -*- coding: utf-8 -*-
"""Utility modules with functions for the chemFilters.chem subpackage."""
import warnings
from functools import partial
from importlib.util import find_spec
from multiprocessing import Pool
from typing import List, Union

from chembl_structure_pipeline import standardizer as chembl_std
from rdkit import Chem
from tqdm import tqdm

from .utils import (
    RDKitVerbosityOFF,
    RDKitVerbosityON,
    molToCanon,
    molToConnectivity,
    molToInchi,
    molToInchiKey,
)

if find_spec("papyrus_structure_pipeline"):
    from papyrus_structure_pipeline import standardizer as papyrus_std

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
    
class ChemStandardizer(MoleculeHandler):
    """A class to standardize molecules/SMILES strings. Initialization allows for the 
    selection of the settings of the standardizer. The object can then be called on a
    iterable containing molecules/SMILES strings to apply the standardization."""

    def __init__(
        self,
        method: str = "chembl",
        n_jobs: int = 5,
        isomeric: bool = True,
        progress: bool = True,
        verbose: bool = False,
        from_smi: bool = True,
    ) -> None:
        """Initializes the ChemStandardizer class.

        Args:
            method: which standardization pipeline to use. Current supports "canon",
                "chembl" and "papyrus", but the latter is an optional dependency.
                Defaults to "chembl". "canon" is rdkit's SMILES canonicalization.
            n_jobs: number of jobs running in parallel. Defaults to 5.
            isomeric: output smiles with isomeric information. Defaults to True.
            progress: display a progress bar with tqdm. Defaults to True.
            verbose: if false, silences rdkit errors. Defaults to False.
            from_smi: if True, the standardizer will expect SMILES strings as input.

        Raises:
            ImportError: if method is "papyrus" but the optional dependency is not
                installed.
            ValueError: when an invalid method is passed.
        """
        if method.lower() == "chembl":
            self.standardizer = partial(
                self.chemblStandardizer, isomeric=isomeric
            )
        elif method.lower() == "canon":
            self.standardizer = partial(molToCanon, isomeric=isomeric)
        elif method.lower() == "papyrus":
            # avoid import since it's not a required dependency
            if find_spec("papyrus_structure_pipeline"):
                self.standardizer = partial(
                    self.papyrusStandardizer, isomeric=isomeric
                )
            else:
                raise ImportError(
                    "Optional dependency not found. Please install it by running:\n"
                    "python -m pip install "
                    "git+https://github.com/OlivierBeq/Papyrus_structure_pipeline.git"
                )
        else:
            raise ValueError(f"Invalid SMILES standardizing method: {method}")
        super().__init__(from_smi)
        self.n_jobs = n_jobs
        self.progress = progress
        self.verbose = verbose

    def __call__(self, stdin: List[Union[str, Chem.Mol]]) -> List[str]:
        """Calls the standardizer on a list of SMILES strings to perform the
        standardization according to the settings set at initialization.

        Args:
            stdin: standard input; a list of SMILES strings or rdkit.Chem.Mol objects
                depending on the value of self.from_smi.

        Returns:
            A list of standardized SMILES strings.
        """
        if not self.verbose:
            RDKitVerbosityOFF()
        with Pool(self.n_jobs) as p:
            if self.progress:
                vals = list(tqdm(p.imap(self.standardizer, stdin), total=len(stdin)))
            else:
                vals = p.map(self.standardizer, stdin)
        # restore the logger level
        RDKitVerbosityON()
        return vals
    
    def papyrusStandardizer(
        self,
        stdin: Union[str, Chem.Mol],
        isomeric: bool = True,
        **kwargs,
    ) -> str:
        """Uses the Papyrus standardizer to standardize a SMILES string. By default,
        this standardization pipeline removes stereocenters, so beware of the isomeric
        flag. Accepts extra keyword arguments that will be passed to the standardizer.

        For more information: https://github.com/OlivierBeq/Papyrus_structure_pipeline

        Args:
            smi: single smiles string.
            isomeric: output isomerc smiles. Defaults to True.

        Returns:
            standardized smiles string
        """
        mol = self._mol_handler(stdin)
        try:
            standard_mol = papyrus_std.standardize(mol, **kwargs)
        except RuntimeError:
            print("Error standardizing molecule: ", stdin)
            standard_mol = None
        if standard_mol is None:
            return None
        standard_smi = Chem.MolToSmiles(
            standard_mol, kekuleSmiles=False, canonical=True, isomericSmiles=isomeric
        )
        return standard_smi

    def chemblStandardizer(self, stdin: str, isomeric: bool = True, **kwargs) -> str:
        """Uses the ChEMBL standardizer to standardize a SMILES string. Accepts extra
        keyword arguments that will be passed to the standardizer

        Args:
            smi: single smiles string.
            isomeric: output isomeric smiles. Defaults to True.

        Returns:
            standardized smiles string
        """
        mol = self._mol_handler(stdin)
        standard_mol = chembl_std.standardize_mol(mol, **kwargs)
        standard_smi = Chem.MolToSmiles(
            standard_mol, kekuleSmiles=False, canonical=True, isomericSmiles=isomeric
        )
        return standard_smi


class InchiHandling(MoleculeHandler):
    """Obtain a list of inchis, inchikeys or connectivities from a list of smiles.
    Initialization allows for the selection of the settings. The object can then be
    called on a iterable containing SMILES strings to obtain the desired identifier."""

    def __init__(
        self,
        convert_to: str,
        n_jobs: int = 5,
        progress: bool = True,
        verbose: bool = False,
        from_smi: bool = True,
    ) -> None:
        """Initialize the InchiHandling class.

        Args:
            convert_to: what to convert the smiles to. Can be "inchi", "inchikey" or
                "connectivity".
            n_jobs: Number of jobs for processing in parallel. Defaults to 5.
            progress: whether to show the progress bar. Defaults to True.
            verbose: if false, will hide the rdkit warnings. Defaults to False.

        Raises:
            ValueError: if the convert_to argument is not one of the three options.
        """
        if convert_to.lower() == "inchi":
            self.converter = partial(molToInchi, verbose=verbose)
        elif convert_to.lower() == "inchikey":
            self.converter = partial(molToInchiKey, verbose=verbose)
        elif convert_to.lower() == "connectivity":
            self.converter = partial(molToConnectivity, verbose=verbose)
        else:
            raise ValueError(f"Invalid convertion method: {self.convert_to}")
        self.n_jobs = n_jobs
        self.progress = progress
        self.verbose = verbose
        super().__init__(from_smi)

    def __call__(self, stdin: list) -> list:
        if not self.verbose:
            RDKitVerbosityOFF()
        with Pool(self.n_jobs) as p:
            mols = p.map(self._mol_handler, stdin)
            if self.progress:
                vals = list(tqdm(p.imap(self.converter, mols), total=len(mols)))
            else:
                vals = p.map(self.converter, mols)
        # restore the logger level
        RDKitVerbosityON()
        return vals
