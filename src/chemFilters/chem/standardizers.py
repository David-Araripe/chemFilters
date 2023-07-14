# -*- coding: utf-8 -*-
"""Utility modules with functions for the chemFilters.chem subpackage."""
from functools import partial
from importlib.util import find_spec
from multiprocessing import Pool
from typing import List

from chembl_structure_pipeline import standardizer as chembl_std
from rdkit import Chem
from tqdm import tqdm

from .utils import (
    RDKitVerbosityOFF,
    RDKitVerbosityON,
    smilesToCanon,
    smilesToConnectivity,
    smilesToInchi,
    smilesToInchiKey,
)

if find_spec("papyrus_structure_pipeline"):
    from papyrus_structure_pipeline import standardizer as papyrus_std


class SmilesStandardizer:
    """A class to standardize SMILES strings. Initialization allows for the selection
    of the settings of the standardizer. The object can then be called on a iterable
    containing SMILES strings to make the standardization."""

    def __init__(
        self,
        method: str = "chembl",
        n_jobs: int = 5,
        isomeric: bool = True,
        progress: bool = True,
        verbose: bool = False,
    ) -> None:
        """Initializes the SmilesStandardizer class.

        Args:
            method: which standardization pipeline to use. Current supports "canon",
                "chembl" and "papyrus", but the latter is an optional dependency.
                Defaults to "chembl". "canon" is rdkit's SMILES canonicalization.
            n_jobs: number of jobs running in parallel. Defaults to 5.
            isomeric: output smiles with isomeric information. Defaults to True.
            progress: display a progress bar with tqdm. Defaults to True.
            verbose: if false, silences rdkit errors. Defaults to False.

        Raises:
            ImportError: if method is "papyrus" but the optional dependency is not
                installed.
            ValueError: when an invalid method is passed.
        """
        if method.lower() == "chembl":
            self.standardizer = partial(
                self.chemblSmilesStandardizer, isomeric=isomeric
            )
        elif method.lower() == "canon":
            self.standardizer = partial(smilesToCanon, isomeric=isomeric)
        elif method.lower() == "papyrus":
            # avoid import since it's not a required dependency
            if find_spec("papyrus_structure_pipeline"):
                self.standardizer = partial(
                    self.papyrusSmilesStandardizer, isomeric=isomeric
                )
            else:
                raise ImportError(
                    "Optional dependency not found. Please install it by running:\n"
                    "python -m pip install "
                    "git+https://github.com/OlivierBeq/Papyrus_structure_pipeline.git"
                )
        else:
            raise ValueError(f"Invalid SMILES standardizing method: {method}")
        self.n_jobs = n_jobs
        self.progress = progress
        self.verbose = verbose

    def __call__(self, smiles: List[str]) -> List[str]:
        """Calls the standardizer on a list of SMILES strings to perform the
        standardization according to the settings set at initialization.

        Args:
            smiles: list of smiles strings.

        Returns:
            _description_
        """
        if not self.verbose:
            RDKitVerbosityOFF()
        with Pool(self.n_jobs) as p:
            if self.progress:
                vals = list(tqdm(p.imap(self.standardizer, smiles), total=len(smiles)))
            else:
                vals = p.map(self.standardizer, smiles)
        # restore the logger level
        RDKitVerbosityON()
        return vals

    @staticmethod
    def papyrusSmilesStandardizer(
        smi: str,
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
        try:
            mol = Chem.MolFromSmiles(smi)
        except TypeError:
            return None
        try:
            standard_mol = papyrus_std.standardize(mol, **kwargs)
        except RuntimeError:
            print("Error standardizing molecule: ", smi)
            standard_mol = None
        if standard_mol is None:
            return None
        standard_smi = Chem.MolToSmiles(
            standard_mol, kekuleSmiles=False, canonical=True, isomericSmiles=isomeric
        )
        return standard_smi

    @staticmethod
    def chemblSmilesStandardizer(smi: str, isomeric: bool = True, **kwargs) -> str:
        """Uses the ChEMBL standardizer to standardize a SMILES string. Accepts extra
        keyword arguments that will be passed to the standardizer

        Args:
            smi: single smiles string.
            isomeric: output isomeric smiles. Defaults to True.

        Returns:
            standardized smiles string
        """
        try:
            mol = Chem.MolFromSmiles(smi)
        except TypeError:
            return None
        standard_mol = chembl_std.standardize_mol(mol, **kwargs)
        standard_smi = Chem.MolToSmiles(
            standard_mol, kekuleSmiles=False, canonical=True, isomericSmiles=isomeric
        )
        return standard_smi


class InchiHandling:
    """Obtain a list of inchis, inchikeys or connectivities from a list of smiles.
    Initialization allows for the selection of the settings. The object can then be
    called on a iterable containing SMILES strings to obtain the desired identifier."""

    def __init__(
        self,
        convert_to: str,
        n_jobs: int = 5,
        progress: bool = True,
        verbose: bool = False,
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
            self.converter = partial(smilesToInchi, verbose=verbose)
        elif convert_to.lower() == "inchikey":
            self.converter = partial(smilesToInchiKey, verbose=verbose)
        elif convert_to.lower() == "connectivity":
            self.converter = partial(smilesToConnectivity, verbose=verbose)
        else:
            raise ValueError(f"Invalid convertion method: {self.convert_to}")
        self.n_jobs = n_jobs
        self.progress = progress
        self.verbose = verbose

    def __call__(self, smiles: list) -> list:
        if not self.verbose:
            RDKitVerbosityOFF()
        with Pool(self.n_jobs) as p:
            if self.progress:
                vals = list(tqdm(p.imap(self.converter, smiles), total=len(smiles)))
            else:
                vals = p.map(self.converter, smiles)
        # restore the logger level
        RDKitVerbosityON()
        return vals
