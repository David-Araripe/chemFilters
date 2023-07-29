# -*- coding: utf-8 -*-
"""Utility modules with functions for the chemFilters.chem subpackage."""
import logging
from functools import partial
from importlib.util import find_spec
from multiprocessing import Pool
from typing import List, Union

from chembl_structure_pipeline import standardizer as chembl_std
from rdkit import Chem
from tqdm import tqdm

from .interface import MoleculeHandler
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


class ChemStandardizer(MoleculeHandler):
    """A class to standardize molecules/SMILES strings. Initialization allows for the
    selection of the settings of the standardizer. The object can then be called on a
    iterable containing molecules/SMILES strings to apply the standardization."""

    def __init__(
        self,
        method: str = "chembl",
        n_jobs: int = 5,
        isomeric: bool = True,
        progress: bool = False,
        verbose: bool = True,
        from_smi: bool = False,
        return_smi: bool = True,
        **kwargs,
    ) -> None:
        """Initializes the ChemStandardizer class.

        Args:
            method: which standardization pipeline to use. Current supports "canon",
                "chembl" and "papyrus", but the latter is an optional dependency.
                Defaults to "chembl". "canon" is rdkit's SMILES canonicalization.
            n_jobs: number of jobs running in parallel. Defaults to 5.
            isomeric: output smiles with isomeric information. Defaults to True.
            progress: display a progress bar with tqdm. Defaults to False.
            verbose: if false, silences rdkit errors. Defaults to True.
            from_smi: if True, the standardizer will expect SMILES strings as input.
                Defaults to False.
            kwargs: additional keyword arguments to pass to the standardizer.

        Raises:
            ImportError: if method is "papyrus" but the optional dependency is not
                installed.
            ValueError: when an invalid method is passed.
        """
        if method.lower() == "chembl":
            self.standardizer = partial(self.chemblStandardizer, **kwargs)
        elif method.lower() == "canon":
            self.standardizer = partial(molToCanon, isomeric=isomeric, **kwargs)
        elif method.lower() == "papyrus":
            # avoid import since it's not a required dependency
            if find_spec("papyrus_structure_pipeline"):
                self.standardizer = partial(self.papyrusStandardizer, **kwargs)
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
        self.return_smi = return_smi
        self.method = method
        self.smi_out_func = partial(
            MoleculeHandler(from_smi=False, isomeric=isomeric)._output_smi
        )

    def __call__(self, stdin: List[Union[str, Chem.Mol]]) -> List[str]:
        """Calls the standardizer on a list of SMILES strings / Chem.Mol objects to
        perform the standardization according to the settings set at initialization.

        Args:
            stdin: standard input; a list of SMILES strings or rdkit.Chem.Mol objects
                depending on the value of self._from_smi.

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
            if self.return_smi and self.method != "canon":  # canon always returns smis
                if self.progress:
                    vals = list(tqdm(p.imap(self.smi_out_func, vals), total=len(vals)))
                else:
                    vals = p.map(self.smi_out_func, vals)
        RDKitVerbosityON()  # restore the logger level
        return vals

    def papyrusStandardizer(
        self,
        stdin: Union[str, Chem.Mol],
        **kwargs,
    ) -> str:
        """Uses the Papyrus standardizer to standardize a SMILES string. By default,
        this standardization pipeline removes stereocenters, so beware of the isomeric
        flag. Accepts extra keyword arguments that will be passed to the standardizer.

        For more information: https://github.com/OlivierBeq/Papyrus_structure_pipeline

        Args:
            stdin: standard input; single SMILES strings or single rdkit.Chem.Mol object
                depending on the value of self._from_smi.
            isomeric: output isomerc smiles. Defaults to True.
            kwargs: aditional keyword arguments to pass to the standardizer.

        Returns:
            standardized smiles string
        """
        mol = self._output_mol(stdin)
        try:
            standard_mol = papyrus_std.standardize(mol, **kwargs)
        except RuntimeError:
            logging.exception("Error standardizing molecule: ", stdin)
            standard_mol = None
        return standard_mol

    def chemblStandardizer(self, stdin: Union[str, Chem.Mol], **kwargs) -> str:
        """Uses the ChEMBL standardizer to standardize a SMILES string. Accepts extra
        keyword arguments that will be passed to the standardizer

        Args:
            stdin: standard input; single SMILES strings or single rdkit.Chem.Mol object
                depending on the value of self._from_smi.
            isomeric: output isomeric smiles. Defaults to True.
            kwargs: additional keyword arguments to pass to the standardizer.

        Returns:
            standardized smiles string
        """
        mol = self._output_mol(stdin)
        if mol is None:
            return None
        try:
            standard_mol = chembl_std.standardize_mol(mol, sanitize=True, **kwargs)
        except Exception as e:
            logging.exception("Error standardizing molecule: ", stdin)
            logging.exception(e)
            standard_mol = None
        return standard_mol


class InchiHandling(MoleculeHandler):
    """Obtain a list of inchis, inchikeys or connectivities from a list of smiles.
    Initialization allows for the selection of the settings. The object can then be
    called on a iterable containing SMILES strings to obtain the desired identifier."""

    def __init__(
        self,
        convert_to: str,
        n_jobs: int = 5,
        progress: bool = False,
        verbose: bool = True,
        from_smi: bool = False,
    ) -> None:
        """Initialize the InchiHandling class.

        Args:
            convert_to: what to convert the smiles to. Can be "inchi", "inchikey" or
                "connectivity".
            n_jobs: Number of jobs for processing in parallel. Defaults to 5.
            progress: whether to show the progress bar. Defaults to False.
            verbose: if false, will hide the rdkit warnings. Defaults to True.
            from_smi: if True, the standardizer will expect SMILES strings as input.
                Defaults to False.

        Raises:
            ValueError: if the convert_to argument is not one of the three options.
        """
        if convert_to.lower() == "inchi":
            self.converter = molToInchi
        elif convert_to.lower() == "inchikey":
            self.converter = molToInchiKey
        elif convert_to.lower() == "connectivity":
            self.converter = molToConnectivity
        else:
            raise ValueError(f"Invalid convertion method: {self.convert_to}")
        self.n_jobs = n_jobs
        self.progress = progress
        self.verbose = verbose
        super().__init__(from_smi)

    def __call__(self, stdin: list) -> list:
        """Calls the standardizer on a list of SMILES strings / Chem.Mol objects to
        perform the standardization according to the settings set at initialization.

        Args:
            stdin: standard input; a list of SMILES strings or rdkit.Chem.Mol objects
                depending on the value of self._from_smi.

        Returns:
            A list of standardized SMILES strings.
        """
        if not self.verbose:
            RDKitVerbosityOFF()
        with Pool(self.n_jobs) as p:
            mols = p.map(self._output_mol, stdin)
            if self.progress:
                vals = list(tqdm(p.imap(self.converter, mols), total=len(mols)))
            else:
                vals = p.map(self.converter, mols)
        # restore the logger level
        RDKitVerbosityON()
        return vals
