# -*- coding: utf-8 -*-
"""Utility modules with functions for the chemFilters.chem subpackage."""
from functools import partial
from importlib.util import find_spec
from typing import Callable, List, Union

from chembl_structure_pipeline import standardizer as chembl_std
from rdkit import Chem

from ..logger import logger
from .interface import MoleculeHandler
from .utils import (
    molToCanon,
    molToConnectivity,
    molToInchi,
    molToInchiKey,
    rdkit_log_controller,
)

if find_spec("papyrus_structure_pipeline"):
    from papyrus_structure_pipeline import standardizer as papyrus_std

if find_spec("molvs"):
    from molvs import Standardizer


class ChemStandardizer(MoleculeHandler):
    """A class to standardize molecules/SMILES strings. Initialization allows for the
    selection of the settings of the standardizer. The object can then be called on a
    iterable containing molecules/SMILES strings to apply the standardization."""

    def __init__(
        self,
        method: Union[str, Callable] = "chembl",
        n_jobs: int = 5,
        isomeric: bool = True,
        progress: bool = False,
        rdkit_loglevel: str = "warning",
        from_smi: bool = False,
        return_smi: bool = True,
        **kwargs,
    ) -> None:
        """Initializes the ChemStandardizer class.

        Args:
            method: standardization pipeline to use. Current supports "canon",
                "chembl", "papyrus", "molvs", or a callable. If callable, ensure it
                takes rdkit.Mol objects as input. Defaults to "chembl". "canon" is
                rdkit's SMILES canonicalization.
            n_jobs: number of jobs running in parallel. Defaults to 5.
            isomeric: output smiles with isomeric information. Defaults to True.
            progress: display a progress bar with tqdm. Defaults to False.
            rdkit_loglevel: one of `debug, info, warning, error, critical`. Defaults to
                "warning".
            from_smi: if True, the standardizer will expect SMILES strings as input.
                Defaults to False.
            kwargs: additional keyword arguments to pass to the standardizer.

        Raises:
            ImportError: if method is "papyrus" but the optional dependency is not
                installed.
            ValueError: when an invalid method is passed.
        """
        if callable(method):
            self.standardizer = method
            logger.info(
                "Custom standardizer detected!! "
                "Ensure it takes Rdkit.Mol objecs as input."
            )
            try:
                import pickle

                _ = pickle.dumps(self.standardizer)
                self._custom_is_pickable = True
            except (pickle.PicklingError, TypeError, AttributeError):
                logger.warning(
                    "Custom standardizer is not pickable. Will use 'n_jobs' "
                    "only to process `SMILES -> Mol` & `Mol -> SMILES` operations"
                )
                self._custom_is_pickable = False
        else:
            self._custom_is_pickable = False
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
                        "python -m pip install git+"
                        "https://github.com/OlivierBeq/Papyrus_structure_pipeline.git"
                    )
            elif method.lower() == "molvs":
                if find_spec("molvs"):
                    self.standardizer = partial(self.molvsStandardizer, **kwargs)
                else:
                    raise ImportError(
                        "Optional dependency not found. Please install it by running:\n"
                        "python -m pip install molvs"
                    )
            else:
                raise ValueError(f"Invalid SMILES standardizing method: {method}")
        super().__init__(from_smi)
        self.n_jobs = n_jobs
        self.progress = progress
        self.rdkit_loglevel = rdkit_loglevel
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
        if callable(self.method):
            if self.method.__name__ == "<lambda>":
                pickable = False
            else:
                pickable = self._custom_is_pickable
        else:
            pickable = True
        with rdkit_log_controller(self.rdkit_loglevel):
            if self._from_smi and self.method == "canon":
                stdin = self.pmap(self.n_jobs, self.progress, stdin, self._output_mol)
            print("NOW PROCESSING")
            vals = self.pmap(
                self.n_jobs,
                self.progress,
                stdin,
                self.standardizer,
                pickable=pickable,
            )
            if self.return_smi and self.method != "canon":  # canon always returns smis
                print("NOW LAST STEP")
                vals = self.pmap(self.n_jobs, self.progress, vals, self.smi_out_func)
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
            isomeric: output isomeric smiles. Defaults to True.
            kwargs: aditional keyword arguments to pass to the standardizer.

        Returns:
            standardized smiles string
        """
        mol = self._output_mol(stdin)
        try:
            standard_mol = papyrus_std.standardize(mol, **kwargs)
        except RuntimeError:
            logger.exception("Error standardizing molecule: ", stdin)
            standard_mol = None
        return standard_mol

    def chemblStandardizer(
        self, stdin: Union[str, Chem.Mol], neutralize: bool = True, **kwargs
    ) -> str:
        """Uses the ChEMBL standardizer to standardize a SMILES string. Accepts extra
        keyword arguments that will be passed to the standardizer

        Args:
            stdin: standard input; single SMILES strings or single rdkit.Chem.Mol object
                depending on the value of self._from_smi.
            isomeric: output isomeric smiles. Defaults to True.
            neutralize: configure `get_parent_mol` to neutralize the molecule. Defaults
                to True.
            kwargs: keyword arguments to pass to the `get_parent_mol` and the
                `standardize_mol` functions.

        Returns:
            standardized smiles string
        """
        mol = self._output_mol(stdin)
        if mol is None:
            return None
        try:
            parent_mol = chembl_std.get_parent_mol(
                mol,
                neutralize=neutralize,
                **{
                    "check_exclusion": kwargs.get("check_exclusion", True),
                    "verbose": kwargs.get("verbose", False),
                },
            )[0]
            standard_mol = chembl_std.standardize_mol(
                parent_mol,
                sanitize=True,
                **{"check_exclusion": kwargs.get("check_exclusion", True)},
            )
        except TypeError as e:
            logger.exception("Error standardizing molecule: ", stdin)
            logger.exception(e)
            standard_mol = None
        return standard_mol

    def molvsStandardizer(
        self,
        stdin: Union[str, Chem.Mol],
        **kwargs,
    ) -> str:
        """Uses molvs to standardize a SMILES string. By default, this standardization
        pipeline applies the functions `canonicalize_tautomer` and `standardize`
        implemented in the package.

        For more information, see the docs: https://molvs.readthedocs.io/en/latest/

        Args:
            stdin: standard input; single SMILES strings or single rdkit.Chem.Mol object
                depending on the value of self._from_smi.
            isomeric: output isomeric smiles. Defaults to True.
            kwargs: aditional `molvs.Standardizer` object.

        Returns:
            standardized smiles string
        """
        mol = self._output_mol(stdin)
        molvs_std = Standardizer(**kwargs)
        try:
            tautomer_mol = molvs_std.canonicalize_tautomer(mol)
            standard_mol = molvs_std.standardize(tautomer_mol)
        except RuntimeError:
            logger.exception("Error standardizing molecule: ", stdin)
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
        rdkit_loglevel: str = "warning",
        from_smi: bool = False,
    ) -> None:
        """Initialize the InchiHandling class.

        Args:
            convert_to: what to convert the smiles to. Can be "inchi", "inchikey" or
                "connectivity".
            n_jobs: Number of jobs for processing in parallel. Defaults to 5.
            progress: whether to show the progress bar. Defaults to False.
            rdkit_loglevel: one of `debug, info, warning, error, critical`. Defaults to
                "warning".
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
        self.rdkit_loglevel = rdkit_loglevel
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
        with rdkit_log_controller(self.rdkit_loglevel):
            mols = self.pmap(self.n_jobs, False, stdin, self._output_mol)
            vals = self.pmap(self.n_jobs, self.progress, mols, self.converter)
        return vals
