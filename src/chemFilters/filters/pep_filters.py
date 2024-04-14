# -*- coding: utf-8 -*-
from itertools import product
from multiprocessing import Pool
from typing import List, Union

import numpy as np
import pandas as pd
from pepsift import PepSift, SiftLevel
from rdkit import Chem

from ..chem.interface import MoleculeHandler
from ..logger import logger

STANDARD_PEP_COLS = [
    "NaturalLAminoAcids",
    "NaturalLDAminoAcids",
    "NaturalAminoAcidDerivatives",
    "NonNaturalAminoAcidDerivatives",
    "AllAmineAndAcid",
]


class PeptideFilters(MoleculeHandler):
    """Wrapper class for PepSift, a tool for identifying peptides and their derivatives
    from small molecule datasets. For the original repo, see:
    https://github.com/OlivierBeq/PepSift/tree/master

    Arguments:
        filter_type: filter type to initialize a PepSift object. See available
            filters on `self.available filters`. Defaults to "all".
        from_smi: treats standard inputs (stdin) as smiles. Defaults to False.
        n_jobs: number of jobs to run in parallel. Defaults to 1.
    """

    def __init__(
        self,
        filter_type: Union[str, int] = "all",
        from_smi: bool = False,
        n_jobs: int = 1,
    ) -> None:
        """Initialize the PeptideFilters class.

        Args:
            filter_type: filter type to initialize a PepSift object. See available
                filters on `self.available filters`. Defaults to "all".
            from_smi: treats standard inputs (stdin) as smiles. Defaults to False.
            n_jobs: number of jobs to run in parallel. Defaults to 1.
        """
        self.filter_type = filter_type
        self.filter = self._get_filter(_type=filter_type)
        self._n_jobs = n_jobs
        super().__init__(from_smi=from_smi)

    @property
    def available_filters(self):
        """List of available filters on pepsift."""
        return dict(SiftLevel.__members__.items())

    def _get_filter(self, _type: Union[str, int] = None):
        """Helper function to get a PepSift object."""
        if isinstance(_type, SiftLevel):
            sift_object = _type
        elif _type in ["NaturalLAminoAcids", 1, "naturall"]:
            sift_object = SiftLevel(1)
        elif _type in ["NaturalLDAminoAcids", 2, "naturald"]:
            sift_object = SiftLevel(2)
        elif _type in ["NaturalAminoAcidDerivatives", 3, "naturalderivative"]:
            sift_object = SiftLevel(3)
        elif _type in ["NonNaturalAminoAcidDerivatives", 4, "nonnaturalderivative"]:
            sift_object = SiftLevel(4)
        elif _type in ["AllAmineAndAcid", 5, "all"]:
            sift_object = SiftLevel(5)
        else:
            raise ValueError(
                f"Filter type {self.filter_type} not available. "
                f"Try one of:\n{self.available_filters}"
            )
        return PepSift(sift_object)

    def _peptide_filter_func(
        self, stdin: Union[str, Chem.Mol], sift_obj: PepSift = None
    ) -> bool:
        """Function to be used in parallelization. Will return True if the molecule
        is a peptide according to the filtering level used to initialize sift_object,
        False otherwise.

        Args:
            stdin: standard input; single SMILES strings or single rdkit.Chem.Mol object
                depending on the value of self._from_smi.
            sift_obj: pepsift object to filter the molecules if `None`, will use
                `self.filter`. Defaults to None.

        Returns:
            bool: True if the molecule is a peptide, False otherwise.
        """
        mol = self._output_mol(stdin)
        if mol is None:
            return None
        if sift_obj is None:
            assert isinstance(self.filter, PepSift), "self.filter must be a PepSift obj"
            sift_obj = self.filter
        try:
            return sift_obj.is_peptide(mol)
        except Exception as e:
            logger.error(
                "PepSift.is_peptide() error - returning np.nan for sift level "
                f"{sift_obj.level}. Message:\n {e} for stdin {stdin}",
            )
            return None

    def filter_mols(self, stdin: List[Chem.Mol]):
        """Filter molecules using the designated pepsift filter. If `sift_level=None`
        as default, will load it from `self.filter`.

        Args:
            stdin: standard input; a list of SMILES strings or rdkit.Chem.Mol objects
                depending on the value of self._from_smi.

        Returns:
            List[bool]: a list of booleans indicating whether the molecule is a peptide
                according to the initialized filter level on `self.filter`.
        """
        with Pool(self._n_jobs) as p:
            bool_mask = p.map(self._peptide_filter_func, stdin)
        return bool_mask

    def get_flagging_df(self, stdin: List[Union[str, Chem.Mol]]):
        """Will flag the molecules according to all filter types avialable in pepsift.

        Args:
            stdin: standard input; a list of SMILES strings or rdkit.Chem.Mol objects
                depending on the value of self._from_smi.

        Returns:
            pd.DataFrame: dataframe with the flags for each filter type.
        """
        all_filters = [self._get_filter(idx) for idx in range(1, 6)]
        all_params = list(product(stdin, all_filters))

        with Pool(self._n_jobs) as p:
            results = p.starmap(self._peptide_filter_func, all_params)
            smiles = p.map(self._output_smi, stdin)

        df = pd.DataFrame(
            np.array(results).reshape(-1, len(all_filters)),
            columns=self.available_filters.keys(),
        )
        df.insert(0, "SMILES", smiles)
        return df.replace({None: np.nan})
