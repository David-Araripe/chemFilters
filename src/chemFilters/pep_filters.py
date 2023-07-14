# -*- coding: utf-8 -*-
from functools import partial
from multiprocessing import Pool
from typing import List, Union

import numpy as np
import pandas as pd
from pepsift import PepSift, SiftLevel
from rdkit import Chem

from .utils import smi_to_mol


def _peptide_filter_func(
    mol: Chem.Mol, sift_obj: PepSift = None, from_smi: bool = False
):
    if from_smi:
        mol = smi_to_mol(mol)
    if mol is None:
        return None
    else:
        return sift_obj.is_peptide(mol)


class PeptideFilters:
    def __init__(
        self,
        filterType: Union[str, int] = "all",
        njobs: int = 1,
        from_smi: bool = False,
    ) -> None:
        """Initialize PeptideFilters object.

        Args:
            filterType: PepSift's filters to be applied. Defaults to "all"
            njobs: number of jobs if wanted to run things in parallel. Defaults to 1.
        """
        self.from_smi = from_smi
        self.filterType = filterType
        self.filter = self._get_filter()
        self.njobs = njobs

    @property
    def available_filters(self):
        """List of available filters on pepsift."""
        return dict(SiftLevel.__members__.items())

    def _get_filter(self):
        """Get the filter from PepSift."""
        filtname = self.filterType
        if filtname in ["NaturalLAminoAcids", 1, "naturall"]:
            _filter = SiftLevel(1)
        elif filtname in ["NaturalLDAminoAcids", 2, "naturald"]:
            _filter = SiftLevel(2)
        elif filtname in ["NaturalAminoAcidDerivatives", 3, "naturalderivative"]:
            _filter = SiftLevel(3)
        elif filtname in ["NonNaturalAminoAcidDerivatives", 4, "nonnaturalderivative"]:
            _filter = SiftLevel(4)
        elif filtname in ["AllAmineAndAcid", 5, "all"]:
            _filter = SiftLevel(5)
        else:
            raise ValueError(
                f"Filter type {self.filterType} not available. "
                f"Try one of:\n{self.available_filters}"
            )
        return _filter

    def filter_mols(self, mols: List[Chem.Mol], sift_level: SiftLevel = None):
        """Filter molecules using the designated pepsift filter. If `sift_level=None`
        as default, will load it from `self.filter`.

        Args:
            mols: list of RDKit Mol objects or SMILES strings if self.from_smi is True.

        Returns:
            filter_names: list of filter names that were matched.
            descriptions: list of filter descriptions that were matched.
            substructs: list of substructures that were matched.
        """
        if sift_level is None:
            sift_level = self.filter
        # This is wrong...
        sift_obj = PepSift(sift_level)
        with Pool(self.njobs) as p:
            bool_mask = p.map(
                partial(
                    _peptide_filter_func, sift_obj=sift_obj, from_smi=self.from_smi
                ),
                mols,
            )
        return bool_mask

    def get_flagging_df(self, mols: List[Chem.Mol]):
        """Will flag the molecules according to all filter types avialable in pepsift.

        Args:
            mols: list of rdkit molecule objects.

        Returns:
            pd.DataFrame: dataframe with the flags for each filter type.
        """
        final_df = pd.DataFrame(columns=self.available_filters.keys())

        for name, level in self.available_filters.items():
            final_df[name] = self.filter_mols(mols, sift_level=level)
        return final_df.replace({None: np.nan})  # In case molecules were None
