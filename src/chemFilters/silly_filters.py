# -*- coding: utf-8 -*-
"""Wrapper function for molspotter, a tool based on Pat Water's silly walks (hence the module name)"""  # noqa: E501
from itertools import product
from multiprocessing import Pool
from typing import List, Union

import numpy as np
import pandas as pd
from molspotter import SillyMolSpotter
from rdkit import Chem

from chemFilters.chem.standardizers import ChemStandardizer

from .chem.interface import MoleculeHandler


class SillyMolSpotterFilter(MoleculeHandler):
    """Wrapper class to molspotter, a tool based on Pat Water's silly walks filter.
    It helps finding unusual molecules in a dataset based the detection of unusual bits
    on a hashed ECFP fingerprint. For more information, see the original repo:
    https://github.com/OlivierBeq/molspotter

    Arguments:
        from_smi: treats standard inputs (stdin) as smiles. Defaults to False.
        standardize: whether to standardize `stdin` or not. Defaults to False.
        std_method: SMILES/mol standardization method. Available: `canon`, `chembl`,
            `papyrus`. Defaults to "chembl".
        n_jobs: number of jobs to run in parallel. Defaults to 1.
    """

    def __init__(
        self, from_smi=False, standardize=False, std_method="chembl", n_jobs=1, **kwargs
    ):
        """Initialize the SillyMolSpotterFilter class.

        Args:
            from_smi: treats standard inputs (stdin) as smiles. Defaults to False.
            standardize: whether to standardize `stdin` or not. Defaults to False.
            std_method: SMILES/mol standardization method. Available: `canon`, `chembl`,
                `papyrus`. Defaults to "chembl".
            n_jobs: number of jobs to run in parallel. Defaults to 1.
        """
        self._spotters = self._get_spotters()
        self._standardize = standardize
        self._smiles_standardizer = self._get_standardizer(
            std_method, from_smi=from_smi, n_jobs=n_jobs, **kwargs
        )
        self._std_method = std_method.lower()
        self._kwargs = kwargs
        self._n_jobs = n_jobs
        super().__init__(from_smi)

    def _get_standardizer(self, std_method, from_smi=True, n_jobs=1, **kwargs):
        return ChemStandardizer(
            method=std_method, from_smi=from_smi, n_jobs=n_jobs, **kwargs
        )

    def _get_spotters(self):
        """Available pretrained spotters from molspotter."""
        return {
            "chembl": SillyMolSpotter.from_pretrained("chembl"),
            "excape": SillyMolSpotter.from_pretrained("excape"),
            "papyrus": SillyMolSpotter.from_pretrained("papyrus"),
        }

    def score_compound(self, stdin: str, spotter_name: str = "chembl"):
        """Score a SMILES string with a pretrained spotter, indicating how silly the
        processed molecule is."""
        spotter = self._spotters[spotter_name]
        try:
            mol = self._output_mol(stdin)
            return spotter.score_mol(mol)
        except Exception as e:
            print(
                f"An error occurred: {str(e)} when scoring "
                f"{stdin} with {spotter_name}."
            )
            return None

    def get_scoring_df(self, stdin: List[Union[str, Chem.Mol]]):
        """Get a dataframe with the scoring results for each spotter.

        Args:
            stdin: standard input; a list of SMILES strings or rdkit.Chem.Mol objects
                depending on the value of self._from_smi.

        Returns:
            pd.DataFrame: a dataframe with the scoring results for each spotter.
        """
        if self._standardize:
            smiles = self._smiles_standardizer(stdin)
            if smiles == []:
                raise ValueError("No valid SMILES found in the input.")
        else:
            with Pool(self._n_jobs) as p:
                smiles = p.map(self._output_smi, stdin)

        with Pool(self._n_jobs) as p:
            all_params = list(product(smiles, self._spotters.keys()))
            results = p.starmap(self.score_compound, all_params)

        num_spotters = len(self._spotters)
        df = pd.DataFrame(
            np.array(results).reshape((-1, num_spotters)), columns=self._spotters.keys()
        )
        df.insert(0, "SMILES", smiles)
        return df.replace({None: np.nan})
