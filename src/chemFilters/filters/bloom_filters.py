# -*- coding: utf-8 -*-
"""Wrapper class for molbloom. Original repo: https://github.com/whitead/molbloom"""
from itertools import product
from multiprocessing import Pool
from pathlib import Path
from typing import List, Union

import numpy as np
import pandas as pd
from molbloom import _DEFAULT_PATH, _load_filter, buy, catalogs
from rdkit import Chem

from .chem.interface import MoleculeHandler
from .chem.standardizers import ChemStandardizer


class MolbloomFilters(MoleculeHandler):
    """Wrapper class for molbloom. Requires molbloom to be installed.

    Arguments:
        from_smi: treats standard inputs (stdin) as smiles. Defaults to False.
        standardize: whether to standardize `stdin` or not. Defaults to False.
        std_method: SMILES/mol standardization method. Available: `canon`, `chembl`,
            `papyrus`. Defaults to "chembl".
        n_jobs: number of jobs to run in parallel. Defaults to 1."""

    def __init__(
        self,
        from_smi: bool = False,
        standardize: bool = False,
        std_method: str = "chembl",
        n_jobs=1,
        **kwargs,
    ):
        """Initalize the MolbloomFilters class.

        Args:
            from_smi: treats standard inputs (stdin) as smiles. Defaults to False.
            standardize: whether to standardize `stdin` or not. Defaults to False.
            std_method: SMILES/mol standardization method. Available: `canon`, `chembl`,
                `papyrus`. Defaults to "chembl".
            n_jobs: number of jobs to run in parallel. Defaults to 1.
        """
        self._catalogs = catalogs()
        self._ensure_filters_downloaded()
        self._standardize = standardize
        self._smiles_standardizer = self._get_standardizer(
            std_method, from_smi, **kwargs
        )
        self._std_method = std_method.lower()
        self._kwargs = kwargs
        self._n_jobs = n_jobs
        super().__init__(from_smi)

    def _get_standardizer(self, std_method: str, from_smi: bool, **kwargs):
        return ChemStandardizer(method=std_method, from_smi=from_smi, **kwargs)

    def get_catalogs(self):
        """Avilable catalogs in molbloom."""
        return self._catalogs

    def _ensure_filters_downloaded(self):
        """Ensures that the molbloom filters are downloaded."""
        for catalog in self._catalogs:
            filter_path = Path(_DEFAULT_PATH) / f"{catalog}.bloom"
            if not filter_path.exists():
                print(f"Starting {catalog} download...")
                _load_filter(catalog)

    def buy_smi(self, smi: str, catalog: str = "zinc-instock"):
        """Wrapper of molbloom.buy. Returns True if the SMILES is probably in the
        catalog, False if it is definitely not."""
        try:
            # canonicalization off as it's handled by the smiles standardizer
            return buy(smi, catalog=catalog, canonicalize=False)
        except Exception as e:
            print(f"An error occurred: {str(e)} on buying {smi} from {catalog}")
            return None

    def get_flagging_df(self, stdin: List[Union[str, Chem.Mol]]):
        """Returns a dataframe with the flagging results for each catalog. Flags will
        be the resutls from molbloom.buy, where True means the SMILES is probably in
        the catalog, False means it is definitely not. For more information, see the
        original repo: https://github.com/whitead/molbloom

        Args:
            stdin: standard input; a list of SMILES strings or rdkit.Chem.Mol objects
                depending on the value of self._from_smi.

        Returns:
            pd.DataFrame: dataframe with the flagging results for each catalog.
        """
        if self._standardize:
            smiles = self._smiles_standardizer(stdin)
        else:
            with Pool(self._n_jobs) as p:
                smiles = p.map(self._output_smi, stdin)

        all_params = list(product(smiles, self._catalogs))
        with Pool(self._n_jobs) as p:
            results = p.starmap(self.buy_smi, all_params)

        num_catalogs = len(self._catalogs)
        df = pd.DataFrame(
            np.array(results).reshape((-1, num_catalogs)), columns=self._catalogs.keys()
        )
        df.insert(0, "SMILES", smiles)
        return df.replace({None: np.nan})
