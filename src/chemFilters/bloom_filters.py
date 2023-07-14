# -*- coding: utf-8 -*-
"""Wrapper class for molbloom. Original repo: https://github.com/whitead/molbloom"""
from itertools import product
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from molbloom import _DEFAULT_PATH, _load_filter, buy, catalogs

from .chem.standardizers import SmilesStandardizer
from .utils import smi_from_mol


class MolbloomFilters:
    def __init__(
        self, from_smi=True, standardize=True, std_method="chembl", n_jobs=1, **kwargs
    ):
        self._catalogs = catalogs()
        self._ensure_filters_downloaded()
        self._from_smi = from_smi
        self._standardize = standardize
        self._std_method = std_method
        self._kwargs = kwargs
        self._n_jobs = n_jobs

    @property
    def _smiles_standardizer(self):
        return SmilesStandardizer(method=self._std_method, **self._kwargs)

    def _process_mols(self, mols):
        if self._from_smi is False:
            with Pool(self._n_jobs) as p:
                return p.map(smi_from_mol, mols)
        else:
            return mols

    def get_catalogs(self):
        return self._catalogs

    def _ensure_filters_downloaded(self):
        for catalog in self._catalogs:
            filter_path = Path(_DEFAULT_PATH) / f"{catalog}.bloom"
            if not filter_path.exists():
                print(f"Starting {catalog} download...")
                _load_filter(catalog)

    def buy_mol(self, smi: str, catalog: str = "zinc-instock"):
        try:
            # canonicalization off as it's handled by the smiles standardizer
            return buy(smi, catalog=catalog, canonicalize=False)
        except Exception as e:
            print(f"An error occurred: {str(e)}")
            return None

    def get_flagging_df(self, smiles):
        if self._standardize:
            smiles = self._smiles_standardizer(smiles)

        all_params = list(product(smiles, self._catalogs))
        with Pool(self._n_jobs) as p:
            results = p.starmap(self.buy_mol, all_params)

        num_catalogs = len(self._catalogs)
        df = pd.DataFrame(
            np.array(results).reshape((-1, num_catalogs)), columns=self._catalogs.keys()
        )
        df.insert(0, "SMILES", smiles)
        return df.replace({None: np.nan})
