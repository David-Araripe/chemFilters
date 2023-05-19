# -*- coding: utf-8 -*-
"""Utility functions to be used in different modules of the chemFilters package."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Pool
from typing import List

import pandas as pd
import pkg_resources
from smallworld_api import SmallWorld


class SmilesSearcher:
    """Class to search for smiles in the Small World API."""

    def __init__(
        self,
        standardize_smiles: bool = True,
        n_threads: int = 5,
        n_jobs: int = 5,
    ) -> None:
        """Initialize SmilesSearcher object. Standardize_smiles is recommended as it
        will use chembl structure pipeline to standardize the smiles before searching
        and apply the same pipeline to the results, allowing to identify the order of
        the query.

        Args:
            standardize_smiles: will toggle the application of the ChEMBL structure
                pipeline to the ['qrySmiles'] returned field and the input smiles for
                the query. Allows identifying the smiles back, after multithreaded
                queries. Defaults to True.
            n_threads: number of threads for the API queries. Defaults to 5.
            n_jobs: number of jobs for the molecular standardization. Defaults to 5.
        """
        self.sw = SmallWorld()
        self.search_params = None
        self.standardize_smiles = standardize_smiles
        self.n_threads = n_threads
        self.n_jobs = n_jobs

    @property
    def smallWorldParameters(self):
        """pd.Dataframe with the parameters for the Small World API.
        Source: https://wiki.docking.org/index.php/How_to_use_SmallWorld_API
        """
        csv_path = pkg_resources.resource_filename(
            "chemFilters", "data/smallworld_parameters.csv"
        )
        return pd.read_csv(csv_path)

    @property
    def _swDBChoices(self) -> List:
        """List of available databases in Small World."""
        return self.sw.db_choices

    @property
    def _myPreferredParams(self) -> dict:
        return {"dist": 5, "db": self.sw.REAL_dataset}

    def setSearchParams(self, params: dict):
        """Set the parameters for the Small World API.
        Source: https://wiki.docking.org/index.php/How_to_use_SmallWorld_API
        """
        self.search_params = params

    def multiThreadSearch(self, smiles_list: list, **search_params) -> pd.DataFrame:
        """Search for a list of SMILES in Small World using multiple threads.

        Args:
            smiles_list: list of SMILES to search in Small World.

        Returns:
            final_df: pd.DataFrame with the results of the search.
        """
        # TODO: implement the smiles standardization / inchi key handling
        # if self.standardize_smiles:
        #     with Pool(self.n_jobs) as pool:
        #         smiles_list = pool.map(standardizer.standardize_smiles, smiles_list)

        with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
            futures = [
                executor.submit(self.sw.search, smi, **search_params)
                for smi in smiles_list
            ]
            results = []
            for future in as_completed(futures):
                results.append(future.result())
        final_df = pd.concat(results, ignore_index=True)
        return final_df
