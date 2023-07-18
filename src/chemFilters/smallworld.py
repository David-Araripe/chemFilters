# -*- coding: utf-8 -*-
"""Utility functions to be used in different modules of the chemFilters package."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from importlib.util import find_spec
from typing import List

import pandas as pd
import pkg_resources
from smallworld_api import SmallWorld

from chemFilters.chem.standardizers import InchiHandling, SmilesStandardizer


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
            "chemFilters", "resources/sw_parameters.csv"
        )
        return pd.read_csv(csv_path)

    @property
    def _swDBChoices(self) -> List:
        """List of available databases in Small World."""
        return self.sw.db_choices

    @property
    def _myPreferredParams(self) -> dict:
        return {
            "dist": 3,  # I want molecules that are fairly similar to the query
            "db": self.sw.REAL_dataset,
            "length": 20,  # I want to get 20 results
        }

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
        if self.standardize_smiles:
            print("Standardizing SMILES...")
            standardizer = SmilesStandardizer(n_jobs=self.n_jobs)
            smiles_list = standardizer(smiles_list)
            print("Calculating InchiKeys...")
            inchi_calculator = InchiHandling("inchikey", n_jobs=self.n_jobs)
            inchikeys = inchi_calculator(smiles_list)
            idx_mapping = {k: idx for idx, k in enumerate(inchikeys)}

        with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
            futures = [
                executor.submit(self.sw.search, smi, **search_params)
                for smi in smiles_list
            ]
            results = []
            for future in as_completed(futures):
                results.append(future.result())

        final_df = pd.concat(results, ignore_index=True)
        if self.standardize_smiles:
            final_df = final_df.assign(
                qrySmiles_std=lambda x: standardizer(x["qrySmiles"]),
                qryInchiKey_std=lambda x: inchi_calculator(x["qrySmiles_std"]),
                submission_idx=lambda x: x["qryInchiKey_std"].map(idx_mapping),
            ).sort_values("submission_idx")
        return final_df

    def lazyMultiThreadSearch(
        self, smiles_list: list, output_file: str, **search_params
    ) -> pd.DataFrame:
        """Search for a list of SMILES in Small World using multiple threads in a lazy
        way. This means that the results will be written to disk as they are obtained to
        prevent the rest of the results being lost from a crash.

        Args:
            smiles_list: list of SMILES to search in Small World.
            output_file: path to the output file where the results will be written.
                Note: this method requires the optional dependency pytables, as .csv
                doesn't directly support lazy writing.

        Returns:
            final_df: pd.DataFrame with the results of the search.
        """
        if not find_spec("tables"):
            raise ImportError(
                "Lazy writing requires the optional dependency pytables."
                "To install it, run: \npython -m pip install tables"
            )
        if self.standardize_smiles:
            print("Standardizing SMILES...")
            standardizer = SmilesStandardizer(n_jobs=self.n_jobs)
            smiles_list = standardizer(smiles_list)
            print("Calculating InchiKeys...")
            inchi_calculator = InchiHandling("inchikey", n_jobs=self.n_jobs)
            inchikeys = inchi_calculator(smiles_list)
            idx_mapping = {k: idx for idx, k in enumerate(inchikeys)}

        with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
            futures = [
                executor.submit(self.sw.search, smi, **search_params)
                for smi in smiles_list
            ]
            results = []
            for future in as_completed(futures):
                result = future.result()
                # append results to the output file
                result.to_hdf(output_file, "results", append=True)

        final_df = pd.concat(results, ignore_index=True)
        if self.standardize_smiles:
            final_df = final_df.assign(
                qrySmiles_std=lambda x: standardizer(x["qrySmiles"]),
                qryInchiKey_std=lambda x: inchi_calculator(x["qrySmiles_std"]),
                submission_idx=lambda x: x["qryInchiKey_std"].map(idx_mapping),
            ).sort_values("submission_idx")
        return final_df
