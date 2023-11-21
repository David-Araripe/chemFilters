# -*- coding: utf-8 -*-
"""Utility functions to be used in different modules of the chemFilters package."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from importlib.util import find_spec
from pathlib import Path
from typing import List

import pandas as pd
import pkg_resources
from smallworld_api import SmallWorld

from .chem.standardizers import ChemStandardizer, InchiHandling


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
        self.inchi_calculator = InchiHandling(
            "inchikey", n_jobs=self.n_jobs, from_smi=True, verbose=False
        )
        self.standardizer = ChemStandardizer(
            n_jobs=self.n_jobs, from_smi=True, verbose=False
        )

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
            "dist": 5,  # I want molecules that are fairly similar to the query
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
            smiles_list = self.standardizer(smiles_list)
            inchikeys = self.inchi_calculator(smiles_list)
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
                qrySmiles_std=lambda x: self.standardizer(x["qrySmiles"]),
                qryInchiKey_std=lambda x: self.inchi_calculator(x["qrySmiles_std"]),
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

        Returns:
            final_df: pd.DataFrame with the results of the search.
        """
        if self.standardize_smiles:
            smiles_list = self.standardizer(smiles_list)
            inchikeys = self.inchi_calculator(smiles_list)
            idx_mapping = {k: idx for idx, k in enumerate(inchikeys)}

        with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
            futures = [
                executor.submit(self.sw.search, smi, **search_params)
                for smi in smiles_list
            ]
            results = []
            for future in as_completed(futures):
                result = future.result()
                if len(result) > 0:
                    result = result.assign(
                        atomMap=lambda x: x["atomMap"].astype(str),
                        atomScore=lambda x: x["atomScore"].astype(str),
                        anonIdx=lambda x: x["anonIdx"].astype(str),
                    )
                    if not Path(output_file).exists():  # In this case, keep the header
                        result.to_csv(output_file, mode="a", index=False)
                    else:
                        result.to_csv(output_file, mode="a", index=False, header=False)
                    results.append(result)
                else:
                    continue

        if results != []:
            final_df = pd.concat(results, ignore_index=True)
            if self.standardize_smiles:
                final_df = final_df.assign(
                    qrySmiles_std=lambda x: self.standardizer(x["qrySmiles"]),
                    qryInchiKey_std=lambda x: self.inchi_calculator(x["qrySmiles_std"]),
                    submission_idx=lambda x: x["qryInchiKey_std"].map(idx_mapping),
                ).sort_values("submission_idx")
            return final_df
        else:
            return
