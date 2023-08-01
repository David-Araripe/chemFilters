# -*- coding: utf-8 -*-
"""A module for filtering molecules using RDKit FilterCatalogs."""

import logging
from functools import partial
from multiprocessing import Pool
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from ..chem.interface import MoleculeHandler
from .utils import get_catalog_match

STANDARD_RD_COLS = [
    "NIH",
    "PAINS_A",
    "PAINS_B",
    "PAINS_C",
    "PAINS_any",
    "ZINC",
    "Brenk",
    "ChEMBL23_Dundee",
    "ChEMBL23_BMS",
    "ChEMBL23_MLSMR",
    "ChEMBL23_Inpharmatica",
    "ChEMBL23_SureChEMBL",
    "ChEMBL23_LINT",
    "ChEMBL23_Glaxo",
    "ChEMBL23_any",
]
FILTER_COLLECTIONS = ["PAINS", "CHEMBL", "BRENK", "ALL"]  # Are collections of filters


class RdkitFilters(MoleculeHandler):
    def __init__(self, filter_type="ALL", n_jobs=1, from_smi: bool = False) -> None:
        """Initiaze RdkitFilters object.

        Args:
            filter_type: type of filter from RDKit FilterCatalogs. Defaults to "ALL".
            n_jobs: number of jobs if wanted to run things in parallel. Defaults to 1.
            from_smi = if True, will do the conversion from SMILES to RDKit Mol object.
        """
        self.filter_type = filter_type
        self.filter = self._get_filter()
        self.n_jobs = n_jobs
        super().__init__(from_smi=from_smi)

    @property
    def available_filters(self):
        """List of available filters from RDKit FilterCatalogs."""
        allFilt = FilterCatalogParams.FilterCatalogs.ALL
        return [m for m in dir(allFilt) if any([m.isupper(), m.startswith("CHEMBL_")])]

    def _get_filter(self):
        """Get the filter from RDKit FilterCatalogs."""
        if self.filter_type not in self.available_filters:
            raise ValueError(f"Filter type {self.filter_type} not available.")
        _filter = getattr(FilterCatalogParams.FilterCatalogs, self.filter_type)
        catalog = FilterCatalogParams()
        catalog.AddCatalog(_filter)
        return FilterCatalog(catalog)

    def filter_mols(
        self,
        stdin: List[Union[Chem.Mol, str]],
        match_type: str = "string",
    ) -> Tuple[List[List[str]], List[List[str]], List[List[Chem.Mol]]]:
        """Filter molecules using RDKit FilterCatalogs.

        Args:
            stdin: list of RDKit Mol objects of SMILES strings if self._from_smi is True
            match_type: values within the flagging dataframe. If `bool`, will spare
                retrieving substructures and descriptions. If `string`, will have the
                description of the filter that was matched. Defaults to `string`.

        Returns:
            filter_names: list of filter names that were matched.
            descriptions: list of filter descriptions that were matched.
            substructs: list of substructures that were matched.
        """
        with Pool(self.n_jobs) as p:
            mols = p.map(self._output_mol, stdin)
            result = p.map(
                partial(
                    get_catalog_match,
                    catalog=self.filter,
                    match_type=match_type,
                ),
                mols,
            )
        filter_names, descriptions, substructs = zip(*result)
        return filter_names, descriptions, substructs

    def get_flagging_df(
        self,
        stdin: List[Union[Chem.Mol, str]],
        match_type: str = "string",
        save_matches: bool = False,
    ) -> pd.DataFrame:
        """Flag molecules using the defined RDKit FilterCatalogs and return a dataframe
        with all the detedcted filters as columns and the molecules as rows. Items
        within the dataframe will be the description of the molecular filter that was
        caught. Will also save the filter names, descriptions, and substructures as
        attributes.

        Args:
            stdin: list of RDKit Mol objects or SMILES strings if self._from_smi is True
            match_type: values within the flagging dataframe. If `bool`, will spare
                retrieving substructures and descriptions. If `string`, will have the
                description of the filter that was matched. Defaults to `string`.
            save_matches: if True, will save the filter names, descriptions, and
                substructures as attributes. Defaults to False.

        Returns:
            pd.DataFrame: dataframe with columns as filter types and rows as molecules.
        """
        if match_type.lower() not in ["string", "bool"]:
            raise ValueError("match_type must be either 'string' or 'bool'.")

        def flatten_labels(labels):
            """Flatten the labels and filter out any None types"""
            return ";".join(filter(None, labels))

        vectorized_flatten = np.vectorize(flatten_labels)
        filter_names, descriptions, substructs = self.filter_mols(
            stdin, match_type=match_type
        )
        if match_type.lower() == "string":
            val_dicts = [
                dict(zip(names, descs))
                if all([names is not None, descs is not None])
                else {}
                for names, descs in zip(filter_names, descriptions)
            ]
            final_df = pd.DataFrame(val_dicts)
            final_df = final_df.applymap(lambda x: [] if pd.isnull(x) else [x])
            final_df = (
                final_df.apply(vectorized_flatten)
                .replace({"": np.nan})
                .reindex(columns=STANDARD_RD_COLS)
            )
        elif match_type.lower() == "bool":
            final_df = pd.DataFrame(columns=STANDARD_RD_COLS, index=range(len(stdin)))
            for col in STANDARD_RD_COLS:
                names_series = pd.Series(
                    [names if names is not None else [] for names in filter_names]
                )
                final_df[col] = names_series.apply(lambda x, col: col in x, col=col)
        # Add smiles to the final dataframe
        with Pool(self.n_jobs) as p:
            smiles = p.map(self._output_smi, stdin)
        final_df.insert(0, "SMILES", smiles)
        # if there are any errors, set the smiles to NaN
        error_idx = [i for i, x in enumerate(filter_names) if x is None]
        if error_idx:
            final_df.loc[error_idx, "SMILES"] = np.nan
            logging.warn(
                f"Failed to get filter names and descriptions for {len(error_idx)} "
                f"molecules in indexes: {error_idx}. SMILES will be set to NaN."
            )
        if save_matches:
            self.filter_names = filter_names
            self.descriptions = descriptions
            self.substructs = substructs
        return final_df
