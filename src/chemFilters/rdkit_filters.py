# -*- coding: utf-8 -*-
"""A module for filtering molecules using RDKit FilterCatalogs."""

import warnings
from functools import partial
from multiprocessing import Pool
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from .utils import get_catalog_match, smi_from_mol


class RdkitFilters:
    def __init__(self, filter_type="ALL", n_jobs=1, from_smi: bool = False) -> None:
        """Initiaze RdkitFilters object.

        Args:
            filter_type: type of filter from RDKit FilterCatalogs. Defaults to "ALL".
            n_jobs: number of jobs if wanted to run things in parallel. Defaults to 1.
            from_smi = if True, will do the conversion from SMILES to RDKit Mol object.
        """
        self._from_smi = from_smi
        self.filter_type = filter_type
        self.filter = self._get_filter()
        self.n_jobs = n_jobs

    @property
    def available_filters(self):
        """List of available filters from RDKit FilterCatalogs."""
        allFilt = FilterCatalogParams.FilterCatalogs.ALL
        return [m for m in dir(allFilt) if any([m.isupper(), m.startswith("CHEMBL_")])]

    @property
    def _filter_params(self):
        """The updated RDKit FilterCatalogParams object used for filtering."""
        catalog = FilterCatalogParams()
        catalog.AddCatalog(self.filter)
        return FilterCatalog(catalog)

    def _get_filter(self):
        """Get the filter from RDKit FilterCatalogs."""
        if self.filter_type not in self.available_filters:
            raise ValueError(f"Filter type {self.filter_type} not available.")
        return getattr(FilterCatalogParams.FilterCatalogs, self.filter_type)

    def filter_mols(
        self, mols: List[Union[Chem.Mol, str]]
    ) -> Tuple[List[List[str]], List[List[str]], List[List[Chem.Mol]]]:
        """Filter molecules using RDKit FilterCatalogs.

        Args:
            mols: list of RDKit Mol objects or SMILES strings if self._from_smi is True.

        Returns:
            filter_names: list of filter names that were matched.
            descriptions: list of filter descriptions that were matched.
            substructs: list of substructures that were matched.
        """
        if "__iter__" not in dir(mols):
            mols = [mols]
        with Pool(self.n_jobs) as p:
            filter_names, descriptions, substructs = zip(
                *p.map(
                    partial(
                        get_catalog_match,
                        catalog=self._filter_params,
                        from_smi=self._from_smi,
                    ),
                    mols,
                )
            )
        return filter_names, descriptions, substructs

    def get_flagging_df(self, mols: List[Union[Chem.Mol, str]]) -> pd.DataFrame:
        """Flag molecules using the defined RDKit FilterCatalogs and return a dataframe
        with all the detedcted filters as columns and the molecules as rows. Items
        within the dataframe will be the description of the molecular filter that was
        caught. Will also save the filter names, descriptions, and substructures as
        attributes.

        Args:
            mols: list of RDKit Mol objects or SMILES strings if self._from_smi is True.

        Returns:
            pd.DataFrame: dataframe with columns as filter types and rows as molecules.
        """
        if self.filter_type != "ALL":
            warnings.warn(
                f"Filter type {self.filter_type} is not 'ALL'. "
                "Some filters may not be applied."
            )
        filter_names, descriptions, substructs = self.filter_mols(mols)
        nope = ["PAINS", "CHEMBL", "BRENK", "ALL"]  # Are collections of filters
        columns = [c
            for c in self.available_filters
            if c not in nope and not c.startswith("CHEMBL_")
        ]
        df = pd.DataFrame(
            list(zip(filter_names, descriptions)),
            columns=["Names", "Descriptions"],
        )
        error_idx = []
        final_df = pd.DataFrame(columns=columns)
        for idx, row in df.iterrows():
            try:
                label_dict = dict(zip(row["Names"], row["Descriptions"]))
            except TypeError:
                label_dict = {}
                error_idx.append(idx)
            final_df = pd.concat(
                [final_df, pd.DataFrame(label_dict, index=[0])], ignore_index=True
            )
        final_df = final_df.applymap(lambda x: [] if pd.isnull(x) else [x])
        for col in final_df.columns:
            final_df[col] = final_df[col].apply(lambda x: ";".join(x))
        if self._from_smi:
            final_df.insert(0, "SMILES", mols)
        else:
            with Pool(self.n_jobs) as p:
                final_df.insert(0, "SMILES", p.map(smi_from_mol, mols))
        if error_idx != []:
            final_df.loc[error_idx, "SMILES"] = np.nan
            warnings.warn(
                f"Failed to get filter names and descriptions for {len(error_idx)} "
                f"molecules in indexes: {error_idx}."
            )
        self.filter_names = filter_names
        self.descriptions = descriptions
        self.substructs = substructs
        return final_df.replace({"": np.nan})
