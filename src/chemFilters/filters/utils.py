# -*- coding: utf-8 -*-
"""Utility functions to be used in different modules of the chemFilters package."""

from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog


def get_catalog_match(
    mol: Chem.Mol, catalog: FilterCatalog, match_type: str = "string"
):
    """Get the filter names, descriptions, and respective substructures from
    rdkit.Chem.FilterCatalog.FilterCatalogParams.GetMatches. Used to get results in
    parallel using multiprocessing.Pool. The match_type parameter is used to save time
    on `RdkitFilters.get_flagging_df` by not having to return the substructures if not
    desred by the user.

    Args:
        mol: rdkitm.Chem.Mol object to be filtered.
        catalog: FilterCatalog object from rdkit.Chem.FilterCatalog.
        match_type: if `string` will return the descriptions and the substructures.
            If `bool` will return only the filter_names. Defaults to "string".

    Returns:
        list[str], list[str], list[rdkit.Chem.Mol]: filter_names, descriptions,
            substructures.
    """
    if mol is None:
        return None, None, None
    matches = catalog.GetMatches(mol)
    filter_names = [m.GetProp("FilterSet") for m in matches]
    if match_type == "bool":
        return filter_names, None, None
    elif match_type == "string":
        descriptions = [m.GetDescription() for m in matches]
        substructs = [
            m.GetFilterMatches(mol)[0].filterMatch.GetPattern() for m in matches
        ]
        return filter_names, descriptions, substructs
