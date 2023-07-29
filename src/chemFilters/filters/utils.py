# -*- coding: utf-8 -*-
"""Utility functions to be used in different modules of the chemFilters package."""

from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog

from ..chem.interface import mol_from_smi


def get_catalog_match(mol: Chem.Mol, catalog: FilterCatalog, from_smi: bool = False):
    """Get the filter names, descriptions, and respective substructures from
    rdkit.Chem.FilterCatalog.FilterCatalogParams.GetMatches.

    Used to get results in parallel using multiprocessing.Pool.
    """
    if from_smi:
        mol = mol_from_smi(mol)
    if mol is None:
        return None, None, None
    matches = catalog.GetMatches(mol)
    filter_names = [m.GetProp("FilterSet") for m in matches]
    descriptions = [m.GetDescription() for m in matches]
    substructs = [m.GetFilterMatches(mol)[0].filterMatch.GetPattern() for m in matches]
    return filter_names, descriptions, substructs
