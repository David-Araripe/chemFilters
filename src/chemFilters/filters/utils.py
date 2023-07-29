# -*- coding: utf-8 -*-
"""Utility functions to be used in different modules of the chemFilters package."""

import numpy as np
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog

from .chem.interface import mol_from_smi


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


def geometric_mean(arr, axis=0):
    """Calculate the geometric mean of an array. Adapted from SciPy:
    https://github.com/scipy/blob/v1.10.1/scipy/stats/_stats.py.py#L199=L269"""
    with np.errstate(divide="ignore"):
        log_a = np.log(arr)
    return np.exp(np.average(log_a, axis=axis))
