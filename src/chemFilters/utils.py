# -*- coding: utf-8 -*-
"""Utility functions to be used in different modules of the chemFilters package."""

import logging
import warnings
from typing import Union

import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.FilterCatalog import FilterCatalog


def smi_to_mol(smi: str):
    """Convert a SMILES string to a rdkit.Chem.Mol object."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        warnings.warn(f"Could not convert SMILES {smi} to rdkit.Chem.Mol")
    return mol


def smi_from_mol(mol: Chem.Mol):
    """Convert a rdkit.Chem.Mol object to a SMILES string."""
    smi = Chem.MolToSmiles(mol)
    if smi is None:
        warnings.warn(f"Could not convert rdkit.Chem.Mol to SMILES {mol}")
    return smi


def get_catalog_match(
    mol: Chem.Mol, catalog: FilterCatalog, from_smi: bool = False
):
    """Get the filter names, descriptions, and respective substructures from
    rdkit.Chem.FilterCatalog.FilterCatalogParams.GetMatches.

    Used to get results in parallel using multiprocessing.Pool.
    """
    if from_smi:
        mol = smi_to_mol(mol)
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
