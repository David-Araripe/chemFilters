# -*- coding: utf-8 -*-
"""Utility functions to be used in different modules of the chemFilters package."""

import warnings

from rdkit import Chem


def smiToMol(smi: str):
    """Convert a SMILES string to a rdkit.Chem.Mol object."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        warnings.warn(f"Could not convert SMILES {smi} to rdkit.Chem.Mol")
    return mol


def getCatalogMatch(mol: Chem.Mol, catalog: Chem.FilterCatalog, from_smi: bool = False):
    """Get the filter names and descriptions from
    rdkit.Chem.FilterCatalog.FilterCatalogParams.GetMatches.

    Used to get results in parallel using multiprocessing.Pool.
    """
    if from_smi:
        mol = smiToMol(mol)
    if mol is None:
        return None, None
    matches = catalog.GetMatches(mol)
    filter_names = [m.GetProp("FilterSet") for m in matches]
    descriptions = [m.GetDescription() for m in matches]
    return filter_names, descriptions


def getCatalogFilterMatch(
    mol: Chem.Mol, catalog: Chem.FilterCatalog, from_smi: bool = False
):
    """Get the substructure patterns from
    rdkit.Chem.FilterCatalog.FilterCatalogParams.GetFilterMatches.

    Used to get results in parallel using multiprocessing.Pool.
    """
    if from_smi:
        mol = smiToMol(mol)
    if mol is None:
        return None
    filter_matches = catalog.GetFilterMatches(mol)
    substructs = [m.filterMatch.GetPattern() for m in filter_matches]
    return substructs
