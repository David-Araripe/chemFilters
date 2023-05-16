# -*- coding: utf-8 -*-
from rdkit import Chem


def getCatalogMatch(mol: Chem.Mol, catalog: Chem.FilterCatalog):
    """Get the filter names and descriptions from
    rdkit.Chem.FilterCatalog.FilterCatalogParams.GetMatches.

    Used to get results in parallel using multiprocessing.Pool.
    """
    matches = catalog.GetMatches(mol)
    filter_names = [m.GetProp("FilterSet") for m in matches]
    descriptions = [m.GetDescription() for m in matches]
    return filter_names, descriptions


def getCatalogFilterMatch(mol: Chem.Mol, catalog: Chem.FilterCatalog):
    """Get the substructure patterns from
    rdkit.Chem.FilterCatalog.FilterCatalogParams.GetFilterMatches.

    Used to get results in parallel using multiprocessing.Pool.
    """
    filter_matches = catalog.GetFilterMatches(mol)
    substructs = [m.filterMatch.GetPattern() for m in filter_matches]
    return substructs
