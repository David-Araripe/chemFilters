# -*- coding: utf-8 -*-
"""Utility functions to be used in different modules of the chemFilters package."""

import logging
import warnings
from typing import Union

from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw


def smi_to_mol(smi: str):
    """Convert a SMILES string to a rdkit.Chem.Mol object."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        warnings.warn(f"Could not convert SMILES {smi} to rdkit.Chem.Mol")
    return mol


def get_catalog_match(mol: Chem.Mol, catalog: Chem.FilterCatalog, from_smi: bool = False):
    """Get the filter names, descriptions, and respective substructures from
    rdkit.Chem.FilterCatalog.FilterCatalogParams.GetMatches.

    Used to get results in parallel using multiprocessing.Pool.
    """
    if from_smi:
        mol = smi_to_mol(mol)
    if mol is None:
        return None, None
    matches = catalog.GetMatches(mol)
    filter_names = [m.GetProp("FilterSet") for m in matches]
    descriptions = [m.GetDescription() for m in matches]
    substructs = [m.GetFilterMatches(mol)[0].filterMatch.GetPattern() for m in matches]
    return filter_names, descriptions, substructs

def smi_to_img(smi: str, size: tuple = (300, 300)) -> Image:
    """
    Wrapper function to convert smiles to images using RDKit.

    Args:
        smi: smiles string
        size: size of the image. Default: (300, 300)

    Returns:
        RDKit image object.
    """
    mol = Chem.MolFromSmiles(smi)
    img = Draw.MolToImage(mol, size=size)
    return img


def smiles_img_matching_pose(
    smi_list: str, match_struct: Union[str, Chem.Mol], size: tuple = (300, 300)
) -> Image:
    """Wrapper function to convert smiles to images using RDKit.

    Args:
        smi: smiles string
        match_struct: used to match the pose of the molecule. Can be a RDKit Mol object,
            SMILES or SMARTS.
        size: size of the image. Default: (300, 300)

    Raises:
        ValueError: if match_struct is not a valid RDKit Mol object/SMILES/SMARTS.

    Returns:
        PIL image object.
    """
    mols = [Chem.MolFromSmiles(smi) for smi in smi_list]
    if isinstance(match_struct, Chem.rdchem.Mol):
        AllChem.Compute2DCoords(match_struct)
    else:
        mol = Chem.MolFromSmiles(match_struct)
        if mol is None:
            logging.warning('Error: "match_struct" from SMILES RDKit Mol is invalid.')
            logging.warning('Trying "match_struct" from SMARTS...')
            mol = Chem.MolFromSmarts(match_struct)
            if mol is None:
                raise ValueError('Error: "match_struct" invalid from SMILES & SMARTS.')
        match_struct = mol
    AllChem.Compute2DCoords(match_struct)
    for idx in range(0, len(smi_list)):
        AllChem.GenerateDepictionMatching2DStructure(
            mols[idx], match_struct, acceptFailure=True
        )
    imgs = [Draw.MolToImage(mol, size=size) for mol in mols]
    return imgs