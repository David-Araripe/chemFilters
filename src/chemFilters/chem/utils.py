# -*- coding: utf-8 -*-
"""Utility functions to be used by the chemFilters.chem subpackage."""

from loguru import logger
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi


def RDKitVerbosityON():
    """Switches on RDKit verbosity"""
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.INFO)
    return lg


def RDKitVerbosityOFF():
    """Switches off RDKit verbosity"""
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    return lg


def molToConnectivity(mol: Chem.Mol):
    """Converts a SMILES string to a connectivity string."""
    try:
        connectivity = Chem.MolToInchiKey(mol).split("-")[0]
    except TypeError as e:
        logger.error(f"Could not convert mol to connectivity. Message:\n {e}")
        connectivity = None
    return connectivity


def molToInchiKey(mol: Chem.Mol):
    """Converts a SMILES string to an InChI string."""
    try:
        inchi_key = Chem.MolToInchiKey(mol)
    except TypeError as e:
        logger.error(f"Could not convert mol to inchikey. Message:\n {e}")
        inchi_key = None
    return inchi_key


def molToInchi(mol: Chem.Mol):
    """Converts a SMILES string to an InChI string."""
    try:
        inchi_str = inchi.MolToInchi(mol)
    except TypeError as e:
        logger.error(f"Could not convert mol to inchi. Message:\n {e}")
        inchi_str = None
    return inchi_str


def molToCanon(mol, isomeric: bool = True):
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=isomeric)
