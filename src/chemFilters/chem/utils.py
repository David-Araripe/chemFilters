# -*- coding: utf-8 -*-
"""Utility functions to be used by the chemFilters.chem subpackage."""

import logging

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


def molToConnectivity(smi, verbose: bool = True):
    """Converts a SMILES string to a connectivity string."""
    try:
        mol = Chem.MolFromSmiles(smi)
        connectivity = Chem.MolToInchiKey(mol).split("-")[0]
        return connectivity
    except TypeError:
        if verbose:
            logging.warn("Error converting SMILES to connectivity: ", smi)
        return None


def molToInchiKey(smi, verbose: bool = True):
    """Converts a SMILES string to an InChI string."""
    try:
        mol = Chem.MolFromSmiles(smi)
        connectivity = Chem.MolToInchiKey(mol)
        return connectivity
    except TypeError:
        if verbose:
            logging.warn("Error converting SMILES to connectivity: ", smi)
        return None


def molToInchi(smi, verbose: bool = True):
    """Converts a SMILES string to an InChI string."""
    try:
        mol = Chem.MolFromSmiles(smi)
        inchi_str = inchi.MolToInchi(mol)
        return inchi_str
    except TypeError:
        if verbose:
            logging.warn("Error converting SMILES to InChI: ", inchi)
        return None


def molToCanon(mol, isomeric: bool = True):
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=isomeric)
