# -*- coding: utf-8 -*-
"""Utility functions to be used by the chemFilters.chem subpackage."""

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
    connectivity = Chem.MolToInchiKey(mol).split("-")[0]
    return connectivity


def molToInchiKey(mol: Chem.Mol):
    """Converts a SMILES string to an InChI string."""
    connectivity = Chem.MolToInchiKey(mol)
    return connectivity


def molToInchi(mol: Chem.Mol):
    """Converts a SMILES string to an InChI string."""
    inchi_str = inchi.MolToInchi(mol)
    return inchi_str


def molToCanon(mol, isomeric: bool = True):
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=isomeric)
