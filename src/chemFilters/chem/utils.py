# -*- coding: utf-8 -*-
"""Utility functions to be used by the chemFilters.chem subpackage."""

from contextlib import contextmanager

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import inchi

from ..logger import logger


@contextmanager
def rdkit_log_controller(level):
    """Context manager for controlling the RDKit logger level.

    Args:
        level: desired logging level. One of `debug, info, warning, error, critical`.
    """
    lg = RDLogger.logger()
    lvl = level.upper()
    # in case no logging is enable, it will be set to critical
    lvl_list = rdBase.LogStatus().split("\n") + ["rdApp.critical:enabled"]
    lg.setLevel(getattr(RDLogger, lvl))
    init_lvl = [lev.split(":")[0] for lev in lvl_list if "enabled" in lev.lower()][0]
    yield
    lg.setLevel(getattr(RDLogger, init_lvl.split(".")[1].upper()))  # restore the level


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
