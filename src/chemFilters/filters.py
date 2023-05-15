# -*- coding: utf-8 -*-
from multiprocessing import Pool
from typing import List

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


class RDKitFilters:
    def __init__(self, filterType="ALL", njobs=1) -> None:
        self.filterType = filterType
        self.filter = self._getFilter()
        self.njobs = njobs
        self.filterParams = FilterCatalog(FilterCatalogParams().AddCatalog(self.filter))
        pass

    @property
    def availableFilters(self):
        return [m for m in dir(FilterCatalogParams.FilterCatalogs) if m.isupper()]

    def _getFilter(self):
        if self.filterType not in self.availableFilters:
            raise ValueError(f"Filter type {self.filterType} not available.")
        return getattr(FilterCatalogParams.FilterCatalogs, self.filterType)

    def _filterFunc(self, mol: rdkit.Chem.Mol):
        matches = self.filterParams.GetMatches(mol)
        filterMatches = self.filterParams.GetFilterMatches(mol)

        description = [m.GetDescription() for m in matches]
        substructs = [m.GetFilterMatches() for m in filterMatches]
        return description, substructs

    def filterCompounds(self, mols: List[rdkit.Chem.Mol]):
        with Pool(self.njobs) as p:
            p.map()

    def lazyFilterCompounds(self, mols: List[rdkit.Chem.Mol]):
        with Pool(self.njobs) as p:
            p.imap()

    @staticmethod
    def lemmeSeeMatches(mol, description, substructs):
        hit_bonds = []
        hit_atoms = []
        for patt in substructs:
            smarts = Chem.MolToSmarts(patt)
            pttrn = Chem.MolFromSmarts(smarts)
            hit_ats = list(mol.GetSubstructMatch(pttrn))
            hit_atoms = hit_atoms + hit_ats

            for bond in pttrn.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())

        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(
            mol, size=(500, 500), highlightAtoms=hit_atoms, highlightBonds=hit_bonds
        )
        return img

    @staticmethod
    def seeSubstructs(mol, description, substructs):
        # This method would just show the substructures as a molstogrid
        pass
