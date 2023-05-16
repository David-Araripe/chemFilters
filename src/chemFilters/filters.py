# -*- coding: utf-8 -*-
from io import BytesIO
from multiprocessing import Pool
from typing import List

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from PIL import Image
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


class RDKitFilters:
    def __init__(self, filterType="ALL", njobs=1) -> None:
        self.filterType = filterType
        self.filter = self._getFilter()
        self.njobs = njobs

    @property
    def availableFilters(self):
        """List of available filters from RDKit FilterCatalogs."""
        return [m for m in dir(FilterCatalogParams.FilterCatalogs) if m.isupper()]

    @property
    def filterParams(self):
        catalog = FilterCatalogParams()
        catalog.AddCatalog(self.filter)
        return FilterCatalog(catalog)

    def _getFilter(self):
        if self.filterType not in self.availableFilters:
            raise ValueError(f"Filter type {self.filterType} not available.")
        return getattr(FilterCatalogParams.FilterCatalogs, self.filterType)

    def _filterFunc(self, mol: Chem.Mol):
        matches = self.filterParams.GetMatches(mol)
        filterMatches = self.filterParams.GetFilterMatches(mol)

        description = [m.GetDescription() for m in matches]
        substructs = [m.filterMatch.GetPattern() for m in filterMatches]
        return description, substructs

    def filterCompounds(self, mols: List[Chem.Mol]):
        with Pool(self.njobs) as p:
            p.map()

    def lazyFilterCompounds(self, mols: List[Chem.Mol]):
        with Pool(self.njobs) as p:
            p.imap()

    @staticmethod
    def renderColoredMatches(
        mol: Chem.Mol,
        description,
        substructs,
        cmap="rainbow",
        alpha=0.5,
        ax=None,
        molSize=(500, 500),
    ):
        """Take descriptions and substructures output from RDKitFilters._filterFunc
        and renders them on the molecular structure, coloring the substructures and
        adding labels according to the descriptions.

        Note: in case of overlap between colors, the final color corresponds to the
        their respective geometric mean.

        Args:
            mol: rdkit molecule object.
            description: description (`output[0]`) from RDKitFilters._filterFunc.
            substructs: substructures (`output[1]`) from RDKitFilters._filterFunc.
            cmap: colormap for the substructures. Defaults to "viridis".
            alpha: transparency of the colors. Defaults to 0.5.
            ax: if desired, add the results to an existing axis. Defaults to None.
            molSize: size of the molecule rendering. Defaults to (500, 500).
                Note: in the end this will be bound to the size of the plot...

        Returns:
            matplotlib.figure.Figure and Axis with the molecule and the
                highlighted substructures.
        """

        def geometric_mean(arr, axis=0):
            """Calculate the geometric mean of an array. Adapted from SciPy:
            https://github.com/scipy/blob/v1.10.1/scipy/stats/_stats.py.py#L199=L269"""
            with np.errstate(divide="ignore"):
                log_a = np.log(arr)
            return np.exp(np.average(log_a, axis=axis))

        qualitative_cmaps = [
            "Pastel1",
            "Pastel2",
            "Paired",
            "Accent",
            "Dark2",
            "Set1",
            "Set2",
            "Set3",
            "tab10",
            "tab20",
            "tab20b",
            "tab20c",
        ]

        smarts = [Chem.MolToSmarts(sub) for sub in substructs]
        unique_smarts, indices = np.unique(smarts, return_index=True)
        unique_descrip = [description[i] for i in indices]

        if cmap in qualitative_cmaps:
            # Unpack colors and add alpha
            colors = [tuple([*col] + [alpha]) for col in cm.get_cmap("tab10").colors]
        else:
            scalar_mappable = cm.ScalarMappable(cmap=cmap)
            colors = scalar_mappable.to_rgba(
                range(len(unique_smarts)), alpha=alpha
            ).tolist()
        color_dict = {}
        patches = []

        drawer = rdMolDraw2D.MolDraw2DCairo(molSize[0], molSize[1])
        if ax is None:
            fig, ax = plt.subplots(figsize=(5, 5))
        else:
            ax = ax.gca()

        for smarts, descr, color in zip(unique_smarts, unique_descrip, colors):
            # Create the patches that are used for the labels
            patches.append(
                mpatches.Patch(facecolor=color, label=descr, edgecolor="black")
            )
            pttrn = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatch(pttrn)
            for match in matches:
                if match in color_dict.keys():
                    color_dict.update({match: np.vstack([color_dict[match], color])})
                else:
                    color_dict.update({match: np.array(color)})
        color_dict = {
            k: tuple(geometric_mean(v)) if v.ndim > 1 else tuple(v)
            for k, v in color_dict.items()
        }  # should be tuples
        drawer.DrawMolecule(
            mol, highlightAtoms=color_dict.keys(), highlightAtomColors=color_dict
        )
        drawer.FinishDrawing()
        img = Image.open(BytesIO(drawer.GetDrawingText()))
        ax.legend(
            handles=patches,
            bbox_to_anchor=(1.04, 0.2),
            loc="lower left",
            borderaxespad=0,
            frameon=False,
        )
        ax.imshow(img)
        ax.set_axis_off()
        plt.close()
        return fig, ax

    @staticmethod
    def renderMatches(mol: Chem.Mol, substructs, molSize=(500, 500)) -> Image:
        """Render the substructures on the molecules using default RDKit coloring.

        Args:
            mol: molecule to be rendered.
            substructs: substructures output from RDKitFilters._filterFunc.
            molSize: final image size. Defaults to (500, 500).

        Returns:
            PIL.Image with the molecule and the highlighted substructures.
        """
        hit_bonds = []
        hit_atoms = []
        if "__iter__" not in dir(substructs):
            substructs = [substructs]
        for patt in substructs:
            smarts = Chem.MolToSmarts(patt)
            pttrn = Chem.MolFromSmarts(smarts)
            hit_ats = list(mol.GetSubstructMatch(pttrn))
            hit_atoms = hit_atoms + hit_ats

            for bond in pttrn.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())

        drawer = rdMolDraw2D.MolDraw2DCairo(molSize[0], molSize[1])
        drawer.DrawMolecule(mol, highlightAtoms=hit_atoms, highlightBonds=hit_bonds)
        drawer.FinishDrawing()
        img = Image.open(BytesIO(drawer.GetDrawingText()))
        return img
