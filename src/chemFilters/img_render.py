# -*- coding: utf-8 -*-
"""Contains classes for finding custom fonts & rendering molecules into figures."""

import os
import re
import sys
from functools import partial
from io import BytesIO
from multiprocessing import Pool
from pathlib import Path
from typing import List, Optional, Tuple, Union

import matplotlib.font_manager as fm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pkg_resources
from loguru import logger
from matplotlib import cm
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from .chem.interface import MoleculeHandler


class FontManager:
    def __init__(self, wsl: Union[bool, str] = "auto"):
        """Initialize the FontManager class. It will search for fonts in the system and
        in matplotlib.

        Args:
            wsl: If you're using windows subsystem for linux. Allows for using fonts
                installed in windows. Defaults to "auto".
            include_plt: Includes fonts available in matplotlib. Defaults to True.
        """
        self.wsl = wsl
        self.fonts = None
        pass

    @property
    def operating_system(self):
        if os.name == "nt":
            return "windows"
        elif os.name == "posix":
            if sys.platform.startswith("linux"):
                if self.wsl:
                    wsl_pattern = re.compile(r"wsl|windows", re.IGNORECASE)
                    with open("/proc/version", "r") as f:
                        content = f.read()
                    if wsl_pattern.findall(content):
                        return "wsl"
                return "linux"
            elif sys.platform.startswith("darwin"):
                return "macos"
        return "Unknown"

    @property
    def available_fonts(self):
        font_dict = dict()
        font_extensions = "[.ttf|.otf|.eot|.woff|.woff2|.svg]"

        if self.operating_system == "windows":
            default_font_dir = Path("C:/Windows/Fonts")
            if default_font_dir.exists():
                font_paths = sorted(list(default_font_dir.glob(f"*{font_extensions}")))
                font_dict.update({_path.stem: _path for _path in font_paths})
            font_dir = Path.home() / "AppData/Local/Microsoft/Windows/Fonts"
            font_paths = sorted(list(font_dir.glob(f"*{font_extensions}")))
            font_dict.update({_path.stem: _path for _path in font_paths})

        elif self.operating_system == "wsl":
            # Search within the default windows path
            default_font_dir = Path("/mnt/c/Windows/Fonts")
            if default_font_dir.exists():
                font_paths = sorted(list(default_font_dir.glob(f"*{font_extensions}")))
                font_dict.update({_path.stem: _path for _path in font_paths})
            user_dir = Path("/mnt/c/Users/")
            font_dir = list(user_dir.glob("*/AppData/Local/Microsoft/Windows/Fonts"))
            if len(font_dir) > 1:
                print("More than one user font path found.")
            for _dir in font_dir:
                font_paths = sorted(list(_dir.glob(f"*{font_extensions}")))
                font_dict.update({_path.stem: _path for _path in font_paths})

        elif self.operating_system == "linux":
            font_dir = [
                Path("/usr/share/fonts"),
                Path("/usr/local/share/fonts"),
                Path("~/.fonts").expanduser(),
            ]
            for _dir in font_dir:
                if _dir.exists():
                    font_paths = sorted(list(_dir.glob(f"*{font_extensions}")))
                    font_dict.update({_path.stem: _path for _path in font_paths})

        elif self.operating_system == "macos":
            font_dirs = [
                Path("/System/Library/Fonts"),
                Path("/Library/Fonts"),
                Path("~/Library/Fonts").expanduser(),
            ]
            for _dir in font_dirs:
                if _dir.exists():
                    font_paths = sorted(list(_dir.glob(f"*{font_extensions}")))
                    font_dict.update({_path.stem: _path for _path in font_paths})

        # add fonts from matplotlib
        font_dict.update({font.name: font.fname for font in fm.fontManager.ttflist})
        # add fonts from rdkit
        font_dir = Path(pkg_resources.resource_filename("rdkit", "Data/Fonts"))
        font_dict.update({_path.stem: _path for _path in font_dir.glob("*.ttf")})
        return font_dict


class MolPlotter(MoleculeHandler):
    """MolPlotter class to render molecules into images. It uses RDKit's
    DrawMolecule method to generate the images. The images can be rendered with
    labels and/or substructure matches.

    Attributes:
        from_smi: output images directly from smiles. Defaults to False.
        size: size of the output images. Defaults to (300, 300).
        cmap: colormap for `render_with_color` method. Defaults to "rainbow".
        font_name: font name found by img_render.FondManager class.
            Defaults to "DejaVu Sans".
        font_size: size of the font patched to the image. Defaults to 15.
        n_jobs: number of jobs to run in parallel. Defaults to 1.
        d2d_params: dictionary with the parameters to be passed to get_d2d method.
        available_fonts: fonts available in the system and in matplotlib.
    """

    def __init__(
        self,
        from_smi: bool = False,
        size: tuple = (300, 300),
        cmap: str = "rainbow",
        font_name: str = "Telex-Regular",
        font_size: int = 15,
        n_jobs: int = 1,
        bond_line_width: float = 2.0,
        add_atom_indices: bool = False,
        add_bond_indices: bool = False,
        explicit_methyl: bool = False,
        unspecified_stereo_unknown: bool = False,
        bw: bool = False,
    ) -> None:
        """Initialize the MolPlotter class. The class is used to render molecules into
        either .sdf or .png files. It uses RDKit's DrawMolecule method to generate the
        images and `self.d2d_params` to set the parameters for the rendering process.

        Args:
            from_smi: if true, methods will treat inputs as smiles. Defaults to False.
            size: size of the image to de displayed. Defaults to (300, 300).
            cmap: colormap for the structure matching. Defaults to "rainbow".
            font_name: font used on molecules' legends. Defaults to "DejaVu Sans".
            font_size: size of the font on the legend. Defaults to 15.
            n_jobs: number of jobs to run in parallel. Defaults to 1.
            bond_line_width (MolDraw2DCairo): width of bond lines. Defaults to 2.0.
            add_atom_indices (MolDraw2DCairo): add atom indices to the rendered mols.
                Defaults to False.
            add_bond_indices (MolDraw2DCairo): add bond indices to the rendered mols.
                Defaults to False.
            explicit_methyl (MolDraw2DCairo): explicitly displays a methil as CH3.
                Defaults to False.
            unspecified_stereo_unknown (MolDraw2DCairo): unspecified stereo atoms/bonds
                are drawn as if they were unknown. Defaults to False.
            bw: render the molecule as black and white. Defaults to False.
        """
        self._cmap = cmap
        self._size = size
        self._font_size = font_size
        self._font_name = font_name
        self._n_jobs = n_jobs
        self.d2d_params = {
            "bondLineWidth": bond_line_width,
            "addAtomIndices": add_atom_indices,
            "addBondIndices": add_bond_indices,
            "explicitMethyl": explicit_methyl,
            "unspecifiedStereoIsUnknown": unspecified_stereo_unknown,
            "bw": bw,
        }
        self._check_font(font_name)
        super().__init__(from_smi)

    @property
    def available_fonts(self):
        fm = FontManager()
        return fm.available_fonts

    def _check_font(self, font_name):
        if self._font_name not in self.available_fonts.keys():
            logger.warning(
                f"Font {self._font_name} not found. Check `available_fonts` attribute."
            )
        return font_name

    @staticmethod
    def geometric_mean(arr, axis=0):
        """Calculate the geometric mean of an array. Adapted from SciPy:
        https://github.com/scipy/blob/v1.10.1/scipy/stats/_stats.py.py#L199=L269"""
        with np.errstate(divide="ignore"):
            log_a = np.log(arr)
        return np.exp(np.average(log_a, axis=axis))

    def get_d2d(
        self,
        bondLineWidth: float = 2.0,
        addAtomIndices: bool = False,
        addBondIndices: bool = False,
        explicitMethyl: bool = False,
        unspecifiedStereoIsUnknown: bool = False,
        bw: bool = False,
    ) -> Draw.MolDraw2DCairo:
        """Function to get the rdkit's MolDraw2DCairo object with the desired options.

        Args:
            bondLineWidth: width of the bonds in the drawing. RDKit defaults to 2.0.
            addAtomIndices: add the atom indices to the rendered mol. Defaults to False.
            addBondIndices: add the bond indices to the rendered mol. Defaults to False.
            explicitMethyl: explicitly displays a methil as CH3. Defaults to False.
            unspecifiedStereoIsUnknown: unspecified stereo atoms/bonds are drawn as if
                they were unknown. Defaults to False.
            bw: render the molecule as black and white. Defaults to False.

        Returns:
            Draw.MolDraw2DCairo object with the desired options.
        """
        d2d = Draw.MolDraw2DCairo(*self._size)
        dopts = d2d.drawOptions()
        dopts.fixedFontSize = self._font_size
        dopts.bondLineWidth = bondLineWidth
        if addAtomIndices:
            dopts.noAtomLabels = True
        if addBondIndices:
            dopts.addBondIndices = True
        if explicitMethyl:
            dopts.explicitMethyl = True
        if unspecifiedStereoIsUnknown:
            dopts.unspecifiedStereoIsUnknown = True
        if bw:
            dopts.useBWAtomPalette()
        font_path = self.available_fonts.get(self._font_name)
        dopts.fontFile = str(font_path)
        return d2d

    def _stdin_to_mol(self, stdin):
        if self._from_smi:
            return self._output_mol(stdin)
        else:
            return stdin

    def _substructure_palette(
        self, n_substruct: int, cmap: str, alpha: float = 1.0
    ) -> List[Tuple[float, float, float, float]]:
        """ """
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
        if cmap in qualitative_cmaps:
            # Unpack colors and add alpha
            colors = [tuple([*col] + [alpha]) for col in cm.get_cmap(cmap).colors]
        else:
            scalar_mappable = cm.ScalarMappable(cmap=cmap)
            colors = scalar_mappable.to_rgba(range(n_substruct), alpha=alpha).tolist()
        return colors

    def _images_to_grid(self, images: List[Image.Image], n_cols: int) -> Image.Image:
        """Helper function to organize a list of images into a single image, as a grid.

        Args:
            images: Images to be added into the grid.
            n_cols: number of columns of the grid.

        Returns:
            The image collage as a PIL.Image.Image object.
        """

        def list_split(list_a, chunk_size):
            """
            Returns a generator with a determined chunk size over list_a.
            Used to organize the figures in the correct `n_cols` format.
            """
            for i in range(0, len(list_a), chunk_size):
                yield list_a[i : i + chunk_size]

        # Splitting the list of images into a list of lists with n_cols
        if n_cols is None:
            n_cols = len(images)
        images = list(list_split(images, n_cols))
        # Appending blank figures so we have the correct vector shapes
        while len(images[-1]) < len(images[0]):
            images[-1].append(Image.new("RGB", self._size, color=(255, 255, 255)))

        list_of_hstacked = list()
        # Creating list of horizontally stacked arrays
        for sublist in images:
            list_of_hstacked.append(np.hstack([np.asarray(img) for img in sublist]))
        # Vertically stacking horizontal arrays
        for _item in list_of_hstacked:
            final_img = np.vstack([hstack for hstack in list_of_hstacked])
        # Creating and returning image from array
        final_img = Image.fromarray(final_img)
        return final_img

    def render_ACS1996(
        self,
        mol: Chem.Mol,
        mean_bond_length: float = None,
        label: str = "",
        return_svg=False,
        **kwargs,
    ) -> Union[Image.Image, str]:
        """Render molecules using the ACS 1996 style. This is a mode which is designed
        to produce images compatible with the drawing standards for American Chemical
        Society (ACS) journals.

        Args:
            mol: rdkit mol object.
            mean_bond_length: the mean bond length of the molecule. Defaults to None.
            label: the label to be added to the molecule. Defaults to "".
            return_svg: whether to return an svg image instead. Defaults to False.
            kwargs: additional keyword arguments to be passed to the DrawMolecule method

        Returns:
            a PIL.Image.Image object or a string with the SVG representation in case
                return_svg is True.
        """
        d2d = self.get_d2d()
        Chem.rdDepictor.Compute2DCoords(mol)
        if mean_bond_length is None:
            mean_bond_length = Draw.MeanBondLength(mol)
        Draw.SetACS1996Mode(d2d.drawOptions(), mean_bond_length)
        d2d.DrawMolecule(mol, legend=label, **kwargs)
        d2d.FinishDrawing()
        if not return_svg:
            img = Image.open(BytesIO(d2d.GetDrawingText()))
        else:
            img = d2d.GetDrawingText()
        return img

    def render_mol(
        self,
        mol: Chem.Mol,
        label: str = None,
        match_pose: Optional[Union[str, Chem.Mol]] = None,
        return_svg: bool = False,
        **kwargs,
    ) -> Union[Image.Image, str]:
        """Plot molecules using rdMolDraw2D.MolDraw2DCairo. Keyword arguments provided
        will be passed into DrawMolecule method of the same class.

        Args:
            mol: rdkit mol object.
            label: the label to be added to the molecule. Defaults to None.
            match_pose: if desired, input a SMILES, a SMARTS or a rdkit mol object to
                render `mol` into a matching position. Defaults to None.
            return_svg: if you'd like to return an svg image instead. Defaults to False.
            kwargs: additional keyword arguments to be passed to the DrawMolecule method

        Raises:
            ValueError: _description_

        Returns:
            a PIL.Image.Image object or a string with the SVG representation in case
                return_svg is True.
        """
        if label is None:
            label = ""
        mol = self._stdin_to_mol(mol)
        if match_pose is not None:
            if isinstance(match_pose, Chem.rdchem.Mol):
                AllChem.Compute2DCoords(match_pose)
            else:
                ref_mol = Chem.MolFromSmiles(match_pose)
                if ref_mol is None:
                    logger.warning(
                        'Error: "match_pose" from SMILES RDKit Mol is invalid.'
                    )
                    logger.warning('Trying "match_pose" from SMARTS...')
                    ref_mol = Chem.MolFromSmarts(match_pose)
                    if ref_mol is None:
                        raise ValueError(
                            'Error: "match_pose" invalid from SMILES & SMARTS.'
                        )
                match_pose = ref_mol
            AllChem.Compute2DCoords(match_pose)
            AllChem.GenerateDepictionMatching2DStructure(
                mol, match_pose, acceptFailure=True
            )
        d2d = self.get_d2d(**self.d2d_params)
        d2d.DrawMolecule(mol, legend=label, **kwargs)
        d2d.FinishDrawing()
        if not return_svg:
            img = Image.open(BytesIO(d2d.GetDrawingText()))
        else:
            img = d2d.GetDrawingText()
        return img

    def render_with_matches(
        self,
        mol: Chem.Mol,
        substructs,
        label: Optional[str] = None,
        match_pose: Optional[Union[str, Chem.Mol]] = None,
        **kwargs,
    ) -> Image.Image:
        """Render the substructures on the molecules using default RDKit coloring.

        Args:
            mol: molecule to be rendered.
            substructs: substructures output from RDKitFilters.filter_mols.
            label: if desired, add a label to the molecule. Defaults to None.
            match_pose: a common structure between `smiles` used to generate the
                matching 2D images. If None, 2D rendering won't match. Defaults to None.
            kwargs: additional keyword arguments to be passed to the DrawMolecule method

        Returns:
            PIL.Image.Image with the molecule and the highlighted substructures.
        """
        mol = self._output_mol(mol)
        hit_bonds = []
        hit_atoms = []
        if "__iter__" not in dir(substructs):
            substructs = [substructs]
        for patt in substructs:
            smarts = Chem.MolToSmarts(patt)
            pttrn = Chem.MolFromSmarts(smarts)
            hit_ats = list(mol.GetSubstructMatch(pttrn))
            if hit_ats == []:
                continue
            hit_atoms = hit_atoms + hit_ats
            for bond in pttrn.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())
        img = self.render_mol(
            mol,
            highlightAtoms=hit_atoms,
            highlightBonds=hit_bonds,
            label=label,
            match_pose=match_pose,
            **kwargs,
        )
        return img

    def render_with_colored_matches(
        self,
        mol: Chem.Mol,
        descriptions: List[str],
        substructs: List[Chem.Mol],
        label: Optional[str] = None,
        cmap: str = "rainbow",
        alpha: float = 0.5,
        match_pose: Optional[Union[str, Chem.Mol]] = None,
        **kwargs,
    ):
        """Take descriptions and substructures output from RDKitFilters.filter_mols
        and renders them on the molecular structure, coloring the substructures and
        adding labels according to the descriptions.

        Note: in case of overlap between colors, the final color corresponds to the
        their respective geometric mean.

        Args:
            mol: rdkit molecule object.
            descriptions: descriptions (`output[0]`) from RDKitFilters.filter_mols.
            substructs: substructures (`output[1]`) from RDKitFilters.filter_mols.
            label: if desired, add a label to the molecule. Defaults to None.
            cmap: colormap for the substructures. Defaults to "viridis".
            alpha: transparency of the colors. Defaults to 0.5.
            match_pose: a common structure between `smiles` used to generate the
                matching 2D images. If None, 2D rendering won't match. Defaults to None.
            kwargs: additional keyword arguments to be passed to the DrawMolecule method

        Returns:
            matplotlib.figure.Figure and Axis with the molecule and the
                highlighted substructures.
        """
        mol = self._output_mol(mol)
        smarts = [Chem.MolToSmarts(sub) for sub in substructs]
        unique_smarts, indices = np.unique(smarts, return_index=True)
        unique_descrip = [descriptions[i] for i in indices]

        color_dict = {}
        colors = self._substructure_palette(len(unique_smarts), cmap=cmap, alpha=alpha)
        for smarts, _descr, color in zip(unique_smarts, unique_descrip, colors):
            # Create the patches that are used for the labels
            pttrn = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatch(pttrn)
            for match in matches:
                if match in color_dict.keys():
                    color_dict.update({match: np.vstack([color_dict[match], color])})
                else:
                    color_dict.update({match: np.array(color)})
        color_dict = {
            k: tuple(self.geometric_mean(v)) if v.ndim > 1 else tuple(v)
            for k, v in color_dict.items()
        }  # should be tuples
        img = self.render_mol(
            mol,
            label=label,
            highlightAtoms=color_dict.keys(),
            highlightAtomColors=color_dict,
            match_pose=match_pose,
            **kwargs,
        )
        return img

    def colored_matches_legend(
        self, descriptions, substructures, cmap="rainbow", alpha=0.5, ax=None
    ):
        if ax is not None:
            ax = plt.gca()
        else:
            fig, ax = plt.subplots(figsize=(5, 5))
        patches = []
        colors = self._substructure_palette(len(substructures), cmap=cmap, alpha=alpha)
        for _smarts, descr, color in zip(substructures, descriptions, colors):
            patches.append(
                mpatches.Patch(facecolor=color, label=descr, edgecolor="black")
            )
        ax.legend(
            handles=patches,
            bbox_to_anchor=(1.05, 0.25),
            loc="lower left",
            borderaxespad=0,
            frameon=False,
        )
        ax.set_axis_off()
        if ax is None:
            return fig, ax


class MolGridPlotter(MolPlotter):
    """MolGridPlotter class to render molecules into grids of images. It uses RDKit's
    DrawMolecule method to generate the images and PIL to patch images together.
    The images can be rendered with labels and/or substructure matches, similar to its
    parent class.

    Attributes:
        from_smi: output images directly from smiles. Defaults to False.
        size: size of the output images. Defaults to (300, 300).
        cmap: colormap for `render_with_color` method. Defaults to "rainbow".
        font_name: font name found by img_render.FondManager class.
            Defaults to "DejaVu Sans".
        font_size: size of the font patched to the image. Defaults to 15.
        n_jobs: number of jobs to run in parallel. Defaults to 1.
        d2d_params: dictionary with the parameters to be passed to get_d2d method.
        available_fonts: fonts available in the system and in matplotlib.
    """

    def __init__(
        self,
        from_smi: bool = False,
        size: Tuple = (300, 300),
        cmap: str = "rainbow",
        font_name: str = "Telex-Regular",
        font_size: int = 15,
        n_jobs: int = 1,
        bond_line_width: float = 2,
        add_atom_indices: bool = False,
        add_bond_indices: bool = False,
        explicit_methyl: bool = False,
        unspecified_stereo_unknown: bool = False,
        bw: bool = False,
    ) -> None:
        """initialize the MolGridPlotter class. The class is used to render molecules
        into either .sdf or .png files. It uses RDKit's DrawMolecule method to generate
        the images and `self.d2d_params` to set the parameters for the rendering process

        Args:
            from_smi: if true, methods will treat inputs as smiles. Defaults to False.
            size: size of the image to de displayed. Defaults to (300, 300).
            cmap: colormap for the structure matching. Defaults to "rainbow".
            font_name: font used on molecules' legends. Defaults to "DejaVu Sans".
            font_size: size of the font on the legend. Defaults to 15.
            n_jobs: number of jobs to run in parallel. Defaults to 1.
            bond_line_width (MolDraw2DCairo): width of bond lines. Defaults to 2.0.
            add_atom_indices (MolDraw2DCairo): add atom indices to the rendered mols.
                Defaults to False.
            add_bond_indices (MolDraw2DCairo): add bond indices to the rendered mols.
                Defaults to False.
            explicit_methyl (MolDraw2DCairo): explicitly displays a methil as CH3.
                Defaults to False.
            unspecified_stereo_unknown (MolDraw2DCairo): unspecified stereo atoms/bonds
                are drawn as if they were unknown. Defaults to False.
            bw: render the molecule as black and white. Defaults to False.
        """
        super().__init__(
            from_smi,
            size,
            cmap,
            font_name,
            font_size,
            n_jobs,
            bond_line_width,
            add_atom_indices,
            add_bond_indices,
            explicit_methyl,
            unspecified_stereo_unknown,
            bw,
        )

    def mol_grid_png(
        self,
        mols: list,
        labels: Optional[list] = None,
        n_cols: int = None,
        match_pose: Optional[Union[str, Chem.Mol]] = None,
        **kwargs,
    ) -> Image.Image:
        """
        Args:
            mols: list of molecules to display.
            labels: list of labels to be written on the mol images. Defaults to None
            n_cols: N columns with molecules. If None, n_cols = len(smiles).
                Defaults to None.
            match_pose: a common structure between `smiles` used to generate the
                matching 2D images. If None, 2D rendering won't match. Defaults to None.
            kwargs: additional keyword arguments to be passed to the DrawMolecule method

        Returns:
            PIL.Image.Image
        """
        if labels is not None:
            assert len(labels) == len(mols), "Labels and mols must have the same length"
            variables = list(zip(mols, labels))
        else:
            variables = list(zip(mols))
        with Pool(self._n_jobs) as p:
            images = p.starmap(
                partial(self.render_mol, match_pose=match_pose, **kwargs),
                variables,
            )
        image = self._images_to_grid(images, n_cols)
        return image

    def mol_structmatch_grid_png(
        self,
        mols: list,
        substructs: list,
        labels: Optional[list] = None,
        n_cols: int = None,
        match_pose: Optional[Union[str, Chem.Mol]] = None,
        **kwargs,
    ) -> Image.Image:
        """
        Args:
            mols: list of molecules to display.
            substructs: list of substructures to be highlighted on the mol images.
            labels: list of labels to be written on the mol images. Defaults to None
            n_cols: N columns with molecules. If None, n_cols = len(smiles).
                Defaults to None.
            match_pose: a common structure between `smiles` used to generate the
                matching 2D images. If None, 2D rendering won't match. Defaults to None.
            kwargs: additional keyword arguments to be passed to the DrawMolecule method

        Returns:
            PIL.Image.Image
        """
        partial_function = partial(
            self.render_with_matches, match_pose=match_pose, **kwargs
        )
        if labels is not None:
            assert len(labels) == len(mols), "Labels and mols must have the same length"
            variables = list(zip(mols, substructs, labels))
        else:
            variables = list(zip(mols, substructs))
        with Pool(self._n_jobs) as p:
            images = p.starmap(partial_function, variables)
        image = self._images_to_grid(images, n_cols)
        return image

    def mol_structmatch_color_grid_png(
        self,
        mols: list,
        descriptions: List[str],
        substructs: List[Chem.Mol],
        labels: Optional[list] = None,
        n_cols: int = None,
        match_pose: Optional[Union[str, Chem.Mol]] = None,
        **kwargs,
    ) -> Image.Image:
        """
        Args:
            mols: list of molecules to display.
            descriptions: descriptions (`output[0]`) from RDKitFilters.filter_mols
            substructs: substructures (`output[1]`) from RDKitFilters.filter_mols.
            labels: list of labels to be written on the mol images. Defaults to None
            n_cols: N columns with molecules. If None, n_cols = len(smiles).
                Defaults to None.
            match_pose: a common structure between `smiles` used to generate the
                matching 2D images. If None, 2D rendering won't match. Defaults to None.
            kwargs: additional keyword arguments to be passed to the DrawMolecule method

        Returns:
            PIL.Image.Image
        """
        mols = self._stdin_to_mol(mols)
        partial_function = partial(
            self.render_with_colored_matches, match_pose=match_pose, **kwargs
        )
        if labels is not None:
            assert len(labels) == len(mols), "Labels and mols must have the same length"
            variables = list(zip(mols, descriptions, substructs, labels))
        else:
            variables = list(zip(mols, descriptions, substructs))
        with Pool(self._n_jobs) as p:
            images = p.starmap(partial_function, variables)
        image = self._images_to_grid(images, n_cols)
        return image
