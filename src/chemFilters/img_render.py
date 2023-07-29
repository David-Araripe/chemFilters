# -*- coding: utf-8 -*-
"""Contains classes for finding custom fonts & rendering molecules into figures."""

import importlib
import logging
import os
import re
import sys
from functools import partial
from io import BytesIO
from multiprocessing import Pool
from pathlib import Path
from typing import List, Optional, Tuple, Union

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import PIL
from matplotlib import cm
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from .chem.interface import MoleculeHandler


class FontManager:
    def __init__(self, wsl: Union[bool, str] = "auto", include_plt: bool = True):
        """Initialize the FontManager class. It will search for fonts in the system and
        in matplotlib.

        Args:
            wsl: If you're using windows subsystem for linux. Allows for using fonts
                installed in windows. Defaults to "auto".
            include_plt: Includes fonts available in matplotlib. Defaults to True.
        """
        self.include_plt = include_plt
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

        if self.include_plt:
            import matplotlib.font_manager as fm

            font_dict.update({font.name: font.fname for font in fm.fontManager.ttflist})
        return font_dict


class MolPlotter(MoleculeHandler):
    """MolPlotter class to render molecules into images. It uses RDKit's
    DrawMolecule method to generate the images. The images can be rendered with
    labels and/or substructure matches.

    Args:
        from_smi: output images directly from smiles. Defaults to False.
        size: size of the output images. Defaults to (500, 500).
        cmap: colormap for `render_with_color` method. Defaults to "rainbow".
        font_name: font name found by img_render.FondManager class.
            Defaults to "DejaVu Sans".
        font_size: size of the font patched to the image. Defaults to 15.
        label_loc: font placement; either `top` or `bottom`. Defaults to "bottom".
        n_jobs: number of jobs to run in parallel. Defaults to 1.
    """

    def __init__(
        self,
        from_smi: bool = False,
        size: tuple = (500, 500),
        cmap: str = "rainbow",
        font_name: str = "DejaVu Sans",
        font_size: int = 15,
        label_loc="bottom",
        n_jobs: int = 1,
    ) -> None:
        self._cmap = cmap
        self._size = size
        self._font_size = font_size
        self._font_name = font_name
        self._label_loc = label_loc
        self._n_jobs = n_jobs
        super().__init__(from_smi)

    @property
    def available_fonts(self):
        try:
            importlib.import_module("matplotlib")
            include_plt = True
        except ModuleNotFoundError:
            include_plt = False
        fm = FontManager(include_plt=include_plt)
        return fm.available_fonts

    @staticmethod
    def geometric_mean(arr, axis=0):
        """Calculate the geometric mean of an array. Adapted from SciPy:
        https://github.com/scipy/blob/v1.10.1/scipy/stats/_stats.py.py#L199=L269"""
        with np.errstate(divide="ignore"):
            log_a = np.log(arr)
        return np.exp(np.average(log_a, axis=axis))

    def _stdin_to_mol(self, stdin):
        if self._from_smi:
            return self._output_mol(stdin)
        else:
            return stdin

    def _add_mol_label(
        self,
        images: List[PIL.Image.Image],
        labels: List[str],
    ) -> Image:
        """Function to add legends to a list of images. The legends are added using PIL
        and placed according to the `label_loc` parameter.

        Args:
            images: list of PIL images.
            labels: list of labels to be added to the images.
            font_name: name of the font to be used. Depend on being found by the
                FontManager class. Defaults to "DejaVu Sans".
            font_size: size of the font. Defaults to 15.
            label_loc: position of the label. Defaults to "bottom".

        Raises:
            ValueError: if the font or the label_loc are not available

        Returns:
            list of the images with the labels added.
        """
        try:
            importlib.import_module("matplotlib")
            include_plt = True
        except ModuleNotFoundError:
            include_plt = False
        fm = FontManager(include_plt=include_plt)

        if self._font_name is None:
            font = ImageFont.load_default()
        elif self._font_name not in fm.available_fonts.keys():
            raise ValueError(
                "Font not available. Try one of the following: \n"
                f"{fm.available_fonts.keys()}"
            )
        else:
            font_path = fm.available_fonts[self._font_name]
            font = ImageFont.truetype(str(font_path), self._font_size)

        img_width, img_height = images[0].size
        # Writing the labels to each of the images
        for img, text in zip(images, labels):
            # Getting the size of the text to center the loation of the text
            font_width, font_height = font.getsize(text)
            adjusted_fontsize = self._font_size
            if font_width > img_width:
                # TODO: improve this so we can instead separate the text in two lines
                while True:
                    adjusted_fontsize -= 1
                    adjusted_font = ImageFont.truetype(
                        str(font_path), adjusted_fontsize
                    )
                    font_width, font_height = font.getsize(text)
                    if font_width < img_width:
                        break
            else:
                font = ImageFont.truetype(str(font_path), self._font_size)

            if self._label_loc not in ["top", "bottom"]:
                raise ValueError("label_loc must be either 'top' or 'bottom'.")
            if self._label_loc == "top":
                centered_w = (img_width - font_width) / 2
                centered_h = (img_height - font_height) / 99
            elif self._label_loc == "bottom":
                centered_w = (img_width - font_width) / 2
                centered_h = (img_height - font_height) / 99 * 97
            draw = ImageDraw.Draw(img)
            if adjusted_fontsize != self._font_size:
                draw.text(
                    (centered_w, centered_h), text, fill="black", font=adjusted_font
                )
            else:
                draw.text((centered_w, centered_h), text, fill="black", font=font)
        return images

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

    def _images_to_grid(self, images: List[Image.Image], n_cols: int) -> Image:
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

    def render_mol(
        self,
        mol: Chem.Mol,
        label: Optional[str] = None,
        match_struct: Optional[Union[str, Chem.Mol]] = None,
        **kwargs,
    ) -> Image:
        """Plot molecules using rdMolDraw2D.MolDraw2DCairo. Keyword arguments provided
        will be passed into DrawMolecule method of the same class.

        Args:
            mol: rdkit mol object.

        Returns:
            a PIL.Image.Image object.
        """
        mol = self._stdin_to_mol(mol)
        if match_struct is not None:
            if isinstance(match_struct, Chem.rdchem.Mol):
                AllChem.Compute2DCoords(match_struct)
            else:
                mol = Chem.MolFromSmiles(match_struct)
                if mol is None:
                    logging.warning(
                        'Error: "match_struct" from SMILES RDKit Mol is invalid.'
                    )
                    logging.warning('Trying "match_struct" from SMARTS...')
                    mol = Chem.MolFromSmarts(match_struct)
                    if mol is None:
                        raise ValueError(
                            'Error: "match_struct" invalid from SMILES & SMARTS.'
                        )
                match_struct = mol
            AllChem.Compute2DCoords(match_struct)
            AllChem.GenerateDepictionMatching2DStructure(
                mol, match_struct, acceptFailure=True
            )
        drawer = rdMolDraw2D.MolDraw2DCairo(*self._size)
        drawer.DrawMolecule(mol, **kwargs)
        drawer.FinishDrawing()
        img = Image.open(BytesIO(drawer.GetDrawingText()))
        if label is not None:
            img = self._add_mol_label([img], [label])[0]
        return img

    def render_with_matches(
        self,
        mol: Chem.Mol,
        substructs,
        label: Optional[str] = None,
        scaff_pose: Optional[Union[str, Chem.Mol]] = None,
    ) -> Image:
        """Render the substructures on the molecules using default RDKit coloring.

        Args:
            mol: molecule to be rendered.
            substructs: substructures output from RDKitFilters.filter_mols.
            molSize: final image size. Defaults to (500, 500).

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
            match_struct=scaff_pose,
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
        scaff_pose: Optional[Union[str, Chem.Mol]] = None,
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
            ax: if desired, add the results to an existing axis. Defaults to None.
            molSize: size of the molecule rendering. Defaults to (500, 500).
                Note: in the end this will be bound to the size of the plot...

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
            match_struct=scaff_pose,
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
        # ax.imshow(img)
        ax.set_axis_off()
        if ax is None:
            return fig, ax


class MolGridPlotter(MolPlotter):
    """MolGridPlotter class to render molecules into grids of images. It uses RDKit's
    DrawMolecule method to generate the images and PIL to patch images together.
    The images can be rendered with labels and/or substructure matches, similar to its
    parent class.

    Args:
        from_smi: output images directly from smiles. Defaults to False.
        size: size of the output images. Defaults to (500, 500).
        cmap: colormap for `render_with_color` method. Defaults to "rainbow".
        font_name: font name found by img_render.FondManager class.
            Defaults to "DejaVu Sans".
        font_size: size of the font patched to the image. Defaults to 15.
        label_loc: font placement; either `top` or `bottom`. Defaults to "bottom".
        n_jobs: number of jobs to run in parallel. Defaults to 1.
    """

    def __init__(self, n_jobs: int = 1, from_smi=False, **kwargs):
        super().__init__(from_smi=from_smi, n_jobs=n_jobs, **kwargs)

    def mol_grid_png(
        self,
        mols: list,
        labels: Optional[list] = None,
        n_cols: int = None,
        scaff_pose: Optional[Union[str, Chem.Mol]] = None,
        bw: bool = False,
    ) -> PIL.Image.Image:
        """
        Args:
            mols: list of molecules to display.
            labels: list of labels to be written on the mol images. Defaults to None
            n_cols: N columns with molecules. If None, n_cols = len(smiles).
                Defaults to None.
            scaff_pose: a common structure between `smiles` used to generate the
                matching 2D images. If None, 2D rendering won't match. Defaults to None.
            bw: if True, final image will be converted to gray scale. Defaults to False.
            label_loc: Location of the label. Defaults to "bottom".

        Returns:
            PIL.Image.Image
        """
        with Pool(self._n_jobs) as p:
            images = p.map(partial(self.render_mol, match_struct=scaff_pose), mols)
        if labels is not None:
            # Putting the fonts on the images
            images = self._add_mol_label(images, labels)
        image = self._images_to_grid(images, n_cols)
        if bw:
            image = image.convert("L")
        return image

    def mol_structmatch_grid_png(
        self,
        mols: list,
        substructs: list,
        labels: Optional[list] = None,
        n_cols: int = None,
        scaff_pose: Optional[Union[str, Chem.Mol]] = None,
        bw: bool = False,
    ) -> PIL.Image.Image:
        """
        Args:
            mols: list of molecules to display.
            labels: list of labels to be written on the mol images. Defaults to None
            n_cols: N columns with molecules. If None, n_cols = len(smiles).
                Defaults to None.
            scaff_pose: a common structure between `smiles` used to generate the
                matching 2D images. If None, 2D rendering won't match. Defaults to None.
            bw: if True, final image will be converted to gray scale. Defaults to False.
            label_loc: Location of the label. Defaults to "bottom".

        Returns:
            PIL.Image.Image
        """
        partial_function = partial(self.render_with_matches, scaff_pose=scaff_pose)
        if labels is not None:
            assert len(labels) == len(mols), "Labels and mols must have the same length"
            variables = list(zip(mols, substructs, labels))
        else:
            variables = list(zip(mols, substructs))
        with Pool(self._n_jobs) as p:
            images = p.starmap(partial_function, variables)
        if labels is not None:
            # Putting the fonts on the images
            images = self._add_mol_label(images, labels)
        image = self._images_to_grid(images, n_cols)
        if bw:
            image = image.convert("L")
        return image

    def mol_structmatch_color_grid_png(
        self,
        mols: list,
        descriptions: List[str],
        substructs: list,
        labels: Optional[list] = None,
        n_cols: int = None,
        scaff_pose: Optional[Union[str, Chem.Mol]] = None,
        bw: bool = False,
    ) -> PIL.Image.Image:
        """
        Args:
            mols: list of molecules to display.
            labels: list of labels to be written on the mol images. Defaults to None
            n_cols: N columns with molecules. If None, n_cols = len(smiles).
                Defaults to None.
            scaff_pose: a common structure between `smiles` used to generate the
                matching 2D images. If None, 2D rendering won't match. Defaults to None.
            bw: if True, final image will be converted to gray scale. Defaults to False.
            label_loc: Location of the label. Defaults to "bottom".

        Returns:
            PIL.Image.Image
        """
        mols = self._stdin_to_mol(mols)
        partial_function = partial(
            self.render_with_colored_matches, scaff_pose=scaff_pose
        )
        if labels is not None:
            assert len(labels) == len(mols), "Labels and mols must have the same length"
            variables = list(zip(mols, descriptions, substructs, labels))
        else:
            variables = list(zip(mols, descriptions, substructs))
        with Pool(self._n_jobs) as p:
            images = p.starmap(partial_function, variables)
        if labels is not None:
            # Putting the fonts on the images
            images = self._add_mol_label(images, labels)
        image = self._images_to_grid(images, n_cols)
        if bw:
            image = image.convert("L")
        return image
