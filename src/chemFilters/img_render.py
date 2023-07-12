# -*- coding: utf-8 -*-
"""Contains a class for managing custom fonts and functions to render molecules into figures."""

import importlib
import os
import re
import sys
from functools import partial
from io import BytesIO
from multiprocessing import Pool
from pathlib import Path
from typing import List, Optional, Union

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import PIL
from matplotlib import cm
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

from .utils import smi_to_img, smiles_img_matching_pose


class FontManager:
    def __init__(self, wsl: bool = "auto", include_plt: bool = True):
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
            default_font_dir = Path("mnt/c/Windows/Fonts")
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

            font_dict.update({font.name: font.fname for font in fm.FontManager.ttflist})

        return font_dict


def render_colored_matches(
    mol: Chem.Mol,
    descriptions: List[str],
    substructs: List[Chem.Mol],
    cmap: str = "rainbow",
    alpha: float = 0.5,
    ax=None,
    molSize=(500, 500),
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
    unique_descrip = [descriptions[i] for i in indices]

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
        patches.append(mpatches.Patch(facecolor=color, label=descr, edgecolor="black"))
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
        bbox_to_anchor=(1.05, 0.25),
        loc="lower left",
        borderaxespad=0,
        frameon=False,
    )
    ax.imshow(img)
    ax.set_axis_off()
    plt.close()
    return fig, ax


class MolPlotter:
    def __init__(
        self,
        from_smi: bool = True,
        size: tuple = (500, 500),
        cmap: str = "rainbow",
        font_name: str = "DejaVu Sans",
        font_size: int = 15,
        label_loc="bottom",
    ) -> None:
        self._drawer = rdMolDraw2D.MolDraw2DCairo(size)
        self._cmap = cmap
        self._size = size
        self._from_smi = from_smi
        self._font_size = font_size
        self._font_name = font_name
        self._label_loc = label_loc

    def _add_legends(
        self,
        images: List[PIL.Image],
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

    def plot_mol(self, mol: Chem.Mol) -> Image:
        self._drawer.DrawMolecule(mol)
        self._drawer.FinishDrawing()
        img = Image.open(BytesIO(self._drawer.GetDrawingText()))
        return img

    def plot_mol_with_matches(self, mol: Chem.Mol, substructs) -> Image:
        """Render the substructures on the molecules using default RDKit coloring.

        Args:
            mol: molecule to be rendered.
            substructs: substructures output from RDKitFilters.filter_mols.
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

        drawer = rdMolDraw2D.MolDraw2DCairo(*self._size)
        drawer.DrawMolecule(mol, highlightAtoms=hit_atoms, highlightBonds=hit_bonds)
        drawer.FinishDrawing()
        img = Image.open(BytesIO(drawer.GetDrawingText()))
        return img

    def molgrid_png(
        self,
        mols: list,
        labels: Optional[list] = None,
        n_cols: int = None,
        scaff_pose: Optional[Union[str, Chem.Mol]] = None,
        n_jobs: int = 1,
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
            font_name: Path to a font to be used by the function. Defaults to None.
            font_size: Font size. Defaults to 15.
            n_jobs: Generate mol images in parallel. Defaults to 1.
            bw: if True, final image will be converted to gray scale. Defaults to False.
            label_loc: Location of the label. Defaults to "bottom".

        Returns:
            PIL.Image.Image
        """

        def list_split(list_a, chunk_size):
            """
            Returns a generator with a determined chunk size over list_a.
            Used to organize the figures in the correct `n_cols` format.
            """
            for i in range(0, len(list_a), chunk_size):
                yield list_a[i : i + chunk_size]

        images = list()
        if scaff_pose:
            images = smiles_img_matching_pose(
                mols, match_struct=scaff_pose, size=self._size
            )
        else:
            with Pool(n_jobs) as pool:
                images = pool.map(partial(smi_to_img, size=self._size), mols)

        if labels is not None:
            # Putting the fonts on the images
            images = self._add_legends(
                self._font_name, self._font_size, self._label_loc
            )

        # Splitting the list of images into a list of lists with n_cols
        if n_cols is None:
            n_cols = len(mols)
        images = list(list_split(images, n_cols))
        # Appending blank figures so we have the correct vector shapes
        while len(images[-1]) < len(images[0]):
            images[-1].append(Image.new("RGB", self._size, color=(255, 255, 255)))

        list_of_hstacked = list()
        # Creating list of horizontally stacked arrays
        for sublist in images:
            list_of_hstacked.append(np.hstack([np.asarray(img) for img in sublist]))
        # Vertically stacking horizontal arrays
        for item in list_of_hstacked:
            final_img = np.vstack([hstack for hstack in list_of_hstacked])
        # Creating and returning image from array
        final_img = Image.fromarray(final_img)
        if bw:
            final_img = final_img.convert("L")
        return final_img
