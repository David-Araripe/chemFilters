# ABOUTME: Generates all visualization figures for the chemFilters documentation.
# ABOUTME: Produces PNGs in docs/_static/figures/visualization/ for embedding in Sphinx docs.

from pathlib import Path

import numpy as np
from PIL import Image, ImageOps
from rdkit import Chem

from chemFilters import RdkitFilters
from chemFilters.img_render import MolGridPlotter, MolPlotter

OUTPUT_DIR = Path(__file__).parent / "_static" / "figures" / "visualization"

# Reference molecules (same as README examples)
SMILES = [
    "CCC1=[O+][Cu-3]2([O+]=C(CC)C1)[O+]=C(CC)CC(CC)=[O+]2",
    "CC1=C2C(=COC(C)C2C)C(O)=C(C(=O)O)C1=O",
    "CCOP(=O)(Nc1cccc(Cl)c1)OCC",
    "Nc1ccc(C=Cc2ccc(N)cc2S(=O)(=O)O)c(S(=O)(=O)O)c1",
]
LABELS = [f"Molecule {i}" for i in range(1, len(SMILES) + 1)]
MOLS = [Chem.MolFromSmiles(smi) for smi in SMILES]

ASPIRIN_SMI = "CC(=O)Oc1ccccc1C(=O)O"
ASPIRIN = Chem.MolFromSmiles(ASPIRIN_SMI)

FONT = "Telex-Regular"
MOL_SIZE = (400, 400)
GRID_SIZE = (250, 250)


def save(img: Image.Image, name: str) -> None:
    """Save a PIL image to the output directory."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    img.save(OUTPUT_DIR / name)
    print(f"  Saved {name}")


def make_comparison(
    img_left: Image.Image,
    img_right: Image.Image,
    label_left: str = "Default",
    label_right: str = "Modified",
) -> Image.Image:
    """Place two images side by side with labels underneath."""
    from PIL import ImageDraw, ImageFont

    gap = 20
    label_height = 30
    w = img_left.width + gap + img_right.width
    h = max(img_left.height, img_right.height) + label_height
    canvas = Image.new("RGB", (w, h), "white")
    canvas.paste(img_left, (0, 0))
    canvas.paste(img_right, (img_left.width + gap, 0))

    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("Arial", 16)
    except OSError:
        font = ImageFont.load_default()

    img_top = max(img_left.height, img_right.height)
    left_center = img_left.width // 2
    right_center = img_left.width + gap + img_right.width // 2
    draw.text(
        (left_center, img_top + 5), label_left, fill="black", font=font, anchor="mt"
    )
    draw.text(
        (right_center, img_top + 5),
        label_right,
        fill="black",
        font=font,
        anchor="mt",
    )
    return canvas


def make_triple(
    imgs: list,
    labels: list,
) -> Image.Image:
    """Place three images side by side with labels underneath."""
    from PIL import ImageDraw, ImageFont

    gap = 20
    label_height = 30
    w = sum(img.width for img in imgs) + gap * (len(imgs) - 1)
    h = max(img.height for img in imgs) + label_height
    canvas = Image.new("RGB", (w, h), "white")

    x_offset = 0
    centers = []
    for img in imgs:
        canvas.paste(img, (x_offset, 0))
        centers.append(x_offset + img.width // 2)
        x_offset += img.width + gap

    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("Arial", 16)
    except OSError:
        font = ImageFont.load_default()

    img_top = max(img.height for img in imgs)
    for center, label in zip(centers, labels):
        draw.text((center, img_top + 5), label, fill="black", font=font, anchor="mt")
    return canvas


# ---------------------------------------------------------------------------
# Figure generators
# ---------------------------------------------------------------------------


def generate_single_mol():
    """Single aspirin molecule."""
    plotter = MolPlotter(from_smi=False, size=MOL_SIZE, font_name=FONT)
    img = plotter.render_mol(ASPIRIN, label="Aspirin")
    save(img, "single_mol.png")


def generate_acs1996():
    """ACS 1996 style rendering with flexible canvas as SVG."""
    plotter = MolPlotter(
        from_smi=False, size=(-1, -1), font_name=FONT, label_font_size=8
    )
    svg = plotter.render_ACS1996(ASPIRIN, label="Aspirin", return_svg=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "acs1996.svg").write_text(svg)
    print("  Saved acs1996.svg")


def generate_flexible_canvas():
    """Fixed size vs flexible canvas side by side."""
    plotter_fixed = MolPlotter(from_smi=False, size=MOL_SIZE, font_name=FONT)
    plotter_flex = MolPlotter(from_smi=False, size=(-1, -1), font_name=FONT)
    img_fixed = plotter_fixed.render_mol(ASPIRIN, label="Aspirin")
    img_flex = plotter_flex.render_mol(ASPIRIN, label="Aspirin")
    save(
        make_comparison(img_fixed, img_flex, "size=(400, 400)", "size=(-1, -1)"),
        "flexible_canvas.png",
    )


def generate_simple_grid():
    """2x2 molecule grid with labels."""
    grid_plotter = MolGridPlotter(from_smi=False, font_name=FONT, size=GRID_SIZE)
    img = grid_plotter.mol_grid_png(MOLS, n_cols=2, labels=LABELS)
    save(img, "simple_grid.png")


def generate_substruct_grid():
    """Grid with highlighted substructure matches."""
    chemFilter = RdkitFilters(filter_type="ALL")
    _, _, substructs = chemFilter.filter_mols(MOLS)

    grid_plotter = MolGridPlotter(from_smi=False, font_name=FONT, size=GRID_SIZE)
    img = grid_plotter.mol_structmatch_grid_png(MOLS, substructs=substructs, n_cols=2)
    save(img, "substruct_grid.png")


def generate_colored_matches():
    """Colored substructure match highlighting."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chemFilter = RdkitFilters(filter_type="NIH")
    _, descriptions, substructs = chemFilter.filter_mols(MOLS)

    plotter = MolPlotter(
        from_smi=False, label_font_size=20, size=(350, 350), font_name=FONT
    )
    img = plotter.render_with_colored_matches(
        MOLS[0],
        descriptions=descriptions[0],
        substructs=substructs[0],
        label="Molecule 1",
        alpha=0.3,
    )

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(img)
    ax.set_axis_off()
    plotter.colored_matches_legend(descriptions[0], substructs[0], ax=ax)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "colored_matches.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved colored_matches.png")


def generate_bw_comparison():
    """Default vs black-and-white rendering."""
    plotter_default = MolPlotter(from_smi=False, size=MOL_SIZE, font_name=FONT)
    plotter_bw = MolPlotter(from_smi=False, size=MOL_SIZE, font_name=FONT, bw=True)
    img_default = plotter_default.render_mol(ASPIRIN, label="Aspirin")
    img_bw = plotter_bw.render_mol(ASPIRIN, label="Aspirin")
    save(make_comparison(img_default, img_bw, "Default", "bw=True"), "bw_comparison.png")


def generate_atom_indices():
    """Default vs atom indices shown."""
    plotter_default = MolPlotter(from_smi=False, size=MOL_SIZE, font_name=FONT)
    plotter_idx = MolPlotter(
        from_smi=False, size=MOL_SIZE, font_name=FONT, add_atom_indices=True
    )
    img_default = plotter_default.render_mol(ASPIRIN, label="Aspirin")
    img_idx = plotter_idx.render_mol(ASPIRIN, label="Aspirin")
    save(
        make_comparison(img_default, img_idx, "Default", "add_atom_indices=True"),
        "atom_indices.png",
    )


def generate_bond_indices():
    """Default vs bond indices shown."""
    plotter_default = MolPlotter(from_smi=False, size=MOL_SIZE, font_name=FONT)
    plotter_idx = MolPlotter(
        from_smi=False, size=MOL_SIZE, font_name=FONT, add_bond_indices=True
    )
    img_default = plotter_default.render_mol(ASPIRIN, label="Aspirin")
    img_idx = plotter_idx.render_mol(ASPIRIN, label="Aspirin")
    save(
        make_comparison(img_default, img_idx, "Default", "add_bond_indices=True"),
        "bond_indices.png",
    )


def generate_explicit_methyl():
    """Default vs explicit methyl groups."""
    plotter_default = MolPlotter(from_smi=False, size=MOL_SIZE, font_name=FONT)
    plotter_ch3 = MolPlotter(
        from_smi=False, size=MOL_SIZE, font_name=FONT, explicit_methyl=True
    )
    img_default = plotter_default.render_mol(ASPIRIN, label="Aspirin")
    img_ch3 = plotter_ch3.render_mol(ASPIRIN, label="Aspirin")
    save(
        make_comparison(img_default, img_ch3, "Default", "explicit_methyl=True"),
        "explicit_methyl.png",
    )


def generate_no_atom_labels():
    """Default vs no atom labels."""
    plotter_default = MolPlotter(from_smi=False, size=MOL_SIZE, font_name=FONT)
    plotter_no = MolPlotter(
        from_smi=False, size=MOL_SIZE, font_name=FONT, no_atom_labels=True
    )
    img_default = plotter_default.render_mol(ASPIRIN, label="Aspirin")
    img_no = plotter_no.render_mol(ASPIRIN, label="Aspirin")
    save(
        make_comparison(img_default, img_no, "Default", "no_atom_labels=True"),
        "no_atom_labels.png",
    )


def generate_bond_line_width():
    """Three bond line widths compared."""
    widths = [1.0, 2.0, 4.0]
    imgs = []
    for w in widths:
        plotter = MolPlotter(
            from_smi=False, size=MOL_SIZE, font_name=FONT, bond_line_width=w
        )
        imgs.append(plotter.render_mol(ASPIRIN, label="Aspirin"))
    save(
        make_triple(imgs, [f"bond_line_width={w}" for w in widths]),
        "bond_line_width.png",
    )


def generate_transparent_bg():
    """Render with transparent background."""
    plotter = MolPlotter(
        from_smi=False, size=MOL_SIZE, font_name=FONT, bg_transparent=True
    )
    img = plotter.render_mol(ASPIRIN, label="Aspirin")
    # Ensure RGBA mode for transparency
    if img.mode != "RGBA":
        img = img.convert("RGBA")
    save(img, "transparent_bg.png")


def generate_cmap_comparison():
    """Same molecule with three different colormaps."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chemFilter = RdkitFilters(filter_type="NIH")
    _, descriptions, substructs = chemFilter.filter_mols(MOLS)

    cmaps = ["rainbow", "Set1", "viridis"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for ax, cmap_name in zip(axes, cmaps):
        plotter = MolPlotter(
            from_smi=False, label_font_size=20, size=(350, 350), font_name=FONT
        )
        img = plotter.render_with_colored_matches(
            MOLS[0],
            descriptions=descriptions[0],
            substructs=substructs[0],
            label="Molecule 1",
            cmap=cmap_name,
            alpha=0.3,
        )
        ax.imshow(img)
        ax.set_axis_off()
        ax.set_title(f'cmap="{cmap_name}"', fontsize=14)
        plotter.colored_matches_legend(
            descriptions[0], substructs[0], cmap=cmap_name, ax=ax
        )

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "cmap_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved cmap_comparison.png")


def generate_chess_grid():
    """2x2 chess-style grid with alternating black/white cells."""
    grid_plotter = MolGridPlotter(
        from_smi=False, font_name=FONT, size=GRID_SIZE, bw=True
    )

    # Render each molecule individually
    cell_imgs = []
    for i, mol in enumerate(MOLS):
        img = grid_plotter.render_mol(mol, label=LABELS[i])
        cell_imgs.append(img)

    # Invert cells where (row + col) is even for chess pattern
    processed = []
    for idx, img in enumerate(cell_imgs):
        row, col = divmod(idx, 2)
        if (row + col) % 2 == 0:
            processed.append(ImageOps.invert(img.convert("RGB")))
        else:
            processed.append(img.convert("RGB"))

    # Stitch into 2x2 grid
    w, h = processed[0].size
    top_row = np.hstack([np.array(processed[0]), np.array(processed[1])])
    bottom_row = np.hstack([np.array(processed[2]), np.array(processed[3])])
    grid = np.vstack([top_row, bottom_row])
    save(Image.fromarray(grid), "chess_grid.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

ALL_GENERATORS = [
    ("Single molecule", generate_single_mol),
    ("ACS 1996", generate_acs1996),
    ("Flexible canvas", generate_flexible_canvas),
    ("Simple grid", generate_simple_grid),
    ("Substructure grid", generate_substruct_grid),
    ("Colored matches", generate_colored_matches),
    ("B&W comparison", generate_bw_comparison),
    ("Atom indices", generate_atom_indices),
    ("Bond indices", generate_bond_indices),
    ("Explicit methyl", generate_explicit_methyl),
    ("No atom labels", generate_no_atom_labels),
    ("Bond line width", generate_bond_line_width),
    ("Transparent background", generate_transparent_bg),
    ("Colormap comparison", generate_cmap_comparison),
    ("Chess grid", generate_chess_grid),
]

if __name__ == "__main__":
    print(f"Output directory: {OUTPUT_DIR}")
    for name, func in ALL_GENERATORS:
        print(f"Generating: {name}")
        func()
    print("Done!")
