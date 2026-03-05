# Visualization

chemFilters includes `MolPlotter` and `MolGridPlotter` for rendering molecules with
RDKit's `MolDraw2DCairo`. Both classes support extensive customization of the rendering
process through constructor parameters and RDKit draw options.

## Single molecule rendering

`MolPlotter` renders individual molecules as PIL images or SVG strings.

```python
from rdkit import Chem
from chemFilters.img_render import MolPlotter

mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")

plotter = MolPlotter(from_smi=False, size=(400, 400))
img = plotter.render_mol(mol, label="Aspirin")
```

```{image} /_static/figures/visualization/single_mol.png
:alt: Single molecule rendering of aspirin
:width: 400px
```

### SVG output

Pass `return_svg=True` to get an SVG string instead of a PIL image:

```python
svg_text = plotter.render_mol(mol, label="Aspirin", return_svg=True)
```

### ACS 1996 style

Render molecules following the American Chemical Society (ACS) 1996 drawing standards.
Because ACS 1996 enforces fixed bond lengths, the molecule size on the canvas is
determined by the standard rather than the `size` parameter. Passing a large fixed size
will result in blank space around the structure. Using `size=(-1, -1)` for a flexible
canvas (see [Flexible canvas](#flexible-canvas)) avoids this, but produces a small
image since bond lengths are fixed. Use `label_font_size` to keep the label
proportional:

```python
plotter = MolPlotter(from_smi=False, size=(-1, -1), label_font_size=8)
img = plotter.render_ACS1996(mol, label="Aspirin")
```

```{image} /_static/figures/visualization/acs1996.svg
:alt: Aspirin rendered in ACS 1996 style
:width: 300px
```

### Pose matching

Align a molecule's 2D depiction to a reference structure using `match_pose`. This
accepts SMILES, SMARTS, or an `rdkit.Chem.Mol` object:

```python
img = plotter.render_mol(mol, match_pose="c1ccccc1")
```

## Molecule grids

`MolGridPlotter` extends `MolPlotter` to arrange multiple molecules into a grid image.

```python
from chemFilters.img_render import MolGridPlotter

mols = [
    Chem.MolFromSmiles("CCC1=[O+][Cu-3]2([O+]=C(CC)C1)[O+]=C(CC)CC(CC)=[O+]2"),
    Chem.MolFromSmiles("CC1=C2C(=COC(C)C2C)C(O)=C(C(=O)O)C1=O"),
    Chem.MolFromSmiles("CCOP(=O)(Nc1cccc(Cl)c1)OCC"),
    Chem.MolFromSmiles("Nc1ccc(C=Cc2ccc(N)cc2S(=O)(=O)O)c(S(=O)(=O)O)c1"),
]
labels = [f"Molecule {i}" for i in range(1, len(mols) + 1)]

grid_plotter = MolGridPlotter(from_smi=False, font_name="Telex-Regular")
img = grid_plotter.mol_grid_png(mols, n_cols=2, labels=labels)
```

```{image} /_static/figures/visualization/simple_grid.png
:alt: 2x2 molecule grid with labels
:width: 500px
```

## Substructure match highlighting

After filtering molecules with `RdkitFilters`, you can highlight the matched
substructures on a grid:

```python
from chemFilters import RdkitFilters

chemFilter = RdkitFilters(filter_type="ALL")
filter_names, descriptions, substructs = chemFilter.filter_mols(mols)

grid_plotter = MolGridPlotter(
    from_smi=False, font_name="Telex-Regular", size=(250, 250)
)
img = grid_plotter.mol_structmatch_grid_png(mols, substructs=substructs, n_cols=2)
```

```{image} /_static/figures/visualization/substruct_grid.png
:alt: Molecule grid with highlighted substructure matches
:width: 500px
```

For a single molecule, use `render_with_matches` directly:

```python
plotter = MolPlotter(from_smi=False)
img = plotter.render_with_matches(
    mols[0], substructs=substructs[0], label="Molecule 1"
)
```

## Colored substructure matches

Each matched substructure can be rendered in a different color. When substructures
overlap, the colors are blended using their geometric mean.

```python
import matplotlib.pyplot as plt
from chemFilters.img_render import MolPlotter

chemFilter = RdkitFilters(filter_type="NIH")
filter_names, descriptions, substructs = chemFilter.filter_mols(mols)

plotter = MolPlotter(
    from_smi=False, label_font_size=20, size=(350, 350), font_name="Telex-Regular"
)
img = plotter.render_with_colored_matches(
    mols[0],
    descriptions=descriptions[0],
    substructs=substructs[0],
    label="Molecule 1",
    alpha=0.3,
)

plt.imshow(img)
ax = plt.gca()
ax.set_axis_off()
plotter.colored_matches_legend(descriptions[0], substructs[0], ax=ax)
```

```{image} /_static/figures/visualization/colored_matches.png
:alt: Colored substructure match highlighting with legend
:width: 500px
```

For grids with colored matches, use `mol_structmatch_color_grid_png`:

```python
grid_plotter = MolGridPlotter(from_smi=False, size=(350, 350))
img = grid_plotter.mol_structmatch_color_grid_png(
    mols, descriptions=descriptions, substructs=substructs, n_cols=2
)
```

### Colormap options

The `cmap` parameter controls the colormap used for highlighting substructures.
Any matplotlib colormap name is accepted:

```{image} /_static/figures/visualization/cmap_comparison.png
:alt: Comparison of rainbow, Set1, and viridis colormaps
:width: 700px
```

## Customizing rendering options

Both `MolPlotter` and `MolGridPlotter` expose several options that map to RDKit's
`MolDraw2DCairo` settings:

| Parameter                    | Default        | Description                                     |
|------------------------------|----------------|-------------------------------------------------|
| `size`                       | `(300, 300)`   | Image dimensions in pixels; `(-1, -1)` for flexible canvas |
| `cmap`                       | `"rainbow"`    | Matplotlib colormap for colored matches         |
| `font_name`                  | `"Telex-Regular"` | Font for molecule labels                    |
| `label_font_size`            | `None`         | Custom font size for labels                     |
| `mol_font_size`              | `None`         | Custom font size on the molecule drawing        |
| `annotation_font_scale`      | `None`         | Scale factor for annotation font size           |
| `bond_line_width`            | `2.0`          | Width of bond lines                             |
| `no_atom_labels`             | `False`        | Suppress atom labels                            |
| `add_atom_indices`           | `False`        | Display atom indices                            |
| `add_bond_indices`           | `False`        | Display bond indices                            |
| `explicit_methyl`            | `False`        | Show methyl groups as CH3                       |
| `unspecified_stereo_unknown` | `False`        | Draw unspecified stereo as unknown               |
| `bg_transparent`             | `False`        | Transparent background                          |
| `bw`                         | `False`        | Black-and-white rendering                       |

Additional `MolDraw2D` options can be passed as keyword arguments:

```python
plotter = MolPlotter(
    from_smi=False,
    size=(500, 500),
    bw=True,
    bg_transparent=True,
    bond_line_width=3.0,
)
```

### Visual parameter effects

(flexible-canvas)=
#### Flexible canvas

Pass `size=(-1, -1)` to let RDKit auto-size the image to fit the molecule, avoiding
blank space around the structure:

```python
plotter = MolPlotter(from_smi=False, size=(-1, -1))
img = plotter.render_mol(mol, label="Aspirin")
```

```{image} /_static/figures/visualization/flexible_canvas.png
:alt: Fixed size vs flexible canvas comparison
:width: 500px
```

#### Black and white mode

```{image} /_static/figures/visualization/bw_comparison.png
:alt: Default vs black and white rendering
:width: 700px
```

#### Atom and bond indices

```{image} /_static/figures/visualization/atom_indices.png
:alt: Default vs atom indices shown
:width: 700px
```

```{image} /_static/figures/visualization/bond_indices.png
:alt: Default vs bond indices shown
:width: 700px
```

#### Explicit methyl groups

```{image} /_static/figures/visualization/explicit_methyl.png
:alt: Default vs explicit methyl groups
:width: 700px
```

#### Suppressing atom labels

```{image} /_static/figures/visualization/no_atom_labels.png
:alt: Default vs no atom labels
:width: 700px
```

#### Bond line width

```{image} /_static/figures/visualization/bond_line_width.png
:alt: Bond line width 1.0, 2.0, and 4.0 comparison
:width: 700px
```

#### Transparent background

```{raw} html
<div class="bg-checkerboard">
```

```{image} /_static/figures/visualization/transparent_bg.png
:alt: Molecule rendered with transparent background
:width: 400px
```

```{raw} html
</div>
```

### Chess grid example

If for some reason (?) you want to create a chess-like rendering of compounds, you can combine
the `bw=True` option with PIL `ImageOps.invert` to produce alternating black-on-white and 
white-on-black cells:

```python
from PIL import ImageOps

grid_plotter = MolGridPlotter(from_smi=False, font_name="Telex-Regular", size=(250, 250), bw=True)

cell_imgs = [grid_plotter.render_mol(mol, label=label) for mol, label in zip(mols, labels)]

processed = []
for idx, img in enumerate(cell_imgs):
    row, col = divmod(idx, 2)
    if (row + col) % 2 == 0:
        processed.append(ImageOps.invert(img.convert("RGB")))
    else:
        processed.append(img.convert("RGB"))
```

```{image} /_static/figures/visualization/chess_grid.png
:alt: Chess-style 2x2 grid with alternating black and white cells
:width: 500px
```

## Font management

`FontManager` discovers fonts available on the system, including fonts from matplotlib
and RDKit. It supports Windows, Linux, macOS, and WSL.

```python
from chemFilters.img_render import FontManager

fm = FontManager()

# List all available fonts
print(fm.list_available_fonts())

# Check if a specific font exists
fm.has_font("Telex-Regular")

# Get the path to a font file
fm.get_font_path("Telex-Regular")
```

## Parallel rendering

Both `MolPlotter` and `MolGridPlotter` support parallel rendering via `n_jobs`:

```python
grid_plotter = MolGridPlotter(from_smi=False, n_jobs=4)
img = grid_plotter.mol_grid_png(mols, n_cols=2)
```
