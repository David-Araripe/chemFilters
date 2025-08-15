<div align="center">

  <img src="logo.svg" alt="" width=360>
  <p><strong>Flag issues, standardize, and visualize molecular structures with ease.</strong></p>

[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Workflow](https://github.com/David-Araripe/chemFilters/actions/workflows/ci.yml/badge.svg?event=push)](https://github.com/David-Araripe/UniProtMapper/actions)

</div>

A collection of chemical filters, with some support for data visualization and analysis. Supported filters include:

- RDKit's [structural alert filters](https://www.rdkit.org/docs/source/rdkit.Chem.rdfiltercatalog.html#rdkit.Chem.rdfiltercatalog.FilterCatalogParams.FilterCatalogs)* including BMS, Dundee, Glaxo, Inpharmatica, LINT, MLSMR, PAINS, and SureChEMBL `FilterCatalogs`;
- Purchasability filters based on [molbloom](https://github.com/whitead/molbloom);
- SMARTS-like Peptide filters as implemented in [PepSift](https://github.com/OlivierBeq/PepSift);
- Silly molecules filters as implemented in [molspotter](https://github.com/OlivierBeq/molspotter);

*Note: RDKit's implementation these chemical filters is only available from rdkit version 2023.03.1 onwards. Check here for the [release notes](https://greglandrum.github.io/rdkit-blog/posts/2023-05-03-new-release-highlights.html).

## Overview:

The different filtering classes are implemented with a similar API, where `get_flagging/scoring_df` run all the filters available for that filtering class and return a dataframe that used to investigate the filters. In case of the RdkitFilters implementation, a few visualization methods are available to render the molecules, substructure matches, and molecular grids.

See available filters and visualization methods below:

- [chemFilters](#chemfilters)
  - [Overview:](#overview)
  - [Installation](#installation)
  - [Filtering Compounds](#filtering-datasets)
    - [RdkitFilters](#rdkitfilters)
    - [Purchasability filters](#purchasability-filters)
    - [Silly molecules filters](#silly-molecules-filters)
    - [Peptide filters](#peptide-filters)
    - [Core filters](#core-filters)
    - [CLI](#cli)
  - [Visualization](#visualization)
    - [Rendering a grid of molecules;](#rendering-a-grid-of-molecules)
    - [Rendering substructure matches:](#rendering-substructure-matches)
    - [Rendering substructure matches with colors:](#rendering-substructure-matches-with-colors)

## Installation

```bash
python -m pip install git+https://github.com/David-Araripe/chemFilters.git
```

## Filtering Compounds

### RdkitFilters
``` Python
from chemFilters import RdkitFilters
from rdkit import Chem

mols = [
    Chem.MolFromSmiles("CCC1=[O+][Cu-3]2([O+]=C(CC)C1)[O+]=C(CC)CC(CC)=[O+]2"),
    Chem.MolFromSmiles("CC1=C2C(=COC(C)C2C)C(O)=C(C(=O)O)C1=O"),
    Chem.MolFromSmiles("CCOP(=O)(Nc1cccc(Cl)c1)OCC"),
    Chem.MolFromSmiles("Nc1ccc(C=Cc2ccc(N)cc2S(=O)(=O)O)c(S(=O)(=O)O)c1"),
]

rdkit_filter = RdkitFilters(filter_type='ALL', from_smi=False)
filtered_df = rdkit_filter.get_flagging_df(mols)
```

### Purchasability filters

``` Python
from chemFilters import MolbloomFilters
bloom_filter = MolbloomFilters(from_smi=False, standardize=False)
bloom_filter.get_flagging_df(mols)
```

### Silly molecules filters

``` Python
from chemFilters import SillyMolFilters
silly_filter = SillyMolFilters(from_smi=False)
silly_filter.get_scoring_df(mols)
```

### Peptide filters

``` Python
from chemFilters import PeptideFilters
pep_filter = PeptideFilters(from_smi=False)
pep_filter.get_flagging_df(mols)
```

### Core filters

The package also has an implementation that allows applying all available filters at once. This implementation is also used in the CLI version of the package. For further configuration options, check the CLI help.

``` Python
from chemFilters.core import CoreFilters

smiles = [
    "CCC1=[O+][Cu-3]2([O+]=C(CC)C1)[O+]=C(CC)CC(CC)=[O+]2",
    "CC1=C2C(=COC(C)C2C)C(O)=C(C(=O)O)C1=O",
    "CCOP(=O)(Nc1cccc(Cl)c1)OCC",
    "Nc1ccc(C=Cc2ccc(N)cc2S(=O)(=O)O)c(S(=O)(=O)O)c1",
]

core_filter = CoreFilters()
filtered_df = core_filter(smiles)
```

### CLI

After installing the package, the CLI can be used to filter datasets. The CLI has the following options:

``` bash
usage: chemFilters [-h] -i INPUT [-c COL_NAME] -o OUTPUT [--rdkit-filter] [--no-rdkit-filter]
                   [--rdkit-subset RDKIT_SUBSET] [--rdkit-valtype RDKIT_VALTYPE] [--pep-filter] [--no-pep-filter]
                   [--silly-filter] [--no-silly-filter] [--bloom-filter] [--no-bloom-filter] [--std-mols]
                   [--no-std-mols] [--std-method STD_METHOD] [--n-jobs N_JOBS] [--chunk-size CHUNK_SIZE]
```

Where `--<name>-filter` and `--no-<name>-filter` enables and disables the implemented filters. Same goes for the parameter `--std-mols`, that enables the molecular standardization according to `--std-method`.

## Visualization

### Rendering a grid of molecules;

``` Python
from rdkit import Chem
from chemFilters.img_render import MolPlotter, MolGridPlotter

mols = [
    Chem.MolFromSmiles("CCC1=[O+][Cu-3]2([O+]=C(CC)C1)[O+]=C(CC)CC(CC)=[O+]2"),
    Chem.MolFromSmiles("CC1=C2C(=COC(C)C2C)C(O)=C(C(=O)O)C1=O"),
    Chem.MolFromSmiles("CCOP(=O)(Nc1cccc(Cl)c1)OCC"),
    Chem.MolFromSmiles("Nc1ccc(C=Cc2ccc(N)cc2S(=O)(=O)O)c(S(=O)(=O)O)c1"),
]
labels = [f"Molecule {i}" for i in range(1, len(mols) + 1)]

# Initialize grid plotter instance
grid_plotter = MolGridPlotter(from_smi=False, font_name="Telex-Regular")

img = grid_plotter.mol_grid_png(mols[:4], n_cols=2, labels=labels)
display(img)
```
<!-- img.save("figures/simple_grid.png") -->

<p align="center">
  <img src="./figures/simple_grid.png" alt="drawing" width="450"/>
</p>

### Rendering substructure matches:

``` Python
chemFilter = RdkitFilters(filter_type="ALL")
filter_names, description, substructs = chemFilter.filter_mols(mols)

grid_plotter = MolGridPlotter(
    from_smi=False, font_name="Telex-Regular", size=(250, 250)
)

img = grid_plotter.mol_structmatch_grid_png(mols, substructs=substructs, n_cols=2)
display(img)
```
<!-- img.save("figures/substruct_grid.png")  # saving the figure -->

<p align="center">
  <img src="./figures/substruct_grid.png" alt="drawing" width="450"/>
</p>

### Rendering substructure matches with colors:

``` Python
from chemFilters import RdkitFilters
import matplotlib.pyplot as plt

chemFilter = RdkitFilters(filter_type="NIH")
filter_names, description, substructs = chemFilter.filter_mols(mols)

plotter = MolPlotter(
    from_smi=False, font_size=20, size=(350, 350), font_name="Telex-Regular"
)
img = plotter.render_with_colored_matches(
    mols[0],
    descriptions=description[0],
    substructs=substructs[0],
    label=labels[0],
    alpha=0.3,
)

plt.imshow(img)
ax = plt.gca()  # get current axis
ax.set_axis_off()
plotter.colored_matches_legend(description[0], substructs[0], ax=ax)
fig = plt.gcf()  # get current figure
fig.savefig(  # save matplotlib figure
    "figures/colored_matches.png", bbox_inches="tight", dpi=150, facecolor="white"
)
```
<p align="center">
  <img src="./figures/colored_matches.png" alt="drawing" width="450"/>
</p>
