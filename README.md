[![Linting Ruff](https://img.shields.io/badge/Linting%20-Ruff-red?style=flat-square)](https://github.com/charliermarsh/ruff)
[![Code style: black](https://img.shields.io/badge/code%20style-black-black?style=flat-square)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat-square&labelColor=ef8336)](https://pycqa.github.io/isort/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)](https://opensource.org/licenses/MIT)

# chemFilters

A collection of chemical filters, with some support for data visualization and analysis. At the moment, the supported filters are:

- RDKit filters based on RDKit implementation of [Pat Walter's version of the ChEMBL filters](https://github.com/PatWalters/rd_filters)*;
- Purchasability filters based on Andrew White's [molbloom](https://github.com/whitead/molbloom);
- Peptide filters based on Olivier Béquignon's [PepSift](https://github.com/OlivierBeq/PepSift);
- Silly molecules filters based Olivier's fork of Pat Water's [silly walks](https://github.com/PatWalters/silly_walks);

*note: RDKit's implementation these chemical filters is only available from rdkit version 2023.03.1 onwards. Check here for the [release notes](https://greglandrum.github.io/rdkit-blog/posts/2023-05-03-new-release-highlights.html).

## Overview:


## Installation

```bash
TODO: make the installation
```

## Usage
    
```python
TODO: make the usage
```

Project to do's:
- [x] Implement bloom filters;
- [x] Implement molspotter from Olivier;
- [ ] Integrate filter rendering according to moleculesToPng;
- [ ] Add support for interactive visualization using streamlit;
- [x] make stuff `snake_case`;
- [ ] Update README.md;