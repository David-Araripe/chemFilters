# Installation

## Basic install

The base package includes RDKit structural alert filters, molecular standardization via
ChEMBL's pipeline, and visualization utilities:

```bash
pip install git+https://github.com/David-Araripe/chemFilters.git
```

## Optional extras

chemFilters uses optional dependency groups so you only install what you need:

| Extra            | What it adds                                      | Install command                                           |
|------------------|---------------------------------------------------|-----------------------------------------------------------|
| `allfilters`     | molbloom, pepsift, molspotter                     | `pip install 'chem-filters[allfilters]'`                  |
| `standardizers`  | papyrus_structure_pipeline, molvs                 | `pip install 'chem-filters[standardizers]'`               |
| `full`           | All of the above                                  | `pip install 'chem-filters[full]'`                        |
| `dev`            | full + pytest, ruff, isort, black                 | `pip install 'chem-filters[dev]'`                         |
| `docs`           | sphinx, furo, myst-parser, sphinx-autodoc-typehints | `pip install 'chem-filters[docs]'`                      |

## From source (development)

```bash
git clone https://github.com/David-Araripe/chemFilters.git
cd chemFilters
pip install -e '.[dev]'
```
