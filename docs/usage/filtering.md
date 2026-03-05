# Filtering compounds

chemFilters provides several filter classes, each wrapping a different filtering tool.
All filters follow a similar API: initialize the filter, then call `get_flagging_df` (or
`get_scoring_df`) on a list of molecules to get a DataFrame of results.

## RDKit structural alert filters

`RdkitFilters` is always available with the base install. It wraps RDKit's
`FilterCatalog` system.

```python
from chemFilters import RdkitFilters
from rdkit import Chem

mols = [
    Chem.MolFromSmiles("CCC1=[O+][Cu-3]2([O+]=C(CC)C1)[O+]=C(CC)CC(CC)=[O+]2"),
    Chem.MolFromSmiles("CC1=C2C(=COC(C)C2C)C(O)=C(C(=O)O)C1=O"),
    Chem.MolFromSmiles("CCOP(=O)(Nc1cccc(Cl)c1)OCC"),
    Chem.MolFromSmiles("Nc1ccc(C=Cc2ccc(N)cc2S(=O)(=O)O)c(S(=O)(=O)O)c1"),
]

rdkit_filter = RdkitFilters(filter_type="ALL", from_smi=False)
flagging_df = rdkit_filter.get_flagging_df(mols)
```

The `filter_type` parameter selects which filter catalog to use. See
`RdkitFilters.available_filters` for the full list. Common choices include `"ALL"`,
`"PAINS"`, `"CHEMBL"`, and `"BRENK"`.

You can also retrieve substructure matches directly:

```python
filter_names, descriptions, substructs = rdkit_filter.filter_mols(mols)
```

## Purchasability filters (molbloom)

:::{note}
Requires the `allfilters` extra: `pip install 'chem-filters[allfilters]'`
:::

```python
from chemFilters import MolbloomFilters

bloom_filter = MolbloomFilters(from_smi=False, standardize=False)
bloom_filter.get_flagging_df(mols)
```

Results indicate whether a molecule is *probably* in a given catalog (`True`) or
*definitely not* (`False`).

## Peptide filters (PepSift)

:::{note}
Requires the `allfilters` extra: `pip install 'chem-filters[allfilters]'`
:::

```python
from chemFilters import PeptideFilters

pep_filter = PeptideFilters(from_smi=False)
pep_filter.get_flagging_df(mols)
```

## Silly molecule filters (molspotter)

:::{note}
Requires the `allfilters` extra: `pip install 'chem-filters[allfilters]'`
:::

```python
from chemFilters import SillyMolSpotterFilter

silly_filter = SillyMolSpotterFilter(from_smi=False)
silly_filter.get_scoring_df(mols)
```

Scores indicate how "unusual" a molecule is based on detection of rare bits in hashed
ECFP fingerprints.

## Running all filters at once (CoreFilter)

`CoreFilter` combines all available filters into a single callable. This is also what
the CLI uses under the hood.

```python
from chemFilters.core import CoreFilter

smiles = [
    "CCC1=[O+][Cu-3]2([O+]=C(CC)C1)[O+]=C(CC)CC(CC)=[O+]2",
    "CC1=C2C(=COC(C)C2C)C(O)=C(C(=O)O)C1=O",
    "CCOP(=O)(Nc1cccc(Cl)c1)OCC",
    "Nc1ccc(C=Cc2ccc(N)cc2S(=O)(=O)O)c(S(=O)(=O)O)c1",
]

core_filter = CoreFilter()
filtered_df = core_filter(smiles)
```

Individual filters can be toggled on/off via constructor arguments (e.g.,
`pep_filter=False`). Data can be processed in chunks with the `chunksize` parameter.
