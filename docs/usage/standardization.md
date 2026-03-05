# Standardization

chemFilters provides molecular standardization through `ChemStandardizer` and identifier
conversion through `InchiHandling`.

## ChemStandardizer

`ChemStandardizer` supports several standardization pipelines:

| Method      | Description                        | Extra required      |
|-------------|------------------------------------|---------------------|
| `"chembl"`  | ChEMBL structure pipeline         | *(base install)*    |
| `"canon"`   | RDKit canonical SMILES             | *(base install)*    |
| `"papyrus"` | Papyrus structure pipeline         | `standardizers`     |
| `"molvs"`   | MolVS standardizer                 | `standardizers`     |
| *callable*  | Any function taking `rdkit.Mol`    | *(depends on func)* |

```python
from chemFilters.chem.standardizers import ChemStandardizer

standardizer = ChemStandardizer(method="chembl", from_smi=True)
standardized = standardizer(["CC(=O)Oc1ccccc1C(=O)O", "c1ccccc1"])
```

### Using papyrus or molvs

:::{note}
Requires the `standardizers` extra: `pip install 'chem-filters[standardizers]'`
:::

```python
standardizer = ChemStandardizer(method="papyrus", from_smi=True)
```

### Custom standardizer

You can pass any callable that accepts `rdkit.Chem.Mol` objects:

```python
def my_standardizer(mol):
    # your logic here
    return mol

standardizer = ChemStandardizer(method=my_standardizer, from_smi=False)
```

## InchiHandling

Convert molecules to InChI, InChIKey, or connectivity layers:

```python
from chemFilters.chem.standardizers import InchiHandling

converter = InchiHandling(convert_to="inchikey", from_smi=True)
keys = converter(["CC(=O)Oc1ccccc1C(=O)O", "c1ccccc1"])
```

The `convert_to` parameter accepts `"inchi"`, `"inchikey"`, or `"connectivity"`.
