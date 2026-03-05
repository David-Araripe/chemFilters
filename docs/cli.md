# Command-line interface

After installing chemFilters, the `chemfilters` command is available for filtering
datasets from the terminal.

## Usage

```text
chemfilters -i INPUT -o OUTPUT [OPTIONS]
```

## Required arguments

| Flag               | Description                                                  |
|--------------------|--------------------------------------------------------------|
| `-i`, `--input`    | Path to the input file containing SMILES strings             |
| `-o`, `--output`   | Path to the output CSV file                                  |

## Optional arguments

| Flag                        | Default    | Description                                                                                                   |
|-----------------------------|------------|---------------------------------------------------------------------------------------------------------------|
| `-c`, `--col-name`          | `None`     | Column name for CSV input. If not provided, input is treated as a text file with one SMILES per line           |
| `--rdkit-filter`            | enabled    | Apply RDKit structural alert filters                                                                          |
| `--no-rdkit-filter`         |            | Disable RDKit filters                                                                                         |
| `--rdkit-subset`            | `"ALL"`    | Subset of RDKit filters to apply (e.g., `PAINS`, `CHEMBL`, `BRENK`, `ALL`)                                   |
| `--rdkit-valtype`           | `"string"` | Output format: `"string"` for descriptions, `"bool"` for boolean flags                                       |
| `--pep-filter`              | enabled    | Apply peptide filters (PepSift)                                                                               |
| `--no-pep-filter`           |            | Disable peptide filters                                                                                       |
| `--silly-filter`            | enabled    | Apply silly molecule filters (molspotter)                                                                     |
| `--no-silly-filter`         |            | Disable silly molecule filters                                                                                |
| `--bloom-filter`            | enabled    | Apply purchasability filters (molbloom)                                                                       |
| `--no-bloom-filter`         |            | Disable purchasability filters                                                                                |
| `--std-mols`                | enabled    | Standardize molecules before filtering                                                                        |
| `--no-std-mols`             |            | Skip standardization                                                                                         |
| `--std-method`              | `"chembl"` | Standardization method (`chembl`, `canon`, `papyrus`, `molvs`)                                                |
| `--n-jobs`                  | `8`        | Number of parallel jobs                                                                                       |
| `--chunk-size`              | `-1`       | Process SMILES in chunks of this size. Negative means process all at once. Reduce if memory is limited         |

## Examples

Filter a text file of SMILES using all filters:

```bash
chemfilters -i molecules.txt -o results.csv
```

Filter a CSV file, disabling peptide and silly filters:

```bash
chemfilters -i data.csv -c smiles_column -o results.csv \
    --no-pep-filter --no-silly-filter
```

Use only PAINS filters with boolean output:

```bash
chemfilters -i molecules.txt -o results.csv \
    --rdkit-subset PAINS --rdkit-valtype bool \
    --no-pep-filter --no-silly-filter --no-bloom-filter
```
