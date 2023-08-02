import argparse
from pathlib import Path

import pandas as pd

from chemFilters.core import CoreFilter
from chemFilters.filters.rdkit_filters import FILTER_COLLECTIONS


def main():
    parser = argparse.ArgumentParser(
        prog="chemFilters",
        description="Filters chemical compounds based on various criteria.",
    )

    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to the input file containing SMILES strings.",
    )
    parser.add_argument(
        "-c",
        "--col-name",
        type=str,
        required=False,
        default=None,
        help=(
            "column name for the csv file containing SMILES strings. "
            "If not provided, will treat input as a txt file with SMILES only"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to the output file to write filtered SMILES strings.",
    )
    parser.add_argument(
        "--rdkit-filter",
        dest="rdkit_filter",
        action="store_true",
        help="Apply RDKit filters. Enabled by default.",
    )
    parser.add_argument(
        "--no-rdkit-filter",
        dest="rdkit_filter",
        action="store_false",
        help="Disable RDKit filters.",
    )
    parser.add_argument(
        "--rdkit-subset",
        type=str,
        default="ALL",
        help=(
            "Subset of the RDKit filters to apply. Example subsets are "
            f"{FILTER_COLLECTIONS}. Defaults to 'ALL'."
        ),
    )
    parser.add_argument(
        "--rdkit-valtype",
        type=str,
        default="string",
        help=(
            "If 'string', will return the description of the match. "
            "If 'bool', will return booleans on whether it matched. "
            " Defaults to `string`."
        ),
    )
    parser.add_argument(
        "--pep-filter",
        dest="pep_filter",
        action="store_true",
        help="Apply peptide filters. Enabled by default.",
    )
    parser.add_argument(
        "--no-pep-filter",
        dest="pep_filter",
        action="store_false",
        help="Disable peptide filters.",
    )
    parser.add_argument(
        "--silly-filter",
        dest="silly_filter",
        action="store_true",
        help="Apply silly filters. Enabled by default.",
    )
    parser.add_argument(
        "--no-silly-filter",
        dest="silly_filter",
        action="store_false",
        help="Disable silly filters.",
    )
    parser.add_argument(
        "--bloom-filter",
        dest="bloom_filter",
        action="store_true",
        help="Apply bloom filters. Enabled by default.",
    )
    parser.add_argument(
        "--no-bloom-filter",
        dest="bloom_filter",
        action="store_false",
        help="Disable bloom filters.",
    )
    parser.add_argument(
        "--std-mols",
        dest="std_mols",
        action="store_true",
        help="Standardize the mols. Enabled by default.",
    )
    parser.add_argument(
        "--no-std-mols",
        dest="std_mols",
        action="store_false",
        help="Do not standardize the mols.",
    )
    parser.add_argument(
        "--std-method",
        type=str,
        default="chembl",
        help='Method to standardize mols. Defaults to "chembl".',
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=8,
        help="Number of jobs to run in parallel. Defaults to 8.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=-1,
        help=(
            "Number of smiles to be processed in chunks. If negative, will process "
            "all smiles at once. If the memory is not enough, reduce this number. "
            "Defaults to -1."
        ),
    )
    parser.set_defaults(
        rdkit_filter=True,
        pep_filter=True,
        silly_filter=True,
        bloom_filter=True,
        std_mols=True,
    )
    args = parser.parse_args()
    core_filter = CoreFilter(
        rdkit_filter=args.rdkit_filter,
        pep_filter=args.pep_filter,
        silly_filter=args.silly_filter,
        bloom_filter=args.bloom_filter,
        rdfilter_subset=args.rdkit_subset,
        rdfilter_output=args.rdkit_valtype,
        std_mols=args.std_mols,
        std_method=args.std_method,
        n_jobs=args.n_jobs,
    )
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file {input_path} does not exist.")

    if args.col_name is None:
        with input_path.open() as f:
            smiles = [line.strip() for line in f]
        result = core_filter(smiles, chunksize=args.chunk_size)
    else:
        smiles = pd.read_csv(input_path, usecols=[args.col_name]).values
        result = core_filter(smiles, chunksize=args.chunk_size)
    result.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
