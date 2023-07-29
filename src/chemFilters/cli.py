import argparse
from pathlib import Path

from chemFilters.core import CoreFilter


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
        help="Standardize the mols. Disabled by default.",
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
            "all smiles at once. If the memory is not enough, reduce this number."
            "Defaults to -1.",
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
        std_mols=args.std_mols,
        std_method=args.std_method,
        n_jobs=args.n_jobs,
    )
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file {input_path} does not exist.")
    with input_path.open() as f:
        smiles = [line.strip() for line in f]
    result = core_filter(smiles, chunksize=args.chunk_size)
    result.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
