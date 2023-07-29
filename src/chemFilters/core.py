from multiprocessing import Pool
from typing import List

import numpy as np
import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from chemFilters import (
    MolbloomFilters,
    PeptideFilters,
    RdkitFilters,
    SillyMolSpotterFilter,
)
from chemFilters.chem.interface import mol_from_smi
from chemFilters.chem.standardizers import ChemStandardizer

STANDARD_RD_COLS = [
    "NIH",
    "PAINS_A",
    "PAINS_B",
    "PAINS_C",
    "PAINS_any",
    "ZINC",
    "Brenk",
    "ChEMBL23_Dundee",
    "ChEMBL23_BMS",
    "ChEMBL23_MLSMR",
    "ChEMBL23_Inpharmatica",
    "ChEMBL23_SureChEMBL",
    "ChEMBL23_LINT",
    "ChEMBL23_Glaxo",
    "ChEMBL23_any",
]

STANDARD_PEP_COLS = [
    "NaturalLAminoAcids",
    "NaturalLDAminoAcids",
    "NaturalAminoAcidDerivatives",
    "NonNaturalAminoAcidDerivatives",
    "AllAmineAndAcid",
]

STANDARD_BLOOM_COLS = [
    "zinc20",
    "zinc-instock",
    "zinc-instock-mini",
    "surechembl",
]

STANDARD_SILLY_COLS = [
    "chembl",
    "excape",
    "papyrus",
]


class CoreFilter:
    """Class implementation to run all filters on a list of smiles, with the option of
    adding a chunk size to process the input in batches. The filtering the dataset is
    done by:

    >>> from chemFilters.core import CoreFilter
    >>> core_filter = CoreFilter() # all filters enabled by default
    >>> filtered_df = core_filter(smiles, chunksize=100)

    Arguments:
        rdkit_filter: toggle applying rdkit filters to smiles. Defaults to True.
        pep_filter: toggle applying peptide filters to smiles. Defaults to True.
        silly_filter: toggle applying silly filters to smiles. Defaults to True.
        bloom_filter: toggle applying bloom filters to smiles. Defaults to True.
        std_mols: whether to standardize the mols. Defaults to False.
        std_method: standardization method to be used. Defaults to "chembl".
        n_jobs: number of jobs to run in parallel. Defaults to 1.
    """

    def __init__(
        self,
        rdkit_filter: bool = True,
        pep_filter: bool = True,
        silly_filter: bool = True,
        bloom_filter: bool = True,
        std_mols: bool = True,
        std_method: str = "chembl",
        n_jobs=1,
    ) -> None:
        """Initialize the CoreFilter class. The filters are initialized with the default
        parameters.

        Args:
            rdkit_filter: toggle applying rdkit filters to smiles. Defaults to True.
            pep_filter: toggle applying peptide filters to smiles. Defaults to True.
            silly_filter: toggle applying silly filters to smiles. Defaults to True.
            bloom_filter: toggle applying bloom filters to smiles. Defaults to True.
            std_mols: whether to standardize the mols. Defaults to True.
            std_method: standardization method to be used. Defaults to "chembl".
            n_jobs: number of jobs to run in parallel. Defaults to 1.
        """
        self.toggle_rdkit_filter = rdkit_filter
        self.toggle_pep_filter = pep_filter
        self.toggle_silly_filter = silly_filter
        self.toggle_bloom_filter = bloom_filter
        self.std_mols = std_mols
        self.n_jobs = n_jobs
        self._configure_filters(std_method, n_jobs)

    def _configure_filters(self, std_method: str, n_jobs: int):
        """Configure the filters. Right now it only sets them up, but could be adjusted
        for more flexibility in the future."""
        if self.toggle_rdkit_filter:
            self.rdFilter = RdkitFilters(
                filter_type="ALL", from_smi=False, n_jobs=n_jobs
            )
        if self.toggle_pep_filter:
            self.pepFilter = PeptideFilters(from_smi=False)
        if self.toggle_bloom_filter:
            self.bloomFilter = MolbloomFilters(from_smi=False, standardize=False)
        if self.toggle_silly_filter:
            self.sillyFilter = SillyMolSpotterFilter(from_smi=False, standardize=False)
        if self.std_mols:
            self.molStandardizer = ChemStandardizer(
                method=std_method,
                from_smi=False,
                return_smi=False,
                n_jobs=n_jobs,
                verbose=False,
            )

    def _get_chunks(self, smiles: List[str], chunksize: int) -> List[List[str]]:
        """Get the smiles in chunks."""
        if chunksize < 0 or chunksize is None:
            chunksize = len(smiles)
        smiles_chunks = [
            smiles[i : i + chunksize] for i in range(0, len(smiles), chunksize)
        ]
        return smiles_chunks

    def _yield_filtered_smiles(self, smiles, chunksize=1000):
        """Process a list of smiles in chunks."""
        smiles_chunks = self._get_chunks(smiles, chunksize)
        for chunk in smiles_chunks:
            with Pool(self.n_jobs) as p:
                mols = p.map(mol_from_smi, chunk)
            if self.std_mols:
                mols = self.molStandardizer(mols)

            filtered_dfs = list()
            if self.toggle_rdkit_filter:
                filtered_dfs.append(self._rdfilterMols(mols))
            if self.toggle_pep_filter:
                filtered_dfs.append(self._pepfilterMols(mols))
            if self.toggle_bloom_filter:
                filtered_dfs.append(self._bloomfilterMols(mols))
            if self.toggle_silly_filter:
                filtered_dfs.append(self._sillyfilterMols(mols))
            if len(filtered_dfs) > 1:
                filtered_dfs = [filtered_dfs[0]] + [
                    df.drop(columns="SMILES") for df in filtered_dfs[1:]
                ]
            concat_df = pd.concat(
                filtered_dfs,
                axis=1,
            )
            yield concat_df

    def _rdfilterMols(self, mols: List[Chem.Mol]):
        """Filter with self.rdFilter."""

        def aggregate_str_vals(df, cols):
            """function to aggregate strings within the columns into the same row. It
            takes into account nan values by replacing with empty strings, joining with
            `; `, and replacing the resulting `; ` * ncols - 1 with nan."""
            n_cols = len(cols)
            null_str = "; " * (
                n_cols - 1
            )  # after join with ';', this is the null string
            aggr_series = (
                df.loc[:, cols]
                .fillna("")
                .apply(lambda x: "; ".join(x), axis=1)
                .replace({null_str: np.nan})
            )
            return aggr_series

        rdfilt_df = self.rdFilter.get_flagging_df(mols)
        pains_cols = [col for col in rdfilt_df.columns if "PAINS" in col]
        chembl_cols = [col for col in rdfilt_df.columns if "ChEMBL" in col]
        any_pains = aggregate_str_vals(rdfilt_df, pains_cols)
        any_chembl = aggregate_str_vals(rdfilt_df, chembl_cols)
        return rdfilt_df.assign(PAINS_any=any_pains, CHEMBL23_any=any_chembl).reindex(
            columns=["SMILES"] + STANDARD_RD_COLS, fill_value=np.nan
        )

    def _pepfilterMols(self, mols: List[Chem.Mol]):
        """Filter with self.pepFilter."""
        pepfilt_df = self.pepFilter.get_flagging_df(mols)
        return pepfilt_df.assign(
            AminoAcid_any=lambda x: x[STANDARD_PEP_COLS].any(axis=1)
        )

    def _bloomfilterMols(self, mols: List[Chem.Mol]):
        """Filter with self.bloomFilter."""
        rename_dict = {"BLOOM_{}".format(name): name for name in STANDARD_BLOOM_COLS}
        return (
            self.bloomFilter.get_flagging_df(mols)
            .assign(BLOOM_any=lambda x: x[STANDARD_BLOOM_COLS].any(axis=1))
            .rename(columns=rename_dict)
        )

    def _sillyfilterMols(self, mols: List[Chem.Mol]):
        """Filter with self.sillyFilter."""
        rename_dict = {"SILLY_{}".format(name): name for name in STANDARD_SILLY_COLS}
        return self.sillyFilter.get_scoring_df(mols).rename(columns=rename_dict)

    def filter_smiles(self, smiles: list, chunksize: int = 1000):
        """Filter a list of smiles based on the filters that are toggled on. Will use
        the chunksize and n_jobs to process the data in chunks & parallel.

        Args:
            smiles: list of smiles to be filtered.
            chunksize: size of chunks to process the data. Set to None or < 0 for
                processing all at once. Defaults to 1000.

        Returns:
            pd.DataFrame: dataframe with the filtered dataframes."""
        n_chunks = len(self._get_chunks(smiles, chunksize))
        return pd.concat(
            tqdm(self._yield_filtered_smiles(smiles, chunksize), total=n_chunks),
            ignore_index=True,
        ).reset_index(drop=True)

    def __call__(self, smiles: list, chunksize: int = 1000):
        """Wrapper around filter_smiles. Will call _yield_filtered_smiles and on chuncks
        of a list of smiles and return a concatenated dataframe of the results.

        Args:
            smiles: list of smiles strings to be filtered.
            chunksize: size of chunks to process the data. Set to None or < 0 for
                processing all at once. Defaults to 1000.

        Returns:
            pd.DataFrame: dataframe with the filtered dataframes."""
        n_chunks = len(self._get_chunks(smiles, chunksize))
        return pd.concat(
            tqdm(self._yield_filtered_smiles(smiles, chunksize), total=n_chunks),
            ignore_index=True,
        ).reset_index(drop=True)
