import logging
import os
import sys
from base64 import b85encode
from functools import lru_cache
from typing import Dict, List, Union

import anndata as ad
import h5py
import numpy as np
from anndata.compat import DaskArray
from anndata.experimental import read_dispatched, read_elem_as_dask
from cellxgene_ontology_guide.ontology_parser import OntologyParser
from cellxgene_schema.matrix_utils import calculate_matrix_nonzero, check_non_csr_matrixes, determine_matrix_format
from xxhash import xxh3_64_intdigest

logger = logging.getLogger(__name__)

KB = 1024
MB = 1024 * KB
GB = 1024 * MB

"""
Ideally, these methods should all only live within matrix_utils. However, we
currently import these into single-cell-data-portal, so we need to keep these
in here for backwards compatibility until we can refactor that.
"""
SPARSE_MATRIX_TYPES = {"csr", "csc", "coo"}
SUPPORTED_SPARSE_MATRIX_TYPES = {"csr"}


def get_matrix_format(matrix: DaskArray) -> str:
    return determine_matrix_format(matrix)


def count_matrix_nonzero(matrix: DaskArray, is_sparse_matrix: bool) -> int:
    return calculate_matrix_nonzero(matrix, is_sparse_matrix)


def check_non_csr_matrices(adata: ad.AnnData):
    return check_non_csr_matrixes(adata)


def replace_ontology_term(dataframe, ontology_name, update_map):
    column_name = f"{ontology_name}_ontology_term_id"
    if dataframe[column_name].dtype != "category":
        dataframe[column_name] = dataframe[column_name].astype("category")
    for old_term, new_term in update_map.items():
        if old_term in dataframe[column_name].cat.categories:
            # add new one if not already in category, else continue
            if new_term not in dataframe[column_name].cat.categories:
                dataframe[column_name] = dataframe[column_name].cat.add_categories(new_term)
            # replace in dataset
            dataframe.loc[dataframe[column_name] == old_term, column_name] = new_term
            # remove deprecated_term from category
            dataframe[column_name] = dataframe[column_name].cat.remove_categories(old_term)


def replace_ontology_term_uns(adata: ad.AnnData, ontology_name, update_map) -> ad.AnnData:
    key_name = f"{ontology_name}_ontology_term_id"

    # TODO: We should remove this after the schema 6.0 migration
    if ontology_name == "organism" and key_name not in adata.uns:
        return move_ontology_term_from_obs_to_uns(adata, key_name, update_map)
    else:
        for old_term, new_term in update_map.items():
            if adata.uns[key_name] == old_term:
                adata.uns[key_name] = new_term
    return adata


def move_ontology_term_from_obs_to_uns(adata: ad.AnnData, key_name, update_map) -> ad.AnnData:
    if key_name not in adata.obs:
        raise KeyError(f"Column '{key_name}' not found in adata.obs, cannot migrate from obs to uns")

    values = adata.obs[key_name].unique()

    if len(values) != 1:
        raise ValueError(f"Cannot migrate from obs to uns because '{key_name}' has multiple values: {values}")

    adata.uns[key_name] = values[0]

    # Map old terms to new terms, if needed
    for old_term, new_term in update_map.items():
        if adata.uns[key_name] == old_term:
            adata.uns[key_name] = new_term

    return adata


def map_ontology_term(dataframe, ontology_name, map_from_column, update_map):
    column_name = f"{ontology_name}_ontology_term_id"
    if dataframe[column_name].dtype != "category":
        dataframe[column_name] = dataframe[column_name].astype("category")
    for map_value, new_term in update_map.items():
        if new_term not in dataframe[column_name].cat.categories:
            dataframe[column_name] = dataframe[column_name].cat.add_categories(new_term)
        dataframe.loc[dataframe[map_from_column] == map_value, column_name] = new_term
    dataframe[column_name] = dataframe[column_name].cat.remove_unused_categories()


def remove_deprecated_features(*, adata: ad.AnnData, deprecated: List[str]) -> ad.AnnData:
    # Filter out genes that don't appear in the approved annotation
    var_to_keep = adata.var.index[~adata.var.index.isin(deprecated)].tolist()
    adata = adata[:, var_to_keep].copy()

    # Repeat much of the same steps for the raw.var, if it exists
    if adata.raw:
        raw_adata = ad.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs)
        var_to_keep = raw_adata.var.index[~raw_adata.var.index.isin(deprecated)].tolist()
        raw_adata = raw_adata[:, var_to_keep].copy()
        adata.raw = raw_adata
    return adata


def remap_deprecated_features(*, adata: ad.AnnData, remapped_features: Dict[str, str]) -> ad.AnnData:
    # Use remapped_terms to map to the new term ids
    adata.var.index = [remapped_features.get(val, val) for val in adata.var.index]

    # Repeat much of the same steps for the raw.var, if it exists
    if adata.raw:
        raw_adata = ad.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs)
        raw_adata.var.index = [remapped_features.get(val, val) for val in raw_adata.var.index]
        adata.raw = raw_adata
    return adata


def getattr_anndata(adata: ad.AnnData, attr: str = None):
    """
    same as getattr but handles the special case of "raw.var" for an anndata.AndData object

    :param anndata.AnnData adata: the anndata.AnnData object from which to extract an attribute
    :param str attr: name of the attribute to extract

    :return the attribute or none if it does not exist
    """

    if attr == "raw.var":
        if adata.raw:
            return adata.raw.var
        else:
            return None
    else:
        return getattr(adata, attr)


def read_backed(f: h5py.File, chunk_size: int) -> ad.AnnData:
    """
    Read an AnnData object from a h5py.File object, reading in matrices (dense or sparse) as dask arrays. Does not
    read full matrices into memory.

    :param f: h5py.File object
    :param chunk_size: size of chunks to read matrices in
    :return: ad.AnnData object
    """

    def callback(func, elem_name: str, elem, iospec):
        if "/layers" in elem_name or "/uns" in elem_name or elem_name == "/X" or elem_name == "/raw/X":
            if iospec.encoding_type == "csr_matrix":
                n_vars = elem.attrs.get("shape")[1]
                return read_elem_as_dask(elem, chunks=(chunk_size, n_vars))
            elif iospec.encoding_type == "csc_matrix":
                n_obs = elem.attrs.get("shape")[0]
                return read_elem_as_dask(elem, chunks=(n_obs, chunk_size))
            elif iospec.encoding_type == "array" and len(elem.shape) == 2:
                n_vars = elem.shape[1]
                return read_elem_as_dask(elem, chunks=(chunk_size, n_vars))
            else:
                return func(elem)
        else:
            return func(elem)

    adata = read_dispatched(f, callback=callback)

    return adata


def read_h5ad(h5ad_path: Union[str, bytes, os.PathLike], chunk_size: int = 5000) -> ad.AnnData:
    """
    Reads h5ad into adata
    :params Union[str, bytes, os.PathLike] h5ad_path: path to h5ad to read

    :rtype None
    """
    try:
        f = h5py.File(h5ad_path)
        adata = read_backed(f, chunk_size)

    except (OSError, TypeError):
        logger.info(f"Unable to open '{h5ad_path}' with AnnData")
        sys.exit(1)

    return adata


def get_hash_digest_column(dataframe):
    """
    Get column with hash digest for each row in dataframe.
    """
    df_index = dataframe.index.to_series()
    assert df_index.is_unique
    return (
        df_index.map(xxh3_64_intdigest)
        .astype(np.uint64)
        .apply(lambda v: b85encode(v.to_bytes(8, "big")).decode("ascii"))
    )


@lru_cache()
def is_ontological_descendant_of(onto: OntologyParser, term: str, target: str, include_self: bool = True) -> bool:
    """
    Determines if :term is an ontological descendant of :target and whether to include :term==:target.

    This function is cached and is safe to call many times.

    #TODO:[EM] needs testing
    """
    return term in set(onto.get_term_descendants(target, include_self))


@lru_cache()
def get_descendants(onto: OntologyParser, term: str, include_self: bool = True) -> List[str]:
    return onto.get_term_descendants(term, include_self=True)
