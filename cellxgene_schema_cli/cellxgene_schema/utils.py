import logging
import os
import sys
from base64 import b85encode
from typing import Dict, List, Union

import anndata as ad
import numpy as np
from scipy import sparse
from xxhash import xxh3_64_intdigest

logger = logging.getLogger(__name__)

SPARSE_MATRIX_TYPES = {"csc", "csr", "coo"}


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
    adata = adata[:, var_to_keep]

    # Repeat much of the same steps for the raw.var, if it exists
    if adata.raw:
        raw_adata = ad.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs)
        var_to_keep = raw_adata.var.index[~raw_adata.var.index.isin(deprecated)].tolist()
        raw_adata = raw_adata[:, var_to_keep]
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


def get_matrix_format(adata: ad.AnnData, matrix: Union[np.ndarray, sparse.spmatrix]) -> str:
    """
    Given a matrix, returns the format as one of: csc, csr, coo, dense
    or unknown.

    This mimics the scipy.sparse `format` property, but extends it to
    support ndarray and other classes AnnData may proxy the matrix with.
    """

    # Note: the AnnData proxy classes DO support the `format_str` property, but
    # doing a slice seemed safer, if less performant.  Using `format_str`, which
    # currently works, uses private API:
    #
    # >>> return getattr(matrix, "format_str", "dense)
    #
    matrix_format = "unknown"
    if adata.n_obs == 0 or adata.n_vars == 0:
        matrix_format = "dense"
    else:
        matrix_slice = matrix[0:1, 0:1]
        if isinstance(matrix_slice, sparse.spmatrix):
            matrix_format = matrix_slice.format
        elif isinstance(matrix_slice, np.ndarray):
            matrix_format = "dense"

    assert matrix_format in ["unknown", "csr", "csc", "coo", "dense"]
    return matrix_format


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


def read_h5ad(h5ad_path: Union[str, bytes, os.PathLike]) -> ad.AnnData:
    """
    Reads h5ad into adata
    :params Union[str, bytes, os.PathLike] h5ad_path: path to h5ad to read

    :rtype None
    """
    try:
        adata = ad.read_h5ad(h5ad_path, backed="r")

        # This code, and AnnData in general, is optimized for row access.
        # Running backed, with CSC, is prohibitively slow. Read the entire
        # AnnData into memory if it is CSC.
        if (get_matrix_format(adata, adata.X) == "csc") or (
            (adata.raw is not None) and (get_matrix_format(adata, adata.raw.X) == "csc")
        ):
            logger.warning("Matrices are in CSC format; loading entire dataset into memory.")
            adata = adata.to_memory()

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
