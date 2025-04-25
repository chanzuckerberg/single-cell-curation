import hashlib

import numpy as np
from scipy import sparse
from utils import (
    SPARSE_MATRIX_TYPES,
    get_matrix_format,
)


def check_duplicate_obs(adata) -> list[str]:
    """
    Checks for duplicate rows in obs by raw count.

    The logic is split by whether it's a sparse or dense matrix.

    :rtype List of errors that were seen. If no errors, validation is good.
    """
    matrix = adata.raw.X if adata.raw else adata.X
    is_sparse_matrix = get_matrix_format(matrix) in SPARSE_MATRIX_TYPES
    if is_sparse_matrix:
        return _check_duplicate_obs_sparse(adata)
    else:
        return _check_duplicate_obs_dense(adata)


def _check_duplicate_obs_dense(adata) -> list[str]:
    errors = []

    if "in_tissue" in adata.obs.columns:
        obs_to_keep = adata.obs[adata.obs["in_tissue"] != 0].index
        adata = adata[obs_to_keep, :]

    matrix = adata.raw.X if adata.raw else adata.X

    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    elif not isinstance(matrix, np.ndarray):
        matrix = np.array(matrix)

    def row_hash(row):
        return hashlib.sha224(row.tobytes()).hexdigest()

    row_hashes = [row_hash(row) for row in matrix]

    hash_df = adata.obs.copy()
    hash_df["row_hash"] = row_hashes

    dup_df = hash_df[hash_df.duplicated(subset="row_hash", keep=False)].copy()

    if not dup_df.empty:
        errors.append(f"Found {len(dup_df)} duplicated raw counts in obs")

    return errors


def _check_duplicate_obs_sparse(adata) -> list[str]:
    """
    Hash sparse csr matrix using np.ndarrays that represent sparse matrix data.
    First pass will hash all rows via slicing the data array and append to copy of obs df
    Second pass will hash only duplicate rows in obs copy via the indices array.
    This will keep only true duplicated matrix rows and not rows with an indicental same
    ordering of their data arrays
    """
    errors = []

    if "in_tissue" in adata.obs.columns:
        obs_to_keep = adata.obs[adata.obs["in_tissue"] != 0].index
        adata = adata[obs_to_keep, :]

    matrix = adata.raw.X if adata.raw else adata.X

    if not isinstance(matrix, sparse.csr_matrix):
        errors.append("Matrix not in sparse csr format, please convert before hashing")
        return errors

    if not matrix.has_canonical_format:
        matrix.sort_indices()
        matrix.sum_duplicates()

    assert matrix.has_canonical_format, "Matrix still in non-canonical format"

    indptr = matrix.indptr
    data = matrix.data
    indices = matrix.indices

    def row_hash(data_slice):
        return hashlib.sha224(data_slice.tobytes()).hexdigest()

    data_hashes = [row_hash(data[indptr[i] : indptr[i + 1]]) for i in range(matrix.shape[0])]
    index_hashes = [row_hash(indices[indptr[i] : indptr[i + 1]]) for i in range(matrix.shape[0])]

    hash_df = adata.obs.copy()
    hash_df["data_hash"] = data_hashes
    hash_df["index_hash"] = index_hashes

    dup_mask = hash_df.duplicated(subset=["data_hash", "index_hash"], keep=False)
    dup_df = hash_df[dup_mask].copy()

    if not dup_df.empty:
        errors.append(f"Found {len(dup_df)} duplicated raw counts in obs")

    return errors
