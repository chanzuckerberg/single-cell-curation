import hashlib
import numpy as np
import pandas as pd
from scipy import sparse

from utils import (
    SPARSE_MATRIX_TYPES,
    get_matrix_format,
)

def check_duplicate_obs(adata) -> list[str]:
    """
    Checks for duplicate rows in obs by raw count.

    :rtype List of errors that were seen. If no errors, validation is good.
    """
    matrix = adata.X[:, column]
    is_sparse_matrix = get_matrix_format(matrix) in SPARSE_MATRIX_TYPES
    if is_sparse_matrix:
        return _check_duplicate_obs_sparse(adata)
    else:
        return _check_duplicate_obs_dense(adata)

def _check_duplicate_obs_dense(adata) -> list[str]:
    errors = []

    if 'in_tissue' in adata.obs.columns:
      obs_to_keep = adata.obs[adata.obs['in_tissue'] != 0].index
      adata = adata[obs_to_keep, :]

    matrix = adata.raw.X if adata.raw else adata.X

    if hasattr(matrix, 'toarray'):
        matrix = matrix.toarray()
    elif not isinstance(matrix, np.ndarray):
        matrix = np.array(matrix)

    def row_hash(row):
        return hashlib.sha1(row.tobytes()).hexdigest()

    row_hashes = [row_hash(row) for row in matrix]

    hash_df = adata.obs.copy()
    hash_df['row_hash'] = row_hashes

    dup_df = hash_df[hash_df.duplicated(subset='row_hash', keep=False)].copy()

    if not dup_df.empty:
        errors.append('duplicated raw counts', 'ERROR')

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

    if 'in_tissue' in adata.obs.columns:
        obs_to_keep = adata.obs[adata.obs['in_tissue'] != 0].index
        adata = adata[obs_to_keep, : ]

    matrix = adata.raw.X if adata.raw else adata.X

    if not isinstance(matrix, sparse.csr_matrix):
        errors.append("Matrix not in sparse csr format, please convert before hashing")
        return

    nnz = matrix.nnz

    if not matrix.has_canonical_format:
        if adata.raw:
            adata.raw.X.sort_indices()
            adata.raw.X.sum_duplicates()
        else:
            adata.X.sort_indices()
            adata.X.sum_duplicates()

    assert matrix.has_canonical_format, "Matrix still in non-canonical format"

    if nnz != matrix.nnz:
        errors.append(f"{nnz - matrix.nnz} duplicates found during canonical conversion")

    data_array = matrix.data
    index_array = matrix.indices
    indptr_array = matrix.indptr

    start, end = 0, matrix.shape[0]
    hashes = []
    while start < end:
        val = hash(data_array[indptr_array[start]:indptr_array[start + 1]].tobytes())
        hashes.append(val)
        start += 1

    def index_hash(index):
        obs_loc = adata.obs.index.get_loc(index)
        val = hash(index_array[indptr_array[obs_loc]:indptr_array[obs_loc + 1]].tobytes())

        return val
    
    hash_df = adata.obs.copy()
    hash_df['data_array_hash'] = hashes
    hash_df = hash_df[hash_df.duplicated(subset='data_array_hash',keep=False) == True]
    hash_df.sort_values('data_array_hash', inplace=True)

    hash_df['index_array_hash'] = [index_hash(row) for row in hash_df.index.to_list()]
    hash_df = hash_df[hash_df.duplicated(subset=['data_array_hash', 'index_array_hash'], keep=False) == True]

    if not hash_df.empty:
        errors.append('duplicated raw counts')
    
    return errors