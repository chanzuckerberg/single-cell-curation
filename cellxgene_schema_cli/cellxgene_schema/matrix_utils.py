from typing import Union

import anndata as ad
import numpy as np
import sparse as sp
from anndata.compat import DaskArray
from dask.array import map_blocks
from scipy import sparse

SPARSE_MATRIX_TYPES = {"csr", "csc", "coo"}
SUPPORTED_SPARSE_MATRIX_TYPES = {"csr"}


def determine_matrix_format(matrix: DaskArray) -> str:
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
    try:
        matrix_slice = matrix[0:1, 0:1].compute()
    except AttributeError:
        # compute() may fail on an unknown matrix value. if so, return "unknown"
        return matrix_format
    if isinstance(matrix_slice, sparse.spmatrix):
        matrix_format = matrix_slice.format
    elif isinstance(matrix_slice, np.ndarray):
        matrix_format = "dense"
    assert matrix_format in SPARSE_MATRIX_TYPES.union({"unknown", "dense"})
    return matrix_format


def compute_column_sums(matrix: Union[DaskArray, np.ndarray, sparse.spmatrix]) -> np.ndarray:
    """
    Compute column-wise sums for a Dask array (dense or sparse), NumPy ndarray, or SciPy sparse matrix.
    Returns a NumPy array of sums.

    For example, this matrix:
    [
        [1, 0, 0],
        [1, 2, 0],
        [1, 0, 1],
    ]

    would return a NumPy array of [3, 12, 1]
    """
    if isinstance(matrix, np.ndarray):
        return matrix.sum(axis=0).ravel()

    if sparse.issparse(matrix):
        return np.array(matrix.sum(axis=0)).ravel()

    # Handle Dask array (could be sparse or dense)
    if isinstance(matrix, DaskArray):

        def to_coo(chunk: Union[np.ndarray, sparse.spmatrix]) -> Union[sp.COO, np.ndarray]:
            # convert to sparse.COO, so we can use dask's internal chunked sum
            if sparse.issparse(chunk):
                return sp.COO.from_scipy_sparse(chunk)
            else:
                return chunk

        x_coo = matrix.map_blocks(to_coo, dtype=matrix.dtype)

        # 'split_every' is used to ensure that the computation is split into memory manageable chunks
        col_sums = x_coo.sum(axis=0, split_every=8).compute()
        return col_sums


def calculate_matrix_nonzero(matrix: DaskArray) -> int:
    def count_nonzeros(matrix_chunk: Union[np.ndarray, sparse.spmatrix], is_sparse_matrix: bool) -> np.array:
        nnz = matrix_chunk.nnz if is_sparse_matrix else np.count_nonzero(matrix_chunk)
        return np.array([nnz])

    is_sparse_matrix = determine_matrix_format(matrix) in SPARSE_MATRIX_TYPES
    if len(matrix.chunks[0]) > 1:
        nonzeros = map_blocks(count_nonzeros, matrix, is_sparse_matrix, drop_axis=1, dtype=int).compute().sum()
    else:
        nonzeros = count_nonzeros(matrix.compute(), is_sparse_matrix)[0]
    return nonzeros


def check_non_csr_matrixes(adata: ad.AnnData):
    """
    Check X, raw.X and layers matrices for having more than 50% zeros and not being csr_matrix

    If found, convert to csr_matrix
    """

    def get_sparsity(matrix: DaskArray, format: str):
        is_sparse_matrix = format in SPARSE_MATRIX_TYPES
        nnz = calculate_matrix_nonzero(matrix, is_sparse_matrix)
        sparsity = 1 - nnz / np.prod(matrix.shape)
        return sparsity

    format = determine_matrix_format(adata.X)
    if format != "csr" and get_sparsity(adata.X, format) >= 0.5:
        adata.X = adata.X.map_blocks(sparse.csr_matrix, dtype=adata.X.dtype)

    if adata.raw is not None:
        format = determine_matrix_format(adata.raw.X)
        if format != "csr" and get_sparsity(adata.raw.X, format) >= 0.5:
            raw_adata = ad.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs)
            raw_adata.X = raw_adata.X.map_blocks(sparse.csr_matrix, dtype=raw_adata.X.dtype)
            adata.raw = raw_adata
            del raw_adata

    for layer in adata.layers:
        format = determine_matrix_format(adata.layers[layer])
        if format != "csr" and get_sparsity(adata.layers[layer], format) >= 0.5:
            adata.layers[layer] = adata.layers[layer].map_blocks(sparse.csr_matrix, dtype=adata.layers[layer].X.dtype)

    return adata


def debug_print_matrix(matrix: ad.AnnData, matrix_name: str, max_rows=20, max_cols=20):
    computed = matrix[:max_rows, :max_cols].compute()

    # If sparse, convert to dense for display
    if hasattr(computed, "toarray"):
        computed = computed.toarray()

    print(f"{matrix_name} (first {max_rows} rows, {max_cols} cols):")
    print(np.round(computed, 2))  # round for readability
