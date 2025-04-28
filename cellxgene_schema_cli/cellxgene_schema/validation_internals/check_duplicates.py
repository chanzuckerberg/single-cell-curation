import hashlib

import anndata
import dask.array as da
import numpy as np
from scipy import sparse
from utils import get_matrix_format


def check_duplicate_obs(adata: anndata.AnnData) -> list[str]:
    """
    Checks for duplicate rows in obs by raw count.

    :rtype List of errors that were seen. If no errors, validation is good.
    """

    errors = []

    if "in_tissue" in adata.obs.columns:
        obs_to_keep = adata.obs[adata.obs["in_tissue"] != 0].index
        adata = adata[obs_to_keep, :]

    to_validate = [(adata.X, "adata.X")]

    if adata.raw:
        to_validate.append((adata.raw.X, "adata.raw.X"))

    for matrix, matrix_name in to_validate:
        matrix_format = get_matrix_format(matrix)
        if matrix_format == "csr" or matrix_format == "dense":
            errors.extend(_check_duplicate_obs_dask(adata, matrix, matrix_name))
        else:
            errors.append(
                f"Unsupported matrix format: {matrix_format} Sparse matrices must be encoded as scipy.sparse.csr_matrix. Dense matrices can be np.ndarray"
            )

    return errors


def _check_duplicate_obs_dask(adata: anndata.AnnData, matrix: da.Array, matrix_name: str) -> list[str]:
    """
    For each row in the matrix, we compute a hash of the row and store it in row_hashes.
    We then check for duplicates in the row_hashes.
    """
    errors = []

    if matrix.ndim != 2:
        return [f"Dask matrix {matrix_name} must be 2D, but got shape {matrix.shape}"]

    def rowwise_hash_block(block):
        # If it's sparse, convert each block to dense array so that we can run tobytes()
        if sparse.issparse(block):
            block = block.toarray()

        return np.array([hashlib.sha224(row.tobytes()).hexdigest() for row in block], dtype=object).reshape(
            -1, 1
        )  # shape (rows, 1) for map_blocks

    row_hashes = (
        matrix.map_blocks(
            rowwise_hash_block,
            dtype=object,
            chunks=(matrix.chunks[0], (1,)),  # one hash per row
        )
        .compute()
        .ravel()
    )

    hash_df = adata.obs.copy()
    hash_df["row_hash"] = row_hashes
    dup_df = hash_df[hash_df.duplicated(subset="row_hash", keep=False)].copy()

    if not dup_df.empty:
        duplicate_rows_to_print = []
        for i, idx in enumerate(hash_df.index[:10]):
            duplicate_rows_to_print.append(f"row {i}: index = {idx}")
        errors.append(
            f"Found {len(dup_df)} duplicated raw counts in obs {matrix_name}. First {len(duplicate_rows_to_print)} duplicate rows found at: {duplicate_rows_to_print}. {hash_df}"
        )

    return errors
