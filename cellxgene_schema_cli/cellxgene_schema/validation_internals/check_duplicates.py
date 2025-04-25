import hashlib

import dask.array as da
import numpy as np
from scipy import sparse
from utils import get_matrix_format


def check_duplicate_obs(adata) -> list[str]:
    """
    Checks for duplicate rows in obs by raw count.

    The logic is split by whether it's a sparse or dense matrix.

    :rtype List of errors that were seen. If no errors, validation is good.
    """

    if "in_tissue" in adata.obs.columns:
        obs_to_keep = adata.obs[adata.obs["in_tissue"] != 0].index
        adata = adata[obs_to_keep, :]

    to_validate = [adata.X]

    if adata.raw:
        to_validate.append(adata.raw.X)
    
    for matrix in to_validate:
        matrix_format = get_matrix_format(matrix)
        if matrix_format == "csr" or matrix_format == "dense":
            return _check_duplicate_obs_dask(adata, matrix)
        else:
            return [
                f"Unsupported matrix format: {matrix_format} Sparse matrices must be encoded as scipy.sparse.csr_matrix. Dense matrices can be np.ndarray"
            ]

def _check_duplicate_obs_dask(adata, matrix: da.Array) -> list[str]:
    errors = []

    if matrix.ndim != 2:
        return [f"Dask matrix must be 2D, but got shape {matrix.shape}"]

    def rowwise_hash_block(block):
        # If it's sparse, convert each block to dense array
        if sparse.issparse(block):
            block = block.toarray()

        return np.array([
            hashlib.sha224(row.tobytes()).hexdigest()
            for row in block
        ], dtype=object).reshape(-1, 1)  # shape (rows, 1) for map_blocks

    row_hashes = matrix.map_blocks(
        rowwise_hash_block,
        dtype=object,
        chunks=(matrix.chunks[0], (1,)),  # one hash per row
    ).compute().ravel()

    hash_df = adata.obs.copy()
    hash_df["row_hash"] = row_hashes
    dup_df = hash_df[hash_df.duplicated(subset="row_hash", keep=False)].copy()

    if not dup_df.empty:
        errors.append(f"Found {len(dup_df)} duplicated raw counts in obs")

    return errors
