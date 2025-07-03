import hashlib
import logging

import anndata
import dask.array as da
import numpy as np
from cellxgene_schema.matrix_utils import determine_matrix_format
from scipy import sparse

logger = logging.getLogger(__name__)


def check_duplicate_obs(adata: anndata.AnnData) -> list[str]:
    """
    Checks for duplicate rows in obs by raw count.

    We only want to check duplicate rows for the 'raw' matrix. If adata.raw is
    not provided, then we assume adata.X is the raw data. If both are provided,
    then we only want to validate against adata.raw.X

    :rtype List of errors that were seen. If no errors, validation is good.
    """

    errors = []

    if "in_tissue" in adata.obs.columns:
        obs_to_keep = adata.obs[adata.obs["in_tissue"] != 0].index
        adata = adata[obs_to_keep, :]

    to_validate = [(adata.X, "adata.X")]

    if adata.raw is not None:
        to_validate = [(adata.raw.X, "adata.raw.X")]

    for matrix, matrix_name in to_validate:
        matrix_format = determine_matrix_format(matrix)
        if matrix_format == "csr" or matrix_format == "dense":
            errors.extend(_check_duplicate_obs_dask(adata, matrix, matrix_name))
        else:
            # We already raise an error for unsupported formats in _validate_sparsity, so no need to
            # check again here.
            continue

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

    logger.debug("Computing row hashes for matrix %s", matrix_name)
    row_hashes = (
        matrix.map_blocks(
            rowwise_hash_block,
            dtype=object,
            chunks=(matrix.chunks[0], (1,)),  # one hash per row
        )
        .compute()
        .ravel()
    )

    logger.debug("Copying hash_df for matrix %s", matrix_name)
    hash_df = adata.obs.copy()
    hash_df["row_hash"] = row_hashes
    logger.debug("Checking for duplicates in row_hashes for matrix %s", matrix_name)
    dup_df = hash_df[hash_df.duplicated(subset="row_hash", keep=False)].copy()

    if not dup_df.empty:
        errors.append(f"Found {len(dup_df)} duplicated raw counts in obs {matrix_name}.")

    logger.debug("Done checking for duplicates in row_hashes for matrix %s", matrix_name)
    return errors
