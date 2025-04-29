from anndata.compat import DaskArray
from dask.array import map_blocks
import sparse 
import numpy as np
from typing import Union
from .utils import (
    SPARSE_MATRIX_TYPES,
    get_matrix_format,
)

@staticmethod
def count_matrix_nonzero(matrix: DaskArray) -> int:
    def count_nonzeros(matrix_chunk: Union[np.ndarray, sparse.spmatrix], is_sparse_matrix: bool) -> np.array:
        nnz = matrix_chunk.nnz if is_sparse_matrix else np.count_nonzero(matrix_chunk)
        return np.array([nnz])

    is_sparse_matrix = get_matrix_format(matrix) in SPARSE_MATRIX_TYPES
    if len(matrix.chunks[0]) > 1:
        nonzeros = map_blocks(count_nonzeros, matrix, is_sparse_matrix, drop_axis=1, dtype=int).compute().sum()
    else:
        nonzeros = count_nonzeros(matrix.compute(), is_sparse_matrix)[0]
    return nonzeros

