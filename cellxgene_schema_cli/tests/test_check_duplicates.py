import numpy
import pytest
from cellxgene_schema.validation_internals.check_duplicates import check_duplicate_obs
from dask.array import from_array
from fixtures.examples_validate import adata
from scipy import sparse


@pytest.fixture
def valid_adata():
    return adata.copy()


class TestCheckDuplicates:
    def test__check_duplicate_obs_ok(self, valid_adata):
        assert check_duplicate_obs(valid_adata) == []

    def test__obs_has_duplicate_rows(self, valid_adata):
        # zeros out the X matrix, resulting in duplicate rows of all 0's
        valid_adata.X = from_array(
            sparse.csr_matrix((valid_adata.obs.shape[0], valid_adata.var.shape[0]), dtype=numpy.float32)
        )
        assert check_duplicate_obs(valid_adata) == [
            "Found 2 duplicated raw counts in obs adata.X. First 2 duplicate rows found at: ['row 0: index = X', 'row 1: index = Y']."
        ]
