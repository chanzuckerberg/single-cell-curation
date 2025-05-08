import anndata
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

    def test__obs_has_duplicate_rows_X(self, valid_adata):
        # zeros out the X matrix, resulting in duplicate rows of all 0's
        # this is valid because X is not the raw matrix
        valid_adata.X = from_array(
            sparse.csr_matrix((valid_adata.obs.shape[0], valid_adata.var.shape[0]), dtype=numpy.float32)
        )
        assert check_duplicate_obs(valid_adata) == []

    def test__obs_has_duplicate_rows_X_raw(self, valid_adata):
        # zeros out raw.X
        raw = anndata.AnnData(X=valid_adata.raw.X, obs=valid_adata.obs, var=valid_adata.raw.var)
        raw.X = from_array(sparse.csr_matrix((valid_adata.obs.shape[0], valid_adata.var.shape[0]), dtype=numpy.float32))
        valid_adata.raw = raw
        assert check_duplicate_obs(valid_adata) == ["Found 2 duplicated raw counts in obs adata.raw.X."]

    def test__obs_has_duplicate_rows_X_no_raw(self, valid_adata):
        # zeros out the X matrix, resulting in duplicate rows of all 0's
        # this is not valid because X is by default the raw matrix
        valid_adata.X = from_array(
            sparse.csr_matrix((valid_adata.obs.shape[0], valid_adata.var.shape[0]), dtype=numpy.float32)
        )
        valid_adata.raw = None
        assert check_duplicate_obs(valid_adata) == ["Found 2 duplicated raw counts in obs adata.X."]
