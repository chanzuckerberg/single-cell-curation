import pytest
from cellxgene_schema.validation_internals.check_duplicates import check_duplicate_obs
from fixtures.examples_validate import adata_minimal, adata


@pytest.fixture
def valid_adata():
    return adata.copy()


class TestCheckDuplicates:
    def test__check_duplicate_obs(self, valid_adata):
        assert check_duplicate_obs(valid_adata) == []
