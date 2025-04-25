
import pytest

from fixtures.examples_validate import adata as adata_valid
from cellxgene_schema.validation_internals.check_duplicates import check_duplicate_obs

@pytest.fixture
def valid_adata():
    return adata_valid.copy()

class TestCheckDuplicates:
    def test__check_duplicate_obs(self, valid_adata):
        assert check_duplicate_obs(valid_adata) == []