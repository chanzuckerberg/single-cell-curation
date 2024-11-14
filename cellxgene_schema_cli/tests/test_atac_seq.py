from pathlib import Path
from typing import Tuple

import pytest
from cellxgene_schema import atac_seq
from fixtures.examples_validate import FIXTURES_ROOT


@pytest.fixture
def expected_artifacts() -> Tuple[Path, Path]:
    bgzip_file = Path(FIXTURES_ROOT + "/atac_seq/fragments_sorted.tsv.bgz")
    index_file = Path(FIXTURES_ROOT + "/atac_seq/fragments_sorted.tsv.bgz.tbi")
    yield (bgzip_file, index_file)
    bgzip_file.unlink(missing_ok=True)
    index_file.unlink(missing_ok=True)


def test_process_fragment(expected_artifacts):
    anndata_file = FIXTURES_ROOT + "/atac_seq/small_atac_seq.h5ad"
    fragments_file = FIXTURES_ROOT + "/atac_seq/fragments_sorted.tsv.gz"
    atac_seq.process_fragment(
        fragments_file,
        anndata_file,
        generate_index=True,
        dask_cluster_config=dict(processes=False),
    )
    for artifact in expected_artifacts:
        assert artifact.exists()
