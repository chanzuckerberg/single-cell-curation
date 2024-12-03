from pathlib import Path
from typing import Tuple
from unittest import mock

import anndata as ad
import dask.dataframe as dd
import pandas as pd
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


@pytest.fixture
def fragment_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "barcode": ["A", "B", "C"],
            "start": [100, 200, 300],
            "end": [200, 300, 400],
            "read support": [1, 2, 3],
            "chromosome": ["chr1", "chr2", "chr3"],
        }
    )


@pytest.fixture
def mock_read_anndata():
    with mock.patch("anndata.read_h5ad") as mock_anndata:
        mock_anndata.return_value = ad.AnnData(obs=pd.DataFrame(index=["A", "B", "C"]), var=pd.DataFrame())
        yield mock_anndata


@pytest.fixture
def mock_read_parquet(fragment_dataframe):
    with mock.patch("dask.dataframe.read_parquet") as mock_read_parquet:
        mock_read_parquet.return_value = dd.from_pandas(fragment_dataframe, npartitions=1)
        yield mock_read_parquet


def test_process_fragment(expected_artifacts):
    anndata_file = FIXTURES_ROOT + "/atac_seq/small_atac_seq.h5ad"
    fragments_file = FIXTURES_ROOT + "/atac_seq/fragments_sortedasd.tsv.gz"
    atac_seq.process_fragment(
        fragments_file,
        anndata_file,
        generate_index=True,
        dask_cluster_config=dict(processes=False),
    )
    for artifact in expected_artifacts:
        assert artifact.exists()


class TestValidateFragmentBarcodeInAdataIndex:
    def test_postive(self, mock_read_parquet, mock_read_anndata):
        result = atac_seq.validate_fragment_barcode_in_adata_index("fake_file", "fake_file")
        assert not result

    def test_missmatch_anndata(self, mock_read_parquet, mock_read_anndata):
        mock_read_anndata.return_value = ad.AnnData(obs=pd.DataFrame(index=["A", "B", "E"]), var=pd.DataFrame())
        result = atac_seq.validate_fragment_barcode_in_adata_index("fake_file", "fake_file")
        assert result == "Barcodes don't match anndata.obs.index"

    def test_missing_in_anndata(self, mock_read_parquet, mock_read_anndata):
        mock_read_anndata.return_value = ad.AnnData(obs=pd.DataFrame(index=["A", "B"]), var=pd.DataFrame())
        result = atac_seq.validate_fragment_barcode_in_adata_index("fake_file", "fake_file")
        assert result == "Barcodes don't match anndata.obs.index"

    def test_missmatch_in_parquet(self, mock_read_parquet, mock_read_anndata):
        mock_read_parquet.return_value = dd.from_pandas(pd.DataFrame({"barcode": ["A", "B", "E"]}), npartitions=1)
        result = atac_seq.validate_fragment_barcode_in_adata_index("fake_file", "fake_file")
        assert result == "Barcodes don't match anndata.obs.index"

    def test_missing_in_parquet(self, mock_read_parquet, mock_read_anndata):
        mock_read_parquet.return_value = dd.from_pandas(pd.DataFrame({"barcode": ["A", "B"]}), npartitions=1)
        result = atac_seq.validate_fragment_barcode_in_adata_index("fake_file", "fake_file")
        assert result == "Barcodes don't match anndata.obs.index"
