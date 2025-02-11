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
            "start coordinate": [100, 200, 300],
            "stop coordinate": [200, 300, 400],
            "read support": [1, 2, 3],
            "chromosome": ["chr1", "chr2", "chr3"],
        }
    )


@pytest.fixture
def mock_read_anndata():
    with mock.patch("anndata.read_h5ad") as mock_anndata:
        obs = pd.DataFrame(index=["A", "B", "C"], columns=["organism_ontology_term_id"], data=["NCBITaxon:9606"] * 3)
        mock_anndata.return_value = ad.AnnData(obs=obs, var=pd.DataFrame())
        yield mock_anndata


def to_parquet(df: pd.DataFrame, tmpdir: str) -> str:
    path = tmpdir + "/fragment"
    dd.from_pandas(df).to_parquet(path, partition_on=["chromosome"])
    return path


@pytest.fixture
def fragment_file(fragment_dataframe, tmpdir):
    return to_parquet(fragment_dataframe, tmpdir)


@pytest.fixture
def mock_read_parquet(fragment_dataframe):
    with mock.patch("dask.dataframe.read_parquet") as mock_read_parquet:
        mock_read_parquet.return_value = dd.from_pandas(fragment_dataframe, npartitions=1)
        yield mock_read_parquet


def test_process_fragment(expected_artifacts):
    # Arrange
    anndata_file = FIXTURES_ROOT + "/atac_seq/small_atac_seq.h5ad"
    fragments_file = FIXTURES_ROOT + "/atac_seq/fragments_sorted.tsv.gz"
    # Act
    atac_seq.process_fragment(
        fragments_file,
        anndata_file,
        generate_index=True,
        dask_cluster_config=dict(processes=False),
    )
    # Assert
    for artifact in expected_artifacts:
        assert artifact.exists()


class TestValidateFragmentBarcodeInAdataIndex:
    def test_postive(self, fragment_file, mock_read_anndata):
        result = atac_seq.validate_fragment_barcode_in_adata_index(fragment_file, "fake_file")
        assert not result

    def test_missmatch_anndata(self, fragment_file, mock_read_anndata, tmpdir):
        # Arrange
        mock_read_anndata.return_value = ad.AnnData(obs=pd.DataFrame(index=["A", "B", "E"]), var=pd.DataFrame())
        # Act
        result = atac_seq.validate_fragment_barcode_in_adata_index(fragment_file, "fake_file")
        # Assert
        assert result == "Barcodes don't match anndata.obs.index"

    def test_missing_in_anndata(self, fragment_file, mock_read_anndata):
        # Arrange
        mock_read_anndata.return_value = ad.AnnData(obs=pd.DataFrame(index=["A", "B"]), var=pd.DataFrame())
        # Act
        result = atac_seq.validate_fragment_barcode_in_adata_index(fragment_file, "fake_file")
        # Assert
        assert result == "Barcodes don't match anndata.obs.index"

    def test_missmatch_in_parquet(self, fragment_dataframe, mock_read_anndata, tmpdir):
        # Arrange
        fragment_dataframe["barcode"] = ["A", "B", "E"]
        fragment_file = to_parquet(fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_barcode_in_adata_index(fragment_file, "fake_file")
        # Assert
        assert result == "Barcodes don't match anndata.obs.index"

    def test_missing_in_parquet(self, fragment_dataframe, mock_read_anndata, tmpdir):
        # Arrange
        fragment_dataframe = fragment_dataframe.drop(index=2)
        fragment_file = to_parquet(fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_barcode_in_adata_index(fragment_file, "fake_file")
        # Assert
        assert result == "Barcodes don't match anndata.obs.index"


class TestValidateFragmentStartCoordianteGreaterThan0:
    def test_positive(self, fragment_file):
        result = atac_seq.validate_fragment_start_coordinate_greater_than_0(fragment_file)
        assert not result

    def test_negative(self, fragment_dataframe, tmpdir):
        # Arrange
        fragment_dataframe["start coordinate"] = -1
        fragment_file = to_parquet(fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_start_coordinate_greater_than_0(fragment_file)
        # Assert
        assert result == "Start coordinate must be greater than 0."


class TestValidateFragmentStopCoordinateGreaterThanStartCoordinate:
    def test_positive(self, fragment_file):
        result = atac_seq.validate_fragment_stop_greater_than_start_coordinate(fragment_file)
        assert not result

    def test_negative(self, fragment_dataframe, tmpdir):
        # Arrange
        fragment_dataframe["stop coordinate"] = 1
        fragment_file = to_parquet(fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_stop_greater_than_start_coordinate(fragment_file)
        # Assert
        assert result == "Stop coordinate must be greater than start coordinate."


class TestValidateFragmentStopCoordinateWithinChromosome:
    def test_positive(self, fragment_file, mock_read_anndata):
        result = atac_seq.validate_fragment_stop_coordinate_within_chromosome(fragment_file, "fake_file")
        assert not result

    def test_negative(self, fragment_dataframe, mock_read_anndata, tmpdir):
        # Arrange
        fragment_dataframe["stop coordinate"] = 10e12
        fragment_file = to_parquet(fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_stop_coordinate_within_chromosome(fragment_file, "fake_file")
        # Assert
        assert result == "Stop coordinate must be less than the chromosome length."


class TestValidateFragmentReadSupport:
    def test_positive(self, fragment_file):
        result = atac_seq.validate_fragment_read_support(fragment_file)
        assert not result

    def test_negative(self, fragment_dataframe, tmpdir):
        # Arrange
        fragment_dataframe["read support"] = 0
        fragment_file = to_parquet(fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_read_support(fragment_file)
        # Assert
        assert result == "Read support must be greater than 0."
