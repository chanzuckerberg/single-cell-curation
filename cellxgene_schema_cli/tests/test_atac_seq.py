from pathlib import Path

import anndata as ad
import dask.dataframe as dd
import pandas as pd
import pysam
import pytest
from cellxgene_schema import atac_seq
from fixtures.examples_validate import FIXTURES_ROOT


@pytest.fixture
def atac_fragment_bgzip_file_path() -> Path:
    bgzip_file = Path(FIXTURES_ROOT + "/atac_seq/new.tsv.bgz")
    yield bgzip_file
    bgzip_file.unlink(missing_ok=True)


@pytest.fixture
def atac_fragment_index_file_path(atac_fragment_bgzip_file_path) -> Path:
    index_file = Path(str(atac_fragment_bgzip_file_path) + ".tbi")
    yield index_file
    index_file.unlink(missing_ok=True)


def to_anndata_file(adata: ad.AnnData, path: str) -> str:
    file_name = path + "/small_atac_seq.h5ad"
    adata.write(file_name)
    return file_name


@pytest.fixture
def atac_anndata():
    obs = pd.DataFrame(index=["A", "B", "C"])
    obs["organism_ontology_term_id"] = ["NCBITaxon:9606"] * 3
    obs["assay_ontology_term_id"] = ["EFO:0030059"] * 3
    return ad.AnnData(obs=obs, var=pd.DataFrame())


@pytest.fixture
def atac_anndata_file(atac_anndata, tmpdir):
    file_name = tmpdir + "/small_atac_seq.h5ad"
    atac_anndata.write(file_name)
    return file_name


def to_parquet_file(df: pd.DataFrame, path: str) -> str:
    file_name = path + "/fragment"
    dd.from_pandas(df).to_parquet(file_name, partition_on=["chromosome"])
    return file_name


@pytest.fixture
def atac_fragment_dataframe() -> pd.DataFrame:
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
def atac_fragment_file(atac_fragment_dataframe, tmpdir):
    return to_parquet_file(atac_fragment_dataframe, tmpdir)


class TestProcessFragment:
    @pytest.mark.parametrize("override_write_algorithm", ["pysam", "cli", None])
    def test_write_algorithms(
        self, atac_fragment_bgzip_file_path, atac_fragment_index_file_path, override_write_algorithm
    ):
        self._test_process_fragment(
            atac_fragment_bgzip_file_path, atac_fragment_index_file_path, override_write_algorithm, "fragments.tsv.bgz"
        )

    @pytest.mark.parametrize("fragment_file", ["fragments.tsv.bgz", "fragments_sorted.tsv.gz"])
    def test_source_file_compression(self, atac_fragment_bgzip_file_path, atac_fragment_index_file_path, fragment_file):
        self._test_process_fragment(
            atac_fragment_bgzip_file_path, atac_fragment_index_file_path, "pysam", fragment_file
        )

    @staticmethod
    def _test_process_fragment(
        atac_fragment_bgzip_file_path, atac_fragment_index_file_path, override_write_algorithm, fragment_file
    ):
        # Arrange
        anndata_file = FIXTURES_ROOT + "/atac_seq/small_atac_seq.h5ad"
        fragments_file = FIXTURES_ROOT + f"/atac_seq/{fragment_file}"
        # Act
        result = atac_seq.process_fragment(
            fragments_file,
            anndata_file,
            generate_index=True,
            dask_cluster_config=dict(processes=False),
            override_write_algorithm=override_write_algorithm,
            output_file=str(atac_fragment_bgzip_file_path),
        )
        # Assert
        assert len(result) == 0
        assert atac_fragment_bgzip_file_path.exists()
        # Testing the bgzip file
        with pysam.libcbgzf.BGZFile(atac_fragment_bgzip_file_path, mode="r") as bgzip:
            previous_stop, previous_start, previous_chomosome = 0, 0, None
            for line in bgzip:
                chromosome, start, stop, _, _ = line.decode("utf-8").split("\t")
                if previous_chomosome != chromosome:
                    previous_stop, previous_start = 0, 0
                    chromosome_length = atac_seq.human_chromosome_by_length[chromosome]
                start, stop = int(start), int(stop)
                assert start <= stop < chromosome_length, (
                    "Stop coordinate must be within the chromosome length and " "greater than start."
                )
                assert start >= 0, "Start coordinate must be greater than 0."
                # Fragment is sorted by start, then stop in ascending order
                assert start >= previous_start
                # The bellow check is failing.
                if (chromosome == previous_chomosome) and (start == previous_start):
                    assert stop >= previous_stop
                previous_stop = stop
                previous_start = start
                previous_chomosome = chromosome

        # Testing index access
        assert atac_fragment_index_file_path.exists()
        with pysam.TabixFile(str(atac_fragment_bgzip_file_path)) as tabix:
            for chromosome in tabix.contigs:
                assert chromosome in atac_seq.human_chromosome_by_length


class TestValidateFragmentBarcodeInAdataIndex:
    def test_postive(self, atac_fragment_file, atac_anndata_file):
        result = atac_seq.validate_fragment_barcode_in_adata_index(atac_fragment_file, atac_anndata_file)
        assert not result

    def test_missmatch_anndata(self, atac_fragment_file, tmpdir):
        # Arrange
        atac_anndata_file = to_anndata_file(
            ad.AnnData(obs=pd.DataFrame(index=["A", "B", "E"]), var=pd.DataFrame()), tmpdir
        )
        # Act
        result = atac_seq.validate_fragment_barcode_in_adata_index(atac_fragment_file, atac_anndata_file)
        # Assert
        assert result == "Barcodes don't match anndata.obs.index"

    def test_missing_in_anndata(self, atac_fragment_file, atac_anndata_file, tmpdir):
        # Arrange
        atac_anndata_file = to_anndata_file(ad.AnnData(obs=pd.DataFrame(index=["A", "B"]), var=pd.DataFrame()), tmpdir)
        # Act
        result = atac_seq.validate_fragment_barcode_in_adata_index(atac_fragment_file, atac_anndata_file)
        # Assert
        assert result == "Barcodes don't match anndata.obs.index"

    def test_missmatch_in_parquet(self, atac_fragment_dataframe, atac_anndata_file, tmpdir):
        # Arrange
        atac_fragment_dataframe["barcode"] = ["A", "B", "E"]
        fragment_file = to_parquet_file(atac_fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_barcode_in_adata_index(fragment_file, atac_anndata_file)
        # Assert
        assert result == "Barcodes don't match anndata.obs.index"

    def test_missing_in_parquet(self, atac_fragment_dataframe, atac_anndata_file, tmpdir):
        # Arrange
        fragment_dataframe = atac_fragment_dataframe.drop(index=2)
        fragment_file = to_parquet_file(fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_barcode_in_adata_index(fragment_file, atac_anndata_file)
        # Assert
        assert result == "Barcodes don't match anndata.obs.index"


class TestValidateFragmentStartCoordianteGreaterThan0:
    def test_positive(self, atac_fragment_file):
        result = atac_seq.validate_fragment_start_coordinate_greater_than_0(atac_fragment_file)
        assert not result

    def test_negative(self, atac_fragment_dataframe, tmpdir):
        # Arrange
        atac_fragment_dataframe["start coordinate"] = -1
        fragment_file = to_parquet_file(atac_fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_start_coordinate_greater_than_0(fragment_file)
        # Assert
        assert result == "Start coordinate must be greater than 0."


class TestValidateFragmentStopCoordinateGreaterThanStartCoordinate:
    def test_positive(self, atac_fragment_file):
        result = atac_seq.validate_fragment_stop_greater_than_start_coordinate(atac_fragment_file)
        assert not result

    def test_negative(self, atac_fragment_dataframe, tmpdir):
        # Arrange
        atac_fragment_dataframe["stop coordinate"] = 1
        fragment_file = to_parquet_file(atac_fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_stop_greater_than_start_coordinate(fragment_file)
        # Assert
        assert result == "Stop coordinate must be greater than start coordinate."


class TestValidateFragmentStopCoordinateWithinChromosome:
    def test_positive(self, atac_fragment_file, atac_anndata_file):
        result = atac_seq.validate_fragment_stop_coordinate_within_chromosome(atac_fragment_file, atac_anndata_file)
        assert not result

    def test_stop_less_than_chromosome_length(self, atac_fragment_dataframe, atac_anndata_file, tmpdir):
        # Arrange
        atac_fragment_dataframe["stop coordinate"] = 10e12
        fragment_file = to_parquet_file(atac_fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_stop_coordinate_within_chromosome(fragment_file, atac_anndata_file)
        # Assert
        assert result == "Stop coordinate must be less than the chromosome length."

    def test_organism_ontology_id_not_unique(self, atac_fragment_dataframe, atac_anndata_file, tmpdir):
        # Arrange
        atac_anndata = ad.read_h5ad(atac_anndata_file)
        atac_anndata.obs["organism_ontology_term_id"] = ["NCBITaxon:10090", "NCBITaxon:9606", "NCBITaxon:9606"]
        atac_anndata.write(tmpdir + "/small_atac_seq.h5ad")
        # Act
        result = atac_seq.validate_fragment_stop_coordinate_within_chromosome(
            to_parquet_file(atac_fragment_dataframe, tmpdir), tmpdir + "/small_atac_seq.h5ad"
        )
        # Assert
        assert result == "Anndata.obs.organism_ontology_term_id must have a unique value."


class TestValidateFragmentReadSupport:
    def test_positive(self, atac_fragment_file):
        result = atac_seq.validate_fragment_read_support(atac_fragment_file)
        assert not result

    def test_negative(self, atac_fragment_dataframe, tmpdir):
        # Arrange
        atac_fragment_dataframe["read support"] = 0
        fragment_file = to_parquet_file(atac_fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_read_support(fragment_file)
        # Assert
        assert result == "Read support must be greater than 0."


class TestCheckAnndataRequiresFragment:

    @pytest.mark.parametrize("assay_ontology_term_id,expected_result", [("EFO:0030059", False), ("EFO:0030007", True)])
    def test_positive(self, atac_anndata, assay_ontology_term_id, expected_result, tmpdir):
        # Arrange
        atac_anndata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        atac_anndata_file = to_anndata_file(atac_anndata, tmpdir)
        # Act
        result = atac_seq.check_anndata_requires_fragment(atac_anndata_file)
        # Assert
        assert result == expected_result

    def test_not_atac(self, atac_anndata, tmpdir):
        # Arrange
        atac_anndata.obs["assay_ontology_term_id"] = "EFO:0030060"
        atac_anndata_file = to_anndata_file(atac_anndata, tmpdir)
        with pytest.raises(
            ValueError, match="Anndata.obs.assay_ontology_term_id are not all descendants of EFO:0010891."  # Assert
        ):
            # Act
            atac_seq.check_anndata_requires_fragment(atac_anndata_file)

    def test_mixed_paired_and_unpaired(self, atac_anndata, tmpdir):
        # Arrange
        atac_anndata.obs["assay_ontology_term_id"] = ["EFO:0030059", "EFO:0030007", "EFO:0030007"]
        atac_anndata_file = to_anndata_file(atac_anndata, tmpdir)
        with pytest.raises(
            ValueError, match="Anndata.obs.assay_ontology_term_id has mixed paired and unpaired assay terms."
        ):  # Assert
            # Act
            atac_seq.check_anndata_requires_fragment(atac_anndata_file)


class TestValidateAnndataIsPrimaryData:
    def test_positive(self, atac_anndata, tmpdir):
        # Arrange
        atac_anndata.obs["is_primary_data"] = True
        atac_anndata_file = to_anndata_file(atac_anndata, tmpdir)
        # Act
        result = atac_seq.validate_anndata_is_primary_data(atac_anndata_file)
        # Assert
        assert not result

    def test_negative(self, atac_anndata, tmpdir):
        # Arrange
        atac_anndata.obs["is_primary_data"] = False
        atac_anndata_file = to_anndata_file(atac_anndata, tmpdir)
        # Act
        result = atac_seq.validate_anndata_is_primary_data(atac_anndata_file)
        # Assert
        assert result == "Anndata.obs.is_primary_data must all be True."


class TestValidateAnndataOrganismOntologyTermId:
    def test_positive(self, atac_anndata, tmpdir):
        # Arrange
        atac_anndata.obs["organism_ontology_term_id"] = ["NCBITaxon:9606"] * 3
        atac_anndata_file = to_anndata_file(atac_anndata, tmpdir)
        # Act
        result = atac_seq.validate_anndata_organism_ontology_term_id(atac_anndata_file)
        # Assert
        assert not result

    def test_organism_ontology_term_id_not_allowed(self, atac_anndata, tmpdir):
        # Arrange
        atac_anndata.obs["organism_ontology_term_id"] = ["NCBITaxon:9607"] * 3
        atac_anndata_file = to_anndata_file(atac_anndata, tmpdir)
        # Act
        result = atac_seq.validate_anndata_organism_ontology_term_id(atac_anndata_file)
        # Assert
        assert result == "Anndata.obs.organism_ontology_term_id must be one of ['NCBITaxon:9606', 'NCBITaxon:10090']."


class TestValidateFragmentNoDuplicateRows:
    def test_positive(self, atac_fragment_file):
        result = atac_seq.validate_fragment_no_duplicate_rows(atac_fragment_file)
        assert not result

    def test_negative(self, atac_fragment_dataframe, tmpdir):
        # Arrange
        atac_fragment_dataframe = pd.concat([atac_fragment_dataframe, atac_fragment_dataframe])
        fragment_file = to_parquet_file(atac_fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_no_duplicate_rows(fragment_file)
        # Assert
        assert result == "Fragment file has duplicate rows."


class TestGetOutputFile:
    @pytest.mark.parametrize(
        "fragment_file,output_file,expected",
        [
            ("fragment.gz", None, "fragment.bgz"),
            ("fragment.gz", "output.bgz", "output.bgz"),
            ("fragment.gz", "output", "output.bgz"),
        ],
    )
    def test_none(self, fragment_file, output_file, expected):
        assert atac_seq.get_output_file(fragment_file, output_file) == expected
