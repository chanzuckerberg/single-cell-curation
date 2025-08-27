import gzip
import os
from pathlib import Path
from unittest import mock

import anndata as ad
import dask.dataframe as dd
import pandas as pd
import pysam
import pytest
from cellxgene_schema import atac_seq
from fixtures.examples_validate import FIXTURES_ROOT

# Test constants
EXPECTED_LINE_COUNT = 37249
TEST_BARCODE = "AAACAAACATTTTATCC-1"


@pytest.fixture
def atac_fragment_bgzip_file_path(tmpdir) -> Path:
    bgzip_file = Path(tmpdir + "new.tsv.bgz")
    yield bgzip_file
    bgzip_file.unlink(missing_ok=True)


@pytest.fixture
def atac_fragment_index_file_path(atac_fragment_bgzip_file_path) -> Path:
    index_file = Path(str(atac_fragment_bgzip_file_path) + ".tbi")
    yield index_file
    index_file.unlink(missing_ok=True)


@pytest.fixture
def test_fragment_files():
    """Returns paths to test fragment files."""
    return {
        "gzip": os.path.join(FIXTURES_ROOT, "atac_seq", "fragments.tsv.gz"),
        "bgzip": os.path.join(FIXTURES_ROOT, "atac_seq", "fragments.tsv.bgz"),
        "anndata": os.path.join(FIXTURES_ROOT, "atac_seq", "small_atac_seq.h5ad"),
    }


@pytest.fixture
def mock_anndata_file(tmpdir):
    """Create a mock anndata file for testing."""
    adata = ad.AnnData(
        obs=pd.DataFrame({"assay_ontology_term_id": ["EFO:0030059"], "is_primary_data": [True]}, index=[TEST_BARCODE]),
        var=pd.DataFrame(index=["GENE1"]),
    )
    adata.uns["organism_ontology_term_id"] = "NCBITaxon:9606"

    anndata_file = os.path.join(tmpdir, "test.h5ad")
    adata.write(anndata_file)
    return anndata_file


def create_fragment_file_from_dataframe(file_path: str, df: pd.DataFrame):
    """Helper to create a gzip fragment file from a DataFrame."""
    with gzip.open(file_path, "wt") as f:
        for _, row in df.iterrows():
            line = "\t".join(map(str, row.tolist())) + "\n"
            f.write(line)
    return file_path


def create_test_fragment_file(tmpdir: Path, filename: str, lines: pd.DataFrame) -> str:
    """Helper to create a gzip fragment file with specified content."""
    file_path = os.path.join(tmpdir, filename)
    with gzip.open(file_path, "wt") as f:
        for line in lines:
            f.write(line)
    return file_path


def to_anndata_file(adata: ad.AnnData, path: str) -> str:
    file_name = os.path.join(path, "small_atac_seq.h5ad")
    adata.write(file_name)
    return file_name


@pytest.fixture
def atac_anndata_file(atac_anndata, tmpdir):
    file_name = os.path.join(tmpdir, "small_atac_seq.h5ad")
    atac_anndata.write(file_name)
    return file_name


def to_parquet_file(df: pd.DataFrame, path: str) -> str:
    file_name = os.path.join(path, "fragment")
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


def count_fragments_per_chromosome(fragment_file):
    fragments = pd.read_csv(
        fragment_file,
        compression="gzip",
        sep="\t",
        header=None,
        names=["chromosome", "start coordinate", "stop coordinate", "barcode", "read support"],
    )
    fragments_per_chromosome = fragments["chromosome"].value_counts().sort_index()
    return fragments_per_chromosome


class TestProcessFragment:
    @pytest.mark.parametrize("fragment_file", ["fragments.tsv.bgz", "fragments.tsv.gz"])
    def test_source_file_compression(self, atac_fragment_bgzip_file_path, atac_fragment_index_file_path, fragment_file):
        self._test_process_fragment(atac_fragment_bgzip_file_path, atac_fragment_index_file_path, fragment_file)

    @staticmethod
    def _test_process_fragment(atac_fragment_bgzip_file_path, atac_fragment_index_file_path, fragment_file):
        # Arrange
        anndata_file = os.path.join(FIXTURES_ROOT, "atac_seq", "small_atac_seq.h5ad")
        fragments_file = os.path.join(FIXTURES_ROOT, "atac_seq", fragment_file)
        # Act
        result = atac_seq.process_fragment(
            str(fragments_file),
            anndata_file,
            generate_index=True,
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

        # check that the number of fragments per chromosome is the same
        orig_chrom_counts = count_fragments_per_chromosome(fragments_file)
        new_chrom_counts = count_fragments_per_chromosome(atac_fragment_bgzip_file_path)
        pd.testing.assert_series_equal(orig_chrom_counts, new_chrom_counts)

        # Testing index access
        assert atac_fragment_index_file_path.exists()
        with pysam.TabixFile(str(atac_fragment_bgzip_file_path)) as tabix:
            for chromosome in tabix.contigs:
                assert chromosome in atac_seq.human_chromosome_by_length
                assert orig_chrom_counts[chromosome] == len(list(tabix.fetch(chromosome)))

    def test_fail(self, atac_fragment_bgzip_file_path, atac_fragment_index_file_path):
        # Arrange
        anndata_file = os.path.join(FIXTURES_ROOT, "atac_seq", "small_atac_seq.h5ad")
        fragments_file = os.path.join(FIXTURES_ROOT, "atac_seq", "fragments_bad.tsv.gz")
        # Act
        result = atac_seq.process_fragment(
            fragments_file,
            anndata_file,
            generate_index=True,
            output_file=str(atac_fragment_bgzip_file_path),
        )
        result = [r for r in result if "Error" in r]
        # Assert
        assert len(result) == 1


class TestPrepareFragment:
    def test_positive(self, atac_fragment_dataframe, tmpdir):
        # Arrange
        input_file = os.path.join(
            tmpdir,
            "fragment.tsv.gz",
        )
        atac_fragment_dataframe.to_csv(
            input_file, sep="\t", index=False, compression="gzip", header=False, columns=atac_seq.column_ordering
        )
        output_file = os.path.join(tmpdir, "fragment.bgz")
        # Act
        atac_seq.prepare_fragment(input_file, output_file)
        # Assert
        assert Path(output_file).exists()

    @pytest.mark.parametrize("tool", ["bgzip", "pigz", "sort"])
    @mock.patch("shutil.which")
    def test_missing_requirements(self, mock, tool, atac_fragment_file, tmpdir):
        # Arrange
        mock.return_value = None  # patch shutil to return None for the tool
        # Act
        with pytest.raises(RuntimeError):
            atac_seq.prepare_fragment(atac_fragment_file, "fragment.bgz")


class TestConvertToParquet:
    def test_positive(self, atac_fragment_dataframe, tmpdir):
        tsv_file = os.path.join(
            tmpdir,
            "fragment.tsv.gz",
        )
        atac_fragment_dataframe.to_csv(
            tsv_file, sep="\t", index=False, compression="gzip", header=False, columns=atac_seq.column_ordering
        )
        parquet_file = Path(atac_seq.convert_to_parquet(tsv_file, tmpdir))
        assert Path(parquet_file).is_dir()

    def test_missing_column(self, tmpdir, atac_fragment_dataframe):
        atac_fragment_dataframe = atac_fragment_dataframe.drop(columns=["read support"])
        tsv_file = os.path.join(tmpdir, "fragment.tsv")
        atac_fragment_dataframe.to_csv(
            tsv_file, sep="\t", index=False, compression="gzip", header=False, columns=atac_seq.column_ordering[:-1]
        )
        with pytest.raises(ValueError):
            atac_seq.convert_to_parquet(tsv_file, tmpdir)

    def test_invalid_column_dtype(self, tmpdir, atac_fragment_dataframe):
        atac_fragment_dataframe["start coordinate"] = "foo"
        tsv_file = os.path.join(tmpdir, "fragment.tsv.gz")
        atac_fragment_dataframe.to_csv(
            tsv_file, sep="\t", index=False, compression="gzip", header=False, columns=atac_seq.column_ordering
        )
        with pytest.raises(ValueError):
            atac_seq.convert_to_parquet(tsv_file, tmpdir)

    def test_with_na_columns(self, atac_fragment_dataframe, tmpdir):
        atac_fragment_dataframe["barcode"] = pd.NA
        tsv_file = os.path.join(tmpdir, "fragment.tsv.gz")
        atac_fragment_dataframe.to_csv(
            tsv_file, sep="\t", index=False, compression="gzip", header=False, columns=atac_seq.column_ordering
        )
        parquet_file = Path(atac_seq.convert_to_parquet(tsv_file, tmpdir))
        parquet_df = dd.read_parquet(parquet_file, columns=["barcode"])
        assert len(parquet_df[parquet_df["barcode"] != ""].compute()) == 0


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
    def test_positive(self, atac_fragment_file):
        result = atac_seq.validate_fragment_stop_coordinate_within_chromosome(atac_fragment_file, "NCBITaxon:9606")
        assert not result

    def test_stop_less_than_chromosome_length(self, atac_fragment_dataframe, tmpdir):
        # Arrange
        atac_fragment_dataframe["stop coordinate"] = 10e12
        fragment_file = to_parquet_file(atac_fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_stop_coordinate_within_chromosome(fragment_file, "NCBITaxon:9606")
        # Assert
        assert result == "Stop coordinate must be less than the chromosome length."

    def test_mismatch_chromosome(self, atac_fragment_dataframe, tmpdir):
        # Arrange
        atac_fragment_dataframe["chromosome"] = ["foo", "chr2", "chr1"]
        fragment_file = to_parquet_file(atac_fragment_dataframe, tmpdir)
        # Act
        result = atac_seq.validate_fragment_stop_coordinate_within_chromosome(fragment_file, "NCBITaxon:9606")
        # Assert
        assert result.startswith("Chromosomes in the fragment do not match the organism")


class TestValidateFragmentReadSupport:
    def test_positive(self, atac_fragment_file):
        result = atac_seq.validate_fragment_read_support(atac_fragment_file)
        assert not result

    @pytest.mark.parametrize("read_support", [0, -1])
    def test_negative(self, atac_fragment_dataframe, tmpdir, read_support):
        # Arrange
        atac_fragment_dataframe["read support"] = read_support
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
        atac_anndata.uns["organism_ontology_term_id"] = "NCBITaxon:9606"
        atac_anndata_file = to_anndata_file(atac_anndata, tmpdir)
        # Act
        result = atac_seq.validate_anndata_organism_ontology_term_id(atac_anndata_file)
        # Assert
        assert not result

    def test_organism_ontology_term_id_not_allowed(self, atac_anndata, tmpdir):
        # Arrange
        atac_anndata.uns["organism_ontology_term_id"] = "NCBITaxon:9607"
        atac_anndata_file = to_anndata_file(atac_anndata, tmpdir)
        # Act
        result = atac_seq.validate_anndata_organism_ontology_term_id(atac_anndata_file)
        # Assert
        assert result.startswith(
            "Anndata.obs.organism_ontology_term_id must be one of ['NCBITaxon:9606', 'NCBITaxon:10090']."
        )


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
        assert result.startswith("Fragment file has duplicate rows.")
        for chrom in atac_fragment_dataframe["chromosome"].unique():
            assert f"Chromosome {chrom} has 2 rows but only 1 are unique" in result


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


class TestCommandGeneration:
    """Test command generation functions for external tools."""

    @pytest.mark.parametrize(
        "func,args,expected_items",
        [
            (atac_seq.get_sort_command, (4,), ["sort", "-S20%", "--parallel=4", "-k1,1", "-k2,2n"]),
            (atac_seq.get_sort_command, (8, 60), ["-S7%", "--parallel=8"]),  # custom memory
            (atac_seq.get_bgzip_command, (6,), ["bgzip", "--threads=6", "-c"]),
            (atac_seq.get_gzip_command, ("/path/to/file.tsv.gz",), ["gzip", "-dc", "/path/to/file.tsv.gz"]),
        ],
    )
    def test_command_generation(self, func, args, expected_items):
        """Test command generation functions return correct parameters."""
        result = func(*args)
        for item in expected_items:
            assert item in result

    def test_get_sort_command_full_structure(self):
        """Test complete sort command structure for critical sorting parameters."""
        command = atac_seq.get_sort_command(4)
        expected_command = [
            "sort",
            "-t",
            "\t",
            "-k1,1",
            "-k2,2n",
            "-k3,3n",
            "-k4,4",
            "-S20%",
            "--parallel=4",
            '--compress-program="pigz -p 4"',
        ]
        assert command == expected_command


class TestCalculateSortMemory:
    """Test memory calculation for sort command."""

    @pytest.mark.parametrize(
        "cores,memory,expected",
        [
            (1, 80, 80),  # single core
            (2, 80, 40),  # basic case
            (4, 80, 20),  # common case
            (8, 80, 10),  # multi-core
            (16, 80, 5),  # many cores
            (32, 80, 2),  # more cores
            (64, 80, 1),  # minimum boundary
            (100, 80, 1),  # excessive cores (minimum case)
            (16, 10, 1),  # edge case - low memory percentage
        ],
    )
    def test_calculate_sort_memory(self, cores, memory, expected):
        """Test memory calculation with various scenarios including edge cases."""
        assert atac_seq.calculate_sort_memory(cores, memory) == expected


class TestCheckExternalRequirements:
    """Test external tool requirements checking."""

    def test_check_external_requirements_all_present(self):
        """Test that check passes when all tools are present."""
        with mock.patch("shutil.which", return_value="/usr/bin/tool"):
            atac_seq.check_external_requirements()  # Should not raise

    @pytest.mark.parametrize(
        "missing_tool,expected_error",
        [
            ("sort", "The 'sort' command is not installed or not found in PATH"),
            ("bgzip", "The 'bgzip' command is not installed or not found in PATH"),
            ("pigz", "The 'pigz' command is not installed or not found in PATH"),
            ("uniq", "The 'uniq' command is not installed or not found in PATH"),
        ],
    )
    def test_check_external_requirements_missing_tool(self, missing_tool, expected_error):
        """Test that check raises RuntimeError when specific tools are missing."""

        def mock_which(tool):
            return None if tool == missing_tool else "/usr/bin/tool"

        with mock.patch("shutil.which", side_effect=mock_which), pytest.raises(RuntimeError, match=expected_error):
            atac_seq.check_external_requirements()

    def test_check_external_requirements_all_missing(self):
        """Test that check raises error for first missing tool when all are missing."""
        with (
            mock.patch("shutil.which", return_value=None),
            pytest.raises(RuntimeError, match="The 'sort' command is not installed"),
        ):
            atac_seq.check_external_requirements()


class TestCountLinesInCompressedFile:
    """Test line counting in compressed files."""

    @pytest.mark.parametrize(
        "file_type,expected_count",
        [
            ("gzip", EXPECTED_LINE_COUNT),
            ("bgzip", EXPECTED_LINE_COUNT),
        ],
    )
    def test_count_lines_compressed_files(self, test_fragment_files, file_type, expected_count):
        """Test counting lines in gzip and bgzip compressed fragment files."""
        line_count = atac_seq.count_lines_in_compressed_file(test_fragment_files[file_type])
        assert line_count == expected_count

    @pytest.mark.parametrize(
        "test_data,expected_count",
        [
            (pd.DataFrame(columns=atac_seq.column_ordering), 0),  # empty file
            (pd.DataFrame([["chr1", 100, 200, "barcode1", 5]], columns=atac_seq.column_ordering), 1),  # single line
        ],
    )
    def test_count_lines_various_sizes(self, tmpdir, test_data, expected_count):
        """Test counting lines in files of various sizes."""
        test_file = create_fragment_file_from_dataframe(os.path.join(tmpdir, "test.tsv.gz"), test_data)
        line_count = atac_seq.count_lines_in_compressed_file(test_file)
        assert line_count == expected_count

    def test_count_lines_nonexistent_file(self):
        """Test handling of non-existent file."""
        with pytest.raises(FileNotFoundError):
            atac_seq.count_lines_in_compressed_file("/nonexistent/file.gz")


class TestIndexFragmentPrepared:
    """Test fragment indexing with prepared/unprepared files."""

    def test_index_fragment_already_prepared(self, test_fragment_files, tmpdir):
        """Test index_fragment with is_prepared=False (file already prepared)."""
        # Arrange
        fragment_file = test_fragment_files["bgzip"]  # Use already prepared bgzip file

        # Act
        atac_seq.index_fragment(
            fragment_file=fragment_file,
            output_file=None,
            is_prepared=False,  # File is already prepared
        )

        # Assert
        index_file = fragment_file + ".tbi"
        assert Path(index_file).exists()

        # Clean up
        Path(index_file).unlink(missing_ok=True)

    def test_index_fragment_needs_preparation(self, tmpdir, mock_anndata_file):
        """Test index_fragment with is_prepared=True (needs preparation)."""
        # Arrange
        test_lines = [f"chr1\t100\t200\t{TEST_BARCODE}\t5\n", f"chr1\t300\t400\t{TEST_BARCODE}\t3\n"]
        test_fragment_file = create_fragment_file_from_dataframe(
            os.path.join(tmpdir, "test_fragments.tsv.gz"),
            pd.DataFrame([line.rstrip("\n").split("\t") for line in test_lines], columns=atac_seq.column_ordering),
        )
        output_file = os.path.join(tmpdir, "output.bgz")
        index_file = output_file + ".tbi"

        # Mock prepare_fragment to create a proper bgzip file
        def mock_prepare_fragment(input_file, output_file):
            with gzip.open(input_file, "rt") as f_in:
                data = f_in.read()
            with gzip.open(output_file, "wt") as f_out:
                f_out.write(data)

        # Mock line_counts_match and tabix_index
        def mock_line_counts_match(file1, file2):
            pass  # Do nothing for test

        def mock_tabix_index(bgzip_file, preset=None, force=None):
            Path(bgzip_file + ".tbi").touch()

        # Act
        with (
            mock.patch("cellxgene_schema.atac_seq.prepare_fragment", side_effect=mock_prepare_fragment),
            mock.patch("cellxgene_schema.atac_seq.line_counts_match", side_effect=mock_line_counts_match),
            mock.patch("pysam.tabix_index", side_effect=mock_tabix_index),
        ):
            atac_seq.index_fragment(
                fragment_file=test_fragment_file,
                output_file=output_file,
                is_prepared=True,  # File needs preparation
            )

        # Assert
        assert Path(output_file).exists()
        assert Path(index_file).exists()

    def test_index_fragment_parameter_logic(self):
        """Test that index_fragment is_prepared parameter logic works correctly."""
        with (
            mock.patch("cellxgene_schema.atac_seq.prepare_fragment") as mock_prepare,
            mock.patch("cellxgene_schema.atac_seq.line_counts_match") as mock_line_counts,
            mock.patch("pysam.tabix_index") as mock_tabix,
        ):
            # Test is_prepared=True (should call prepare_fragment)
            atac_seq.index_fragment("test.gz", "output.bgz", is_prepared=True)
            mock_prepare.assert_called_once_with("test.gz", "output.bgz")
            mock_line_counts.assert_called_once_with("test.gz", "output.bgz")
            mock_tabix.assert_called_once()

            # Reset mocks
            mock_prepare.reset_mock()
            mock_line_counts.reset_mock()
            mock_tabix.reset_mock()

            # Test is_prepared=False (should NOT call prepare_fragment)
            atac_seq.index_fragment("test.bgz", None, is_prepared=False)
            mock_prepare.assert_not_called()
            mock_line_counts.assert_not_called()
            mock_tabix.assert_called_once_with("test.bgz", preset="bed", force=True)


class TestLineCountsMatch:
    """Test line count validation between files."""

    @pytest.mark.parametrize(
        "data1,data2,should_match",
        [
            # Both empty
            (pd.DataFrame(columns=atac_seq.column_ordering), pd.DataFrame(columns=atac_seq.column_ordering), True),
            # Same data
            (
                pd.DataFrame(
                    {
                        "chromosome": ["chr1", "chr2"],
                        "start coordinate": [100, 200],
                        "stop coordinate": [200, 300],
                        "barcode": ["A", "B"],
                        "read support": [1, 2],
                    }
                ),
                pd.DataFrame(
                    {
                        "chromosome": ["chr1", "chr2"],
                        "start coordinate": [100, 200],
                        "stop coordinate": [200, 300],
                        "barcode": ["A", "B"],
                        "read support": [1, 2],
                    }
                ),
                True,
            ),
            # Different lengths
            (
                pd.DataFrame(
                    {
                        "chromosome": ["chr1", "chr2"],
                        "start coordinate": [100, 200],
                        "stop coordinate": [200, 300],
                        "barcode": ["A", "B"],
                        "read support": [1, 2],
                    }
                ),
                pd.DataFrame(
                    {
                        "chromosome": ["chr1"],
                        "start coordinate": [100],
                        "stop coordinate": [200],
                        "barcode": ["A"],
                        "read support": [1],
                    }
                ),
                False,
            ),
        ],
    )
    def test_line_counts_match(self, tmpdir, data1, data2, should_match):
        """Test line count matching with various data scenarios."""
        file1 = create_fragment_file_from_dataframe(os.path.join(tmpdir, "file1.tsv.gz"), data1)
        file2 = create_fragment_file_from_dataframe(os.path.join(tmpdir, "file2.tsv.gz"), data2)

        if should_match:
            atac_seq.line_counts_match(file1, file2)  # Should not raise
        else:
            with pytest.raises(ValueError, match="Line count validation failed"):
                atac_seq.line_counts_match(file1, file2)


class TestDeduplicateFragmentRows:
    def test_deduplicate_rows(self, atac_fragment_dataframe, tmpdir):
        # Arrange
        atac_fragment_dataframe = pd.concat([atac_fragment_dataframe, atac_fragment_dataframe])
        input_file = create_fragment_file_from_dataframe(
            os.path.join(tmpdir, "fragment.tsv.gz"), atac_fragment_dataframe
        )
        output_file = os.path.join(tmpdir, "deduplicated.tsv.bgz")
        # Act
        atac_seq.deduplicate_fragment_rows(input_file, output_file)
        # Assert
        assert Path(output_file).exists()
        df = pd.read_csv(
            output_file,
            compression="gzip",
            sep="\t",
            header=None,
            names=["chromosome", "start coordinate", "stop coordinate", "barcode", "read support"],
        )
        assert len(df) == len(atac_fragment_dataframe) // 2

    def test_no_duplicates(self, atac_fragment_dataframe, tmpdir):
        # Arrange
        input_file = create_fragment_file_from_dataframe(
            os.path.join(tmpdir, "fragment.tsv.gz"), atac_fragment_dataframe
        )
        output_file = os.path.join(tmpdir, "deduplicated.tsv.bgz")
        # Act
        atac_seq.deduplicate_fragment_rows(input_file, output_file)
        # Assert
        assert Path(output_file).exists()
        df = pd.read_csv(
            output_file,
            compression="gzip",
            sep="\t",
            header=None,
            names=["chromosome", "start coordinate", "stop coordinate", "barcode", "read support"],
        )
        assert len(df) == len(atac_fragment_dataframe)
