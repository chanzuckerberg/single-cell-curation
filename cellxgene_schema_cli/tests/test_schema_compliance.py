"""
Tests for schema compliance of an AnnData object
"""

import unittest
from copy import deepcopy

import anndata
import fixtures.examples_validate as examples
import numpy
import pandas as pd
import pytest
import scipy.sparse
from cellxgene_schema.schema import get_schema_definition
from cellxgene_schema.utils import getattr_anndata, read_h5ad
from cellxgene_schema.validate import (
    ASSAY_VISIUM_11M,
    ERROR_SUFFIX_IS_SINGLE,
    ERROR_SUFFIX_SPATIAL,
    ERROR_SUFFIX_VISIUM,
    ERROR_SUFFIX_VISIUM_11M,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_IN_TISSUE_0,
    SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE,
    SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE_VISIUM_11MM,
    VISIUM_11MM_AND_IS_SINGLE_TRUE_MATRIX_SIZE,
    VISIUM_AND_IS_SINGLE_TRUE_MATRIX_SIZE,
    Validator,
)
from cellxgene_schema.write_labels import AnnDataLabelAppender
from dask.array import from_array
from fixtures.examples_validate import visium_library_id

schema_def = get_schema_definition()

# Number of genes in valid adata
NUMBER_OF_GENES = examples.NUMBER_OF_GENES


@pytest.fixture(scope="module")
def validator() -> Validator:
    validator = Validator()

    # Override the schema definition here
    validator._set_schema_def()

    # lower threshold for low gene count warning
    validator.schema_def["components"]["var"]["warn_if_less_than_rows"] = 1
    return validator


@pytest.fixture
def validator_with_adata(validator) -> Validator:
    validator.adata = examples.adata.copy()
    return validator


@pytest.fixture
def validator_with_mouse_adata(validator) -> Validator:
    validator.adata = examples.adata_mouse.copy()
    return validator


@pytest.fixture
def validator_with_tissue_type(validator_with_adata) -> Validator:
    obs = validator_with_adata.adata.obs
    obs.loc[obs.index[0], "tissue_type"] = "tissue"
    obs.loc[obs.index[0], "tissue_ontology_term_id"] = "UBERON:0001062"
    return validator


@pytest.fixture()
def validator_with_validated_adata(validator_with_adata) -> Validator:
    validator_with_adata.validate_adata()
    return validator_with_adata


@pytest.fixture
def adata_with_labels() -> anndata.AnnData:
    # Manually created  data (positive control)
    return examples.adata_with_labels.copy()


@pytest.fixture
def validator_with_adata_missing_raw(validator) -> Validator:
    validator.adata = examples.adata_non_raw.copy()
    return validator


@pytest.fixture
def validator_with_spatial_and_is_single_false(validator) -> Validator:
    validator.adata = examples.adata_spatial_is_single_false.copy()
    return validator


@pytest.fixture
def validator_with_visium_assay(validator) -> Validator:
    validator.adata = examples.adata_visium.copy()
    validator.reset(None, None)

    return validator


@pytest.fixture
def validator_with_slide_seq_v2_assay(validator) -> Validator:
    validator.adata = examples.adata_slide_seqv2.copy()
    return validator


@pytest.fixture
def label_writer() -> AnnDataLabelAppender:
    """
    Fixture that returns an AnnDataLabelAppender object
    """
    label_writer = AnnDataLabelAppender(examples.adata.copy())
    label_writer._add_labels()
    label_writer.adata.uns["schema_version"] = "4.0.0"
    label_writer.adata.uns["schema_reference"] = (
        "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/4.0.0/schema.md"
    )
    label_writer.adata.uns["citation"] = (
        "Publication: <doi> Dataset Version: https://datasets.cellxgene.cziscience.com/<dataset_version_id>.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/<collection_id>"
    )
    return label_writer


class TestValidAnndata:
    """
    Tests a valid AnnData object. Most other tests below modify this AnnData object and test for failure cases.

    The valid AnnData object has all valid cases described in the schema.
    """

    def test_valid_anndata(self, validator_with_adata):
        validator = validator_with_adata
        validator.validate_adata()
        assert not validator.errors


class TestH5adValidation:
    """
    Checks that validation from h5ad works, only does one invalid example as extensive testing is done in the classes
    below
    """

    def test_validate(self, validator):
        h5ad_valid_file = examples.h5ad_valid
        assert validator.validate_adata(h5ad_valid_file)

    def test_invalidate(self, validator):
        h5ad_invalid_file = examples.h5ad_invalid
        assert not validator.validate_adata(h5ad_invalid_file)


class TestExpressionMatrix:
    """
    Fail cases for expression matrices (anndata.X and anndata.raw.X)
    """

    def test_shapes(self, validator_with_adata):
        """
        All matrix layers MUST have the same shape, and have the same cell labels and gene labels.
        """
        validator = validator_with_adata
        # Creates a raw layer
        validator.adata.raw = validator.adata
        validator.adata.raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        validator.adata.X = examples.adata_non_raw.X.copy()

        # remove one gene
        validator.adata = validator.adata[:, 1:]
        validator.validate_adata()
        assert (
            f"ERROR: Number of genes in X ({NUMBER_OF_GENES - 1}) is different than raw.X ({NUMBER_OF_GENES})."
            in validator.errors
        )

    def test_csc_matrix_invalid(self, validator_with_adata):
        is_valid_before = validator_with_adata.validate_adata()
        assert is_valid_before

        # Update X to be a CSC matrix
        X = from_array(scipy.sparse.csc_matrix([[1, 2, 3, 4, 5, 6, 7], [0, 1, 2, 3, 4, 5, 6]], dtype=numpy.float32))
        validator_with_adata.adata.X = X
        is_valid_after = validator_with_adata.validate_adata()
        assert not is_valid_after

        assert validator_with_adata.errors == [
            "ERROR: Invalid sparse encoding for X with encoding csc. Only csr sparse encodings are supported.",
        ]

    def test_sparsity(self, validator_with_adata):
        """
        In any layer, if a matrix has 50% or more values that are zeros, the matrix MUST be encoded as a scipy.sparse.csr_matrix
        """
        validator = validator_with_adata
        sparse_X = numpy.zeros([validator.adata.obs.shape[0], validator.adata.var.shape[0]], dtype=numpy.float32)
        sparse_X[0, 1] = 1
        sparse_X[1, 1] = 1
        validator.adata.X = from_array(sparse_X)
        validator.validate_adata()
        assert (
            " which is greater than 0.5, and it is not a 'scipy.sparse.csr_matrix'. The matrix MUST use this type of matrix for the given sparsity."
            in validator.errors[0]
        )

    @pytest.mark.parametrize("invalid_value", [1.5, -1])
    def test_raw_values__invalid(self, validator_with_adata, invalid_value):
        """
        When both `adata.X` and `adata.raw.X` are present, but `adata.raw.X` contains negative or non-integer values,
        an error is raised.
        """

        validator = validator_with_adata
        validator.adata.raw.X[0, 1] = invalid_value
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: All non-zero values in raw matrix must be positive integers of type numpy.float32.",
            "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements.",
        ]

    @pytest.mark.parametrize("invalid_value", [1.5, -1])
    def test_raw_values__invalid_spatial(self, validator_with_visium_assay, invalid_value):
        """
        When both `adata.X` and `adata.raw.X` are present, but `adata.raw.X` contains negative or non-integer values,
        an error is raised.
        """

        validator = validator_with_visium_assay
        validator.adata.raw.X[0, 1] = invalid_value
        validator.reset(None, 2)
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: All non-zero values in raw matrix must be positive integers of type numpy.float32.",
            "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements.",
        ]

    @pytest.mark.parametrize("datatype", [int, "float64"])
    def test_raw_values__wrong_datatype(self, validator_with_adata, datatype):
        """
        When both `adata.X` and `adata.raw.X` are present, but `adata.raw.X` values are stored as the wrong datatype
        """

        validator = validator_with_adata
        raw = anndata.AnnData(X=validator.adata.raw.X, obs=validator.adata.obs, var=validator.adata.raw.var)
        raw.X = raw.X.astype(datatype)
        validator.adata.raw = raw
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Raw matrix values must have type numpy.float32.",
            "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements.",
        ]

    def test_raw_values__contains_zero_row(self, validator_with_adata):
        """
        When both `adata.X` and `adata.raw.X` are present, but `adata.raw.X` contains a row with all zeros
        """

        validator = validator_with_adata
        validator.adata.raw.X[0] = numpy.zeros(validator.adata.var.shape[0], dtype=numpy.float32)
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Each cell must have at least one non-zero value in its row in the raw matrix.",
            "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements.",
        ]

    def test_raw_values__contains_zero_row_in_tissue_1(self, validator_with_visium_assay):
        """
        Raw Matrix contains a row with all zeros and in_tissue is 1, but no values are in_tissue 0.
        """

        validator: Validator = validator_with_visium_assay
        validator.reset(None, 2)
        validator.adata.obs["in_tissue"] = 1
        validator.adata.X[0] = numpy.zeros(validator.adata.var.shape[0], dtype=numpy.float32)
        validator.adata.raw.X[0] = numpy.zeros(validator.adata.var.shape[0], dtype=numpy.float32)
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Each cell must have at least one non-zero value in its row in the raw matrix.",
            "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements.",
        ]

    def test_raw_values__contains_zero_row_in_tissue_1_mixed_in_tissue_values(self, validator_with_visium_assay):
        """
        Raw Matrix contains a row with all zeros and in_tissue is 1, and there are also values with in_tissue 0.
        """

        validator: Validator = validator_with_visium_assay
        validator.adata.X[1] = numpy.zeros(validator.adata.var.shape[0], dtype=numpy.float32)
        validator.adata.raw.X[1] = numpy.zeros(validator.adata.var.shape[0], dtype=numpy.float32)
        validator.reset(None, 2)
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Each observation with obs['in_tissue'] == 1 must have at least one "
            "non-zero value in its row in the raw matrix.",
            "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements.",
        ]

    def test_raw_values__contains_all_zero_rows_in_tissue_0(self, validator_with_visium_assay):
        """
        When both `adata.X` and `adata.raw.X` are present, but `adata.raw.X` contains all rows with all zeros
        and in_tissue is 0
        """

        validator = validator_with_visium_assay
        validator.adata.obs["in_tissue"] = 0
        validator.adata.obs["cell_type_ontology_term_id"] = "unknown"
        validator.adata.X = from_array(
            numpy.zeros([validator.adata.obs.shape[0], validator.adata.var.shape[0]], dtype=numpy.float32)
        )

        # make sure it's encoded as a sparse matrix to isolate the validation to the contents of the matrix.
        validator.adata.X = validator.adata.X.map_blocks(scipy.sparse.csr_matrix)

        validator.adata.raw = validator.adata.copy()
        validator.adata.raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        validator.reset(None, 2)
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: If obs['in_tissue'] contains at least one value 0, then there must be at least "
            "one row with obs['in_tissue'] == 0 that has a non-zero value in the raw matrix.",
            "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements.",
        ]

    def test_raw_values__contains_some_zero_rows_in_tissue_0(self, validator_with_visium_assay):
        """
        When both `adata.X` and `adata.raw.X` are present, but `adata.raw.X` contains some rows with all zeros
        and in_tissue is 0. Success case.
        """

        validator = validator_with_visium_assay
        validator.adata.obs["in_tissue"] = 0
        validator.adata.obs["cell_type_ontology_term_id"] = "unknown"
        validator.adata.X[0] = numpy.zeros(validator.adata.var.shape[0], dtype=numpy.float32)
        validator.adata.raw.X[0] = numpy.zeros(validator.adata.var.shape[0], dtype=numpy.float32)
        validator.reset(None, 2)
        validator.validate_adata()
        assert validator.errors == []

    @pytest.mark.parametrize(
        "assay_ontology_term_id, req_matrix_size, image_size",
        [
            ("EFO:0022858", VISIUM_AND_IS_SINGLE_TRUE_MATRIX_SIZE, SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE),
            (
                "EFO:0022860",
                VISIUM_11MM_AND_IS_SINGLE_TRUE_MATRIX_SIZE,
                SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE_VISIUM_11MM,
            ),
        ],
    )
    def test_raw_values__invalid_visium_and_is_single_true_row_length(
        self, validator_with_visium_assay, assay_ontology_term_id, req_matrix_size, image_size
    ):
        """
        Dataset is visium and uns['is_single'] is True, but raw.X is the wrong length.
        """
        validator: Validator = validator_with_visium_assay
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id

        # hires image size must be present in order to validate the raw.
        validator.adata.uns["spatial"][visium_library_id]["images"]["hires"] = numpy.zeros(
            (1, image_size, 3), dtype=numpy.uint8
        )

        validator.validate_adata()
        if assay_ontology_term_id == ASSAY_VISIUM_11M:
            _errors = [
                f"ERROR: When {ERROR_SUFFIX_VISIUM_11M} and {ERROR_SUFFIX_IS_SINGLE}, the raw matrix must be the "
                "unfiltered feature-barcode matrix 'raw_feature_bc_matrix'. It must have exactly "
                f"{validator.visium_and_is_single_true_matrix_size} rows. Raw matrix row count is 2.",
                "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements.",
            ]
        else:
            _errors = [
                f"ERROR: When {ERROR_SUFFIX_VISIUM} and {ERROR_SUFFIX_IS_SINGLE}, the raw matrix must be the "
                "unfiltered feature-barcode matrix 'raw_feature_bc_matrix'. It must have exactly "
                f"{validator.visium_and_is_single_true_matrix_size} rows. Raw matrix row count is 2.",
                "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements.",
            ]

        assert validator.errors == _errors

    @pytest.mark.parametrize(
        "assay_ontology_term_id, req_matrix_size, image_size",
        [
            ("EFO:0022858", VISIUM_AND_IS_SINGLE_TRUE_MATRIX_SIZE, SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE),
            (
                "EFO:0022860",
                VISIUM_11MM_AND_IS_SINGLE_TRUE_MATRIX_SIZE,
                SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE_VISIUM_11MM,
            ),
        ],
    )
    def test_raw_values__multiple_invalid_in_tissue_errors(
        self, validator_with_visium_assay, assay_ontology_term_id, req_matrix_size, image_size
    ):
        """
        Dataset is visium and uns['is_single'] is True, in_tissue has both 0 and 1 values and there
        are issues validating rows of both in the matrix.
        """

        validator = validator_with_visium_assay

        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        # hires image size must be present in order to validate the raw.
        validator._visium_and_is_single_true_matrix_size = None
        validator._hires_max_dimension_size = image_size
        validator.adata.uns["spatial"][visium_library_id]["images"]["hires"] = numpy.zeros(
            (1, image_size, 3), dtype=numpy.uint8
        )
        validator.adata.X = from_array(
            numpy.zeros([validator.adata.obs.shape[0], validator.adata.var.shape[0]], dtype=numpy.float32)
        )

        # encode X with csr
        validator.adata.X = validator.adata.X.map_blocks(scipy.sparse.csr_matrix)

        validator.adata.raw = validator.adata.copy()
        validator.adata.raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        validator.validate_adata()
        if assay_ontology_term_id == ASSAY_VISIUM_11M:
            assert (
                validator.errors[0]
                == f"ERROR: When {ERROR_SUFFIX_VISIUM_11M} and {ERROR_SUFFIX_IS_SINGLE}, the raw matrix must be the "
                "unfiltered feature-barcode matrix 'raw_feature_bc_matrix'. It must have exactly "
                f"{validator.visium_and_is_single_true_matrix_size} rows. Raw matrix row count is 2."
            )
        else:
            assert (
                validator.errors[0]
                == f"ERROR: When {ERROR_SUFFIX_VISIUM} and {ERROR_SUFFIX_IS_SINGLE}, the raw matrix must be the "
                "unfiltered feature-barcode matrix 'raw_feature_bc_matrix'. It must have exactly "
                f"{validator.visium_and_is_single_true_matrix_size} rows. Raw matrix row count is 2."
            )

        assert validator.errors[1:] == [
            "ERROR: If obs['in_tissue'] contains at least one value 0, then there must be at least "
            "one row with obs['in_tissue'] == 0 that has a non-zero value in the raw matrix.",
            "ERROR: Each observation with obs['in_tissue'] == 1 must have at least one "
            "non-zero value in its row in the raw matrix.",
            "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements.",
        ]

    def test_raw_values__multiple_errors(self, validator_with_adata):
        """
        When both `adata.X` and `adata.raw.X` are present, but `adata.raw.X` contains multiple errors and all are
        reported
        """

        validator = validator_with_adata
        validator.adata.raw.X[0] = numpy.zeros(validator.adata.var.shape[0], dtype=numpy.float32)
        validator.adata.raw.X[1, 1] = 1.5
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Each cell must have at least one non-zero value in its row in the raw matrix.",
            "ERROR: All non-zero values in raw matrix must be positive integers of type numpy.float32.",
            "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements.",
        ]

    def test_raw_values__non_rna(self, validator_with_adata):
        """
        Except for ATAC-seq and methylation data, raw data is REQUIRED
        """

        # ATAC - raw layer not required
        # The assignment above makes X to not be raw: validator.adata.uns["X_normalization"] = "CPM"
        # The following line makes it to be scATAC-seq data (EFO:0010891)
        # Missing raw data in atac-seq data is allowed, thus the following should not return an error message
        validator = validator_with_adata
        obs = validator.adata.obs
        validator.errors = []
        obs["assay_ontology_term_id"] = "EFO:0010891"
        obs["suspension_type"] = "nucleus"
        obs["suspension_type"] = obs["suspension_type"].astype("category")
        validator.validate_adata()
        assert validator.errors == []

    def test_raw_values__matrix_chunks(self, validator_with_adata):
        """
        Test adata is validated correctly when matrix is larger than the chunk size
        """
        with unittest.mock.patch.object(read_h5ad, "__defaults__", (1,)):
            validator = validator_with_adata
            validator.validate_adata()
            assert validator.errors == []

    def test_missing_raw_matrix(self, validator_with_adata_missing_raw):
        """
        Test error message appears if dataset with RNA assay is missing raw matrix
        """
        validator = validator_with_adata_missing_raw
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: All non-zero values in raw matrix must be positive integers of type numpy.float32.",
            "ERROR: Raw data is missing: there is only a normalized matrix in X and no raw.X",
        ]

    def test_final_strongly_recommended(self, validator_with_adata):
        """
        Except for ATAC-seq and methylation data, final matrix is STRONGLY RECOMMENDED
        """

        # move raw to X amd: i.e. there is no final
        validator = validator_with_adata
        validator.adata.X = validator.adata.raw.X
        del validator.adata.raw
        validator.validate_adata()
        assert (
            "WARNING: Only raw data was found, i.e. there is no 'raw.X'. It is STRONGLY RECOMMENDED that 'final' (normalized) data is provided."
            in validator.warnings
        )


class TestObs:
    """
    Fail cases in adata.obs
    """

    def get_format_error_message(self, error_message_suffix: str, detail: str) -> str:
        return detail + " " + error_message_suffix

    @pytest.mark.parametrize(
        "column",
        [
            "development_stage_ontology_term_id",
            "disease_ontology_term_id",
            "self_reported_ethnicity_ontology_term_id",
            "is_primary_data",
            "sex_ontology_term_id",
            "tissue_ontology_term_id",
            "donor_id",
            "suspension_type",
        ],
    )
    def test_column_presence(self, validator_with_adata, column):
        """
        obs is a pandas.DataFrame. Curators MUST annotate the following columns in the obs dataframe.
        """
        validator = validator_with_adata
        validator.adata.obs.drop(column, axis=1, inplace=True)
        # Remove batch condition because it has a dependency with is_primary_data
        validator.adata.uns.pop("batch_condition")

        validator.validate_adata()
        assert f"ERROR: Dataframe 'obs' is missing " f"column '{column}'." in validator.errors

    def test_key_presence_organism_ontology_term_id(self, validator_with_adata):
        """
        organism_ontology_term_id must be defined in uns
        """
        validator = validator_with_adata
        del validator.adata.uns["organism_ontology_term_id"]
        validator.validate_adata()
        assert len(validator.errors) > 0
        assert validator.errors[0] == "ERROR: 'organism_ontology_term_id' in 'uns' is not present."

    def test_key_presence_organism(self, validator_with_adata):
        """
        organism_ontology_term_id must be defined in uns
        """
        validator = validator_with_adata
        validator.adata.uns["organism"] = "Homo sapiens"
        validator.validate_adata()
        assert len(validator.errors) > 0
        assert (
            validator.errors[0]
            == "ERROR: Add labels error: Column 'organism' is a reserved column name of 'uns'. Remove it from h5ad and try again."
        )

    @pytest.mark.parametrize(
        "deprecated_column",
        [
            "organism_ontology_term_id",
            "organism",
        ],
    )
    def test_obs_presence_organism(self, validator_with_adata, deprecated_column):
        """
        organism_ontology_term_id cannot be in obs since it's a deprecated column
        """
        validator = validator_with_adata
        validator.adata.obs[deprecated_column] = "NCBITaxon:9606"
        validator.validate_adata()
        assert f"ERROR: The field '{deprecated_column}' is present in 'obs', but it is deprecated." in validator.errors

    @pytest.mark.parametrize(
        "organism_ontology_term_id",
        [
            "NCBITaxon:2697049",  # Severe acute respiratory syndrome coronavirus 2
            "NCBITaxon:32630",  # synthetic construct
        ],
    )
    def test_exempt_organism_values(self, validator_with_adata, organism_ontology_term_id):
        validator = validator_with_adata
        validator.adata.uns["organism_ontology_term_id"] = organism_ontology_term_id
        validator.validate_adata()
        assert (
            f"ERROR: '{organism_ontology_term_id}' in 'organism_ontology_term_id' is not allowed." in validator.errors
        )

    def test_invalid_organism_genes(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["organism_ontology_term_id"] = "NCBITaxon:10090"
        validator.validate_adata()
        assert (
            "ERROR: uns['organism_ontology_term_id'] is 'NCBITaxon:10090' but feature_ids are from ['NCBITaxon:9606']."
        )

    def test_column_presence_assay(self, validator_with_adata):
        """
        obs is a pandas.DataFrame. Curators MUST annotate the following columns in the obs dataframe.

        A separate check is need for assay_ontology_term_id because removing from anndata results in multiple
        errors given that other columns depend on its presence
        """
        validator = validator_with_adata
        validator.adata.obs.drop("assay_ontology_term_id", axis=1, inplace=True)
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Dataframe 'obs' is missing column " "'assay_ontology_term_id'.",
        ]

    @pytest.mark.parametrize(
        "assay_ontology_term_id, is_descendant",
        [("EFO:0022859", True), ("EFO:0022858", True), ("EFO:0030029", False), ("EFO:0002697", False)],
    )
    def test_column_presence_in_tissue(self, validator_with_visium_assay, assay_ontology_term_id, is_descendant):
        """
        Spatial assays that are descendants of visium must have a valid "in_tissue" column.
        """
        validator: Validator = validator_with_visium_assay

        # reset and test
        validator.reset()
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        validator.adata.obs["assay_ontology_term_id"] = validator.adata.obs["assay_ontology_term_id"].astype("category")
        validator._validate_spatial_tissue_position("in_tissue", 0, 1)
        if is_descendant:
            assert validator.errors == []
        else:
            assert validator.errors == [
                f"obs['in_tissue'] is only allowed for {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE}."
            ]

    @pytest.mark.parametrize("reserved_column", schema_def["components"]["obs"]["reserved_columns"])
    def test_obs_reserved_columns_presence(self, validator_with_adata, reserved_column):
        """
        Reserved columns must NOT be used in obs
        """
        validator = validator_with_adata
        validator.adata.obs[reserved_column] = "dummy_value"
        validator.validate_adata()
        assert validator.errors == [
            f"ERROR: Column '{reserved_column}' is a reserved column name "
            f"of 'obs'. Remove it from h5ad and try again."
        ]

    def test_obsolete_term_id(self, validator_with_adata):
        """
        Terms documented as obsolete in an ontology MUST NOT be used. For example, EFO:0009310
        for obsolete_10x v2 was marked as obsolete in EFO version 3.31.0 and replaced by
        EFO:0009899 for 10x 3' v2.

        https://www.ebi.ac.uk/ols/ontologies/efo/terms?short_form=EFO_0009310
        """

        # Not a valid term
        validator = validator_with_adata
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"]["assay_ontology_term_id"][
            "error_message_suffix"
        ]
        validator.adata.obs.loc[validator.adata.obs.index[0], "assay_ontology_term_id"] = "EFO:0009310"
        validator.validate_adata()
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'EFO:0009310' in 'assay_ontology_term_id' is a deprecated term id of 'EFO'.",
            )
        ]

    @pytest.mark.parametrize(
        "assay_ontology_term_id,error",
        [
            ("CL:000001", "ERROR: 'CL:000001' in 'assay_ontology_term_id' is not a valid ontology term id of 'EFO'."),
            ("EFO:0000001", "ERROR: 'EFO:0000001' in 'assay_ontology_term_id' is not an allowed term id."),
            ("EFO:0002772", "ERROR: 'EFO:0002772' in 'assay_ontology_term_id' is not an allowed term id."),
            ("EFO:0010183", "ERROR: 'EFO:0010183' in 'assay_ontology_term_id' is not an allowed term id."),
            (
                "EFO:0010183 (sci-plex)",
                "ERROR: 'EFO:0010183 (sci-plex)' in 'assay_ontology_term_id' is not a valid ontology term id of 'EFO'.",
            ),
        ],
    )
    def test_assay_ontology_term_id(self, validator_with_adata, assay_ontology_term_id, error):
        """
        assay_ontology_term_id categorical with str categories.
        This MUST be an EFO term that is a descendant of either "EFO:0002772" or "EFO:0010183"
        """
        validator = validator_with_adata
        validator.adata.obs.loc[validator.adata.obs.index[0], "assay_ontology_term_id"] = assay_ontology_term_id
        validator.validate_adata()
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"]["assay_ontology_term_id"][
            "error_message_suffix"
        ]
        assert validator.errors == [self.get_format_error_message(error_message_suffix, error)]

    def test_assay_ontology_term_id__as_categorical(self, validator_with_visium_assay):
        """
        Formally, assay_ontology_term_id is expected to be a categorical variable of type string. However, it should work for categorical dtypes as well.
        """
        validator: Validator = validator_with_visium_assay

        # check encoding as string
        validator.reset(None, 2)
        validator._check_spatial()
        validator._validate_raw()
        assert validator.errors == []

        # force encoding as 'categorical'
        validator.reset(None, 2)
        validator.adata.obs["assay_ontology_term_id"] = validator.adata.obs["assay_ontology_term_id"].astype("category")
        validator._check_spatial()
        validator._validate_raw()
        assert validator.errors == []

    @pytest.mark.parametrize(
        "assay_ontology_term_id, all_same",
        [("EFO:0022859", True), ("EFO:0030062", True), ("EFO:0022860", True), ("EFO:0008995", False)],
    )
    def test_assay_ontology_term_id__all_same(self, validator_with_visium_assay, assay_ontology_term_id, all_same):
        """
        Spatial assays (descendants of Visium Spatia Gene Expression, or Slide-SeqV2) require all values in the column to be identical.
        """
        validator: Validator = validator_with_visium_assay

        # mix values (with otherwise allowed values)
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        validator.adata.obs["assay_ontology_term_id"].iloc[0] = "EFO:0010183"

        # check that unique values are allowed
        validator._check_spatial_obs()
        EXPECTED_ERROR = f"When {ERROR_SUFFIX_SPATIAL}, all observations must contain the same value."
        if all_same:
            assert EXPECTED_ERROR in validator.errors
        else:
            assert validator.errors not in validator.errors

    def test_cell_type_ontology_term_id_invalid_term(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.obs.loc[validator.adata.obs.index[0], "cell_type_ontology_term_id"] = "EFO:0000001"
        validator.validate_adata()
        assert len(validator.errors) > 0

    @pytest.mark.parametrize("term", ["CL:0000000"])
    def test_cell_type_ontology_term_id(self, validator_with_adata, term):
        validator = validator_with_adata
        validator.adata.obs["cell_type_ontology_term_id"] = term
        validator.validate_adata()
        assert len(validator.errors) == 0

    @pytest.mark.parametrize(
        "term",
        schema_def["components"]["obs"]["columns"]["cell_type_ontology_term_id"]["curie_constraints"]["forbidden"][
            "terms"
        ],
    )
    def test_cell_type_ontology_term_id__forbidden(self, validator_with_adata, term):
        """
        cell_type_ontology_term_id categorical with str categories. This MUST be a CL term, and must NOT match forbidden
        columns defined in schema
        """
        validator = validator_with_adata
        validator.adata.obs.loc[validator.adata.obs.index[0], "cell_type_ontology_term_id"] = term
        validator.validate_adata()
        # Forbidden columns may be marked as either "not allowed" or "deprecated"
        not_allowed_error_message = f"ERROR: '{term}' in 'cell_type_ontology_term_id' is not allowed."
        deprecated_error_message = f"ERROR: '{term}' in 'cell_type_ontology_term_id' is a deprecated term id of 'CL'."
        assert not_allowed_error_message in validator.errors or deprecated_error_message in validator.errors

    @pytest.mark.parametrize(
        "organism_ontology_term_id",
        [
            "NCBITaxon:10090",  # Mus musculus (ancestor term, always allowed)
            "NCBITaxon:947985",  # Mus musculus albula
            "NCBITaxon:35531",  # Mus musculus bactrianus
        ],
    )
    def test_development_stage_ontology_term_id_mouse_descendant(
        self, validator_with_mouse_adata, organism_ontology_term_id
    ):
        """
        If organism_ontology_term_id is "NCBITaxon:10090" for Mus musculus OR its descendants,
        this MUST be the most accurate MmusDv:0000001 descendant.
        """
        validator = validator_with_mouse_adata
        validator.adata.uns["organism_ontology_term_id"] = organism_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    @pytest.mark.parametrize(
        "development_stage_ontology_term_id,error",
        [
            (
                "CL:000001",
                "ERROR: 'CL:000001' in 'development_stage_ontology_term_id' is not a valid ontology term id of 'HsapDv'.",
            ),
            (
                "HsapDv:0000000",
                "ERROR: 'HsapDv:0000000' in 'development_stage_ontology_term_id' is not an allowed term id.",
            ),
            (
                "HsapDv:0000001",
                "ERROR: 'HsapDv:0000001' in 'development_stage_ontology_term_id' is not an allowed term id.",
            ),
        ],
    )
    def test_development_stage_ontology_term_id_human(
        self, validator_with_adata, development_stage_ontology_term_id, error
    ):
        """
        development_stage_ontology_term_id categorical with str categories. If unavailable, this MUST be "unknown".
        If organism_ontology_term_id is "NCBITaxon:9606" for Homo sapiens,
        this MUST be the most accurate HsapDv:0000001 descendant.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = development_stage_ontology_term_id
        validator.validate_adata()
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "development_stage_ontology_term_id"
        ]["dependencies"][0]["error_message_suffix"]
        assert validator.errors == [self.get_format_error_message(error_message_suffix, error)]

    @pytest.mark.parametrize(
        "development_stage_ontology_term_id,error",
        [
            (
                "CL:000001",
                "ERROR: 'CL:000001' in 'development_stage_ontology_term_id' is not a valid ontology term id of 'MmusDv'.",
            ),
            (
                "MmusDv:0000000",
                "ERROR: 'MmusDv:0000000' in 'development_stage_ontology_term_id' is not an allowed term id.",
            ),
            (
                "MmusDv:0000001",
                "ERROR: 'MmusDv:0000001' in 'development_stage_ontology_term_id' is not an allowed term id.",
            ),
        ],
    )
    def test_development_stage_ontology_term_id_mouse(
        self, validator_with_mouse_adata, development_stage_ontology_term_id, error
    ):
        """
        If organism_ontology_term_id is "NCBITaxon:10090" for Mus musculus,
        this MUST be the most accurate MmusDv:0000001 descendant.
        """
        validator = validator_with_mouse_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = development_stage_ontology_term_id
        validator.validate_adata()
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "development_stage_ontology_term_id"
        ]["dependencies"][1]["error_message_suffix"]
        assert validator.errors == [self.get_format_error_message(error_message_suffix, error)]

    def test_development_stage_ontology_term_id_all_species(self, validator_with_adata):
        """
        All other it MUST be descendants of UBERON:0000105 and not UBERON:0000071
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        # Fail case not an UBERON term
        validator.adata.uns["organism_ontology_term_id"] = "NCBITaxon:9541"
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = "EFO:0000001"
        obs.loc[
            obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "na"
        validator.validate_adata()
        assert (
            "ERROR: 'EFO:0000001' in 'development_stage_ontology_term_id' is not a valid ontology term id of 'UBERON'. When 'organism_ontology_term_id'-specific requirements are not defined in the schema definition, 'development_stage_ontology_term_id' MUST be a descendant term id of 'UBERON:0000105' excluding 'UBERON:0000071', or unknown."
            in validator.errors
        )

        # All other it MUST be descendants of UBERON:0000105 and not UBERON:0000071
        # Fail case UBERON:0000071
        validator.errors = []
        validator.adata.uns["organism_ontology_term_id"] = "NCBITaxon:9541"
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = "UBERON:0000071"
        obs.loc[
            obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "na"
        validator.validate_adata()
        assert (
            "ERROR: 'UBERON:0000071' in 'development_stage_ontology_term_id' is not allowed. When 'organism_ontology_term_id'-specific requirements are not defined in the schema definition, 'development_stage_ontology_term_id' MUST be a descendant term id of 'UBERON:0000105' excluding 'UBERON:0000071', or unknown."
            in validator.errors
        )

    def test_disease_ontology_term_id(self, validator_with_adata):
        """
        disease_ontology_term_id categorical with str categories. This MUST be one of:
        - PATO:0000461 for normal or healthy
        - one or more valid MONDO terms

        MONDO terms must be one of:
        - descendant of MONDO:0000001 for disease
        - self or descendant of MONDO:0021178 for injury
        """
        validator = validator_with_adata
        obs = validator.adata.obs

        # Invalid ontology
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "EFO:0000001"
        validator.validate_adata()
        assert (
            "ERROR: 'EFO:0000001' in 'disease_ontology_term_id' is not a valid ontology term id of 'MONDO, PATO'."
            in validator.errors[0]
        )

        # Invalid PATO term id
        validator.errors = []
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "PATO:0001894"
        validator.validate_adata()
        assert "ERROR: 'PATO:0001894' in 'disease_ontology_term_id' is not an allowed term id." in validator.errors[0]

        # Invalid MONDO term id - disease characteristic
        validator.errors = []
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "MONDO:0021125"
        validator.validate_adata()
        assert "ERROR: 'MONDO:0021125' in 'disease_ontology_term_id' is not an allowed term id." in validator.errors[0]

        # Invalid MONDO term id - disease parent term
        validator.errors = []
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "MONDO:0000001"
        validator.validate_adata()
        assert "ERROR: 'MONDO:0000001' in 'disease_ontology_term_id' is not an allowed term id." in validator.errors[0]

        # Valid PATO term id - healthy
        validator.errors = []
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "PATO:0000461"
        validator.validate_adata()
        assert validator.errors == []

        # Valid MONDO term id - disease descendant term
        validator.errors = []
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "MONDO:0005491"
        validator.validate_adata()
        assert validator.errors == []

        # Valid MONDO term id - injury parent term
        validator.errors = []
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "MONDO:0021178"
        validator.validate_adata()
        assert validator.errors == []

        # Valid MONDO term id - injury descendant term
        validator.errors = []
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "MONDO:0015796"
        validator.validate_adata()
        assert validator.errors == []

        # Valid MONDO term ids - disease descendant term, injury descendant term
        validator.errors = []
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "MONDO:0005491 || MONDO:0015796"
        validator.validate_adata()
        assert validator.errors == []

        # Valid PATO term id and valid MONDO term id. However, multiple terms are only allowed for MONDO terms
        validator.errors = []
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "MONDO:0005491 || PATO:0000461"
        validator.validate_adata()
        assert (
            "ERROR: 'PATO:0000461' in 'disease_ontology_term_id' is not a valid ontology term id of 'MONDO'."
            in validator.errors[0]
        )

    def test_self_reported_ethnicity_ontology_term_id__unknown(self, validator_with_adata):
        """
        Test 'unknown' self_reported_ethnicity_ontology_term is valid when used as single term
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[
            obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "unknown"
        assert validator.validate_adata()
        assert validator.errors == []

    def test_cell_type_ontology_term_id__unknown(self, validator_with_adata):
        """
        Test 'unknown' cell_type_ontology_term_id is valid
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.at["Y", "cell_type_ontology_term_id"] = "unknown"
        assert validator.validate_adata()
        assert validator.errors == []

    def test_tissue_ontology_term_id__unknown(self, validator_with_adata):
        """
        Test 'unknown' tissue_ontology_term_id is valid if tissue_type is 'cell culture'
        """
        validator = validator_with_adata
        obs = validator.adata.obs

        # Arrange -- relies on "tissue_type" value for index "Y" being "cell culture", set explicitly
        obs.at["Y", "tissue_type"] = "cell culture"
        obs.at["Y", "tissue_ontology_term_id"] = "unknown"

        assert validator.validate_adata()
        assert validator.errors == []

    def test_tissue_ontology_term_id__unknown_invalid(self, validator_with_adata):
        """
        Test 'unknown' tissue_ontology_term_id is invalid if tissue_type is NOT 'cell culture'
        """
        validator = validator_with_adata
        obs = validator.adata.obs

        # Arrange -- 'tissue_ontology_term_id' cannot be "unknown" when 'tissue_type is "tissue"
        obs.at["Y", "tissue_type"] = "tissue"
        obs.at["Y", "tissue_ontology_term_id"] = "unknown"

        assert not validator.validate_adata()
        assert len(validator.errors) > 0

    def test_self_reported_ethnicity_ontology_term_id__unknown_in_multi_term(self, validator_with_adata):
        """
        Test 'unknown' self_reported_ethnicity_ontology_term is invalid when used in multi-term comma-delimited str
        """
        validator = validator_with_adata
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "self_reported_ethnicity_ontology_term_id"
        ]["dependencies"][0]["error_message_suffix"]

        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "HANCESTRO:0005 || HANCESTRO:0014 || unknown"
        validator.validate_adata()
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'unknown' in 'self_reported_ethnicity_ontology_term_id' is not a valid ontology term id "
                "of 'HANCESTRO'.",
            )
        ]

    def test_self_reported_ethnicity_ontology_term_id__invalid_ontology(self, validator_with_adata):
        """
        Test self_reported_ethnicity_ontology_term error message when passed a valid term from an invalid ontology
        """
        validator = validator_with_adata
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "self_reported_ethnicity_ontology_term_id"
        ]["dependencies"][0]["error_message_suffix"]
        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "EFO:0000001"
        validator.validate_adata()
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'EFO:0000001' in 'self_reported_ethnicity_ontology_term_id' is not a valid ontology term "
                "id of 'HANCESTRO'.",
            )
        ]

    def test_self_reported_ethnicity_ontology_term_id__forbidden_term(self, validator_with_adata):
        """
        Test self_reported_ethnicity_ontology_term error message when passed an explicitly forbidden ontology term that
        is otherwise valid
        """
        validator = validator_with_adata
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "self_reported_ethnicity_ontology_term_id"
        ]["dependencies"][0]["error_message_suffix"]

        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "HANCESTRO:0003"
        validator.validate_adata()
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'HANCESTRO:0003' in 'self_reported_ethnicity_ontology_term_id' is not allowed.",
            )
        ]

    def test_self_reported_ethnicity_ontology_term_id__forbidden_term_ancestor(self, validator_with_adata):
        """
        Test self_reported_ethnicity_ontology_term error message when passed an ontology term that has
        both itself and its descendants forbidden
        """
        validator = validator_with_adata
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "self_reported_ethnicity_ontology_term_id"
        ]["dependencies"][0]["error_message_suffix"]

        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "HANCESTRO:0002"
        validator.validate_adata()
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'HANCESTRO:0002' in 'self_reported_ethnicity_ontology_term_id' is not allowed.",
            )
        ]

    def test_self_reported_ethnicity_ontology_term_id__forbidden_term_descendant(self, validator_with_adata):
        """
        Test self_reported_ethnicity_ontology_term error message when passed the descendant term of an ontology term that has
        both itself and its descendants forbidden
        """
        validator = validator_with_adata
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "self_reported_ethnicity_ontology_term_id"
        ]["dependencies"][0]["error_message_suffix"]

        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "HANCESTRO:0306"
        validator.validate_adata()
        assert (
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'HANCESTRO:0306' in 'self_reported_ethnicity_ontology_term_id' is not allowed. Descendant terms "
                "of 'HANCESTRO:0304' are not allowed.",
            )
            in validator.errors
        )

    def test_self_reported_ethnicity_ontology_term_id__forbidden_term_in_multi_term(self, validator_with_adata):
        """
        Test error message for self_reported_ethnicity_ontology_term_id involving a forbidden term among an otherwise
        valid comma-delimited str of multiple terms
        """
        validator = validator_with_adata
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "self_reported_ethnicity_ontology_term_id"
        ]["dependencies"][0]["error_message_suffix"]

        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "HANCESTRO:0005 || HANCESTRO:0014 || HANCESTRO:0018"
        validator.validate_adata()
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'HANCESTRO:0018' in 'self_reported_ethnicity_ontology_term_id' is not allowed.",
            )
        ]

    def test_self_reported_ethnicity_ontology_term_id__non_human(self, validator_with_adata):
        """
        Test self_reported_ethnicity_ontology_term error message if term is not 'na' for a non-human organism
        """
        validator = validator_with_adata
        # Mouse organism ID
        validator.adata.uns["organism_ontology_term_id"] = "NCBITaxon:10090"
        # Required to set to avoid development_stage_ontology_term_id errors
        validator.adata.obs.loc[validator.adata.obs.index[0], "development_stage_ontology_term_id"] = "MmusDv:0000003"
        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "HANCESTRO:0005"
        validator.validate_adata()
        assert (
            "ERROR: 'HANCESTRO:0005' in 'self_reported_ethnicity_ontology_term_id' is not a valid value of 'self_reported_ethnicity_ontology_term_id'. When 'organism_ontology_term_id' is NOT 'NCBITaxon:9606' (Homo sapiens), self_reported_ethnicity_ontology_term_id MUST be 'na'."
            in validator.errors
        )

    def test_self_reported_ethnicity_ontology_term_id__unsorted(self, validator_with_adata):
        """
        Test error message for self_reported_ethnicity_ontology_term_id with valid comma-delimited terms in a str,
         but NOT in ascending lexical order
        """
        validator = validator_with_adata
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "self_reported_ethnicity_ontology_term_id"
        ]["dependencies"][0]["error_message_suffix"]

        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "HANCESTRO:0014 || HANCESTRO:0005"
        validator.validate_adata()
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'HANCESTRO:0014 || HANCESTRO:0005' in 'self_reported_ethnicity_ontology_term_id' is not in "
                "ascending lexical order.",
            )
        ]

    def test_self_reported_ethnicity_ontology_term_id__invalid_delimiters(self, validator_with_adata):
        """
        Test error message for self_reported_ethnicity_ontology_term_id involving
        delimiters that are not specified in the schema definition yaml, such as whitespace
        """
        validator = validator_with_adata

        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "HANCESTRO:0005,HANCESTRO:0014"
        validator.validate_adata()
        assert (
            "'HANCESTRO:0005,HANCESTRO:0014' in 'self_reported_ethnicity_ontology_term_id' is not a valid ontology term id of 'HANCESTRO'"
            in validator.errors[0]
        )

    def test_self_reported_ethnicity_ontology_term_id__multiple_errors_in_multi_term(self, validator_with_adata):
        """
        Test that multiple distinct error messages are reported for self_reported_ethnicity_ontology_term_id with
        multiple different error types in a comma-delimited multi term str
        """
        validator = validator_with_adata
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "self_reported_ethnicity_ontology_term_id"
        ]["dependencies"][0]["error_message_suffix"]

        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "EFO:0000001 || HANCESTRO:0004 || HANCESTRO:0014 || HANCESTRO:1"
        validator.validate_adata()
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'EFO:0000001' in 'self_reported_ethnicity_ontology_term_id' is not a valid ontology term "
                "id of 'HANCESTRO'.",
            ),
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'HANCESTRO:0004' in 'self_reported_ethnicity_ontology_term_id' is not allowed.",
            ),
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'HANCESTRO:1' in 'self_reported_ethnicity_ontology_term_id' is not a valid ontology term "
                "id of 'HANCESTRO'.",
            ),
        ]

    def test_self_reported_ethnicity_ontology_term_id__duplicates_in_multi_term(self, validator_with_adata):
        """
        Test error message for self_reported_ethnicity_ontology_term_id with duplicate valid ontology terms in a
        comma-delimited multi term str
        """
        validator = validator_with_adata
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "self_reported_ethnicity_ontology_term_id"
        ]["dependencies"][0]["error_message_suffix"]

        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "HANCESTRO:0014 || HANCESTRO:0014"
        validator.validate_adata()
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'HANCESTRO:0014 || HANCESTRO:0014' in 'self_reported_ethnicity_ontology_term_id' contains "
                "duplicates.",
            )
        ]

    def test_self_reported_ethnicity_ontology_term_id__multi_term_list(self, validator_with_adata):
        """
        Test error message for self_reported_ethnicity_ontology_term_id with list type instead of str
        """
        validator = validator_with_adata
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "self_reported_ethnicity_ontology_term_id"
        ]["dependencies"][0]["error_message_suffix"]

        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = ["HANCESTRO:0005 || HANCESTRO:0014"]
        validator.validate_adata()
        assert validator.errors[1] == self.get_format_error_message(
            error_message_suffix,
            "ERROR: '['HANCESTRO:0005 || HANCESTRO:0014']' in 'self_reported_ethnicity_ontology_term_id' is not "
            "a valid ontology term value, it must be a string.",
        )

    @pytest.mark.parametrize(
        "organism_ontology_term_id",
        ["EFO:0000001", "CL:4023077", "NCBITaxon:1234", "NCBITaxon:9615"],
    )
    def test_organism_ontology_term_id__invalid(self, validator_with_adata, organism_ontology_term_id):
        """
        organism_ontology_term_id categorical with str categories. This MUST be one of approved enumerated species.
        """
        validator = validator_with_adata
        validator.adata.uns["organism_ontology_term_id"] = organism_ontology_term_id
        validator.validate_adata()
        assert len(validator.errors) > 0

    def test_tissue_ontology_term_id_base(self, validator_with_adata):
        """
        tissue_ontology_term_id categorical with str categories. This MUST be the term that best describes the tissue
        that this cell was derived from, depending on the type of biological sample:
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = "EFO:0000001"
        obs.loc[obs.index[0], "tissue_type"] = "tissue"
        validator.validate_adata()
        assert len(validator.errors) > 0

    def test_tissue_ontology_term_id_cell_culture__suffix_in_term_id(self, validator_with_adata):
        """
        Cell Culture - Cannot include suffixes.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = "cell culture"
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = "CL:0000057 (cell culture)"
        validator.validate_adata()
        assert len(validator.errors) > 0

    def test_tissue_ontology_term_id_cell_culture__not_a_CL_term(self, validator_with_adata):
        """
        Cell Culture - MUST be CL term
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = "cell culture"
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = "EFO:0000001"
        validator.validate_adata()
        assert len(validator.errors) > 0

    @pytest.mark.parametrize(
        "term",
        schema_def["components"]["obs"]["columns"]["cell_type_ontology_term_id"]["curie_constraints"]["forbidden"][
            "terms"
        ],
    )
    def test_tissue_ontology_term_id_cell_culture_3(self, validator_with_adata, term):
        """
        Cell Culture - must be a valid CL term other than forbidden columns in schema definition.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = "cell culture"
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = term
        validator.validate_adata()

        # Forbidden columns may be marked as either "not allowed" or "deprecated"
        not_allowed_error = f"ERROR: '{term}' in 'tissue_ontology_term_id' is not allowed. When 'tissue_type' is 'cell culture', 'tissue_ontology_term_id' MUST follow the validation rules for cell_type_ontology_term_id."
        deprecated_error = f"ERROR: '{term}' in 'tissue_ontology_term_id' is a deprecated term id of 'CL'. When 'tissue_type' is 'cell culture', 'tissue_ontology_term_id' MUST follow the validation rules for cell_type_ontology_term_id."
        invalid_ontology_error = f"ERROR: '{term}' in 'tissue_ontology_term_id' is not a valid ontology term id of 'CL'. When 'tissue_type' is 'cell culture', 'tissue_ontology_term_id' MUST follow the validation rules for cell_type_ontology_term_id."
        assert (
            not_allowed_error in validator.errors
            or deprecated_error in validator.errors
            or invalid_ontology_error in validator.errors
        )

    def test_tissue_ontology_term_id_organoid(self, validator_with_adata):
        """
        Organoid - must not accept suffixes like "(organoid)"
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = "UBERON:0000057 (organoid)"
        obs.loc[obs.index[0], "tissue_type"] = "organoid"
        validator.validate_adata()
        assert (
            "ERROR: 'UBERON:0000057 (organoid)' in 'tissue_ontology_term_id' is not a valid ontology term id"
            in validator.errors[0]
        )

    def test_tissue_ontology_term_id_descendant_of_anatomical_entity__tissue(self, validator_with_adata):
        """
        Tissue ontology term ID must be a descendant term of 'UBERON:0001062' (anatomical entity) if tissue_type is
        organoid or tissue.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = "UBERON:0001062"
        obs.loc[obs.index[0], "tissue_type"] = "tissue"
        validator.validate_adata()
        assert "ERROR: 'UBERON:0001062' in 'tissue_ontology_term_id' is not an allowed term id." in validator.errors[0]

    def test_tissue_ontology_term_id_descendant_of_anatomical_entity__organoid(self, validator_with_adata):
        """
        Tissue ontology term ID must be a descendant term of 'UBERON:0001062' (anatomical entity) if tissue_type is
        organoid or tissue.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = "UBERON:0001062"
        obs.loc[obs.index[0], "tissue_type"] = "organoid"
        validator.validate_adata()
        assert "ERROR: 'UBERON:0001062' in 'tissue_ontology_term_id' is not an allowed term id." in validator.errors[0]

    def test_tissue_type(self, validator_with_adata):
        """
        tissue_type must be one of 'cell culture', 'tissue', or 'organoid'
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.tissue_type = obs.tissue_type.cat.add_categories(["organ"])
        obs.loc[obs.index[0], "tissue_type"] = "organ"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Column 'tissue_type' in dataframe 'obs' contains invalid values "
            "'['organ']'. Values must be one of ['cell culture', 'organoid', 'tissue']"
        ]

    def test_sex_ontology_term_id(self, validator_with_adata):
        """
        sex_ontology_term_id categorical with str categories.
        This MUST be a descendant of PATOPATO:0001894 for phenotypic sex or "unknown" if unavailable
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "sex_ontology_term_id"] = "EFO:0000001"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'EFO:0000001' in 'sex_ontology_term_id' is "
            "not a valid ontology term id of 'PATO'. Only 'PATO:0000383', 'PATO:0000384', 'PATO:0001340', "
            "or 'unknown' are allowed."
        ]

    def test_is_primary_data(self, validator_with_adata):
        """
        is_primary_data	bool. This MUST be True if this is the canonical instance of this cellular
        observation and False if not. This is commonly False
        for meta-analyses reusing data or for secondary views of data.
        """
        validator = validator_with_adata
        validator.adata.obs["is_primary_data"] = "FALSE"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Column 'is_primary_data' in dataframe 'obs' " "must be boolean, not 'object'."
        ]

    def test_is_primary_data__spatial(self, validator_with_spatial_and_is_single_false):
        """
        is_primary_data	bool. This MUST be False if dataset has uns['spatial']['is_single'] == False
        """
        validator = validator_with_spatial_and_is_single_false
        validator.adata.obs["is_primary_data"][0] = True
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: When uns['spatial']['is_single'] is False, obs['is_primary_data'] must be False for all rows."
        ]

    def test_donor_id_must_be_categorical(self, validator_with_adata):
        """
        donor_id categorical with str categories. This MUST be free-text that identifies
        a unique individual that data were derived from. It is STRONGLY RECOMMENDED
        that this identifier be designed so that it is unique to:
        - a given individual within the collection of datasets that includes this dataset
        - a given individual across all collections in the cellxgene Data Portal
        """
        validator = validator_with_adata
        validator.adata.obs["donor_id"] = "NA"
        validator.validate_adata()
        assert validator.errors == ["ERROR: Column 'donor_id' in dataframe 'obs' " "must be categorical, not object."]

    def test_donor_id_must_not_be_empty(self, validator_with_adata):
        validator = validator_with_adata
        obs = validator.adata.obs
        obs["donor_id"] = obs["donor_id"].cat.add_categories("")
        obs["donor_id"].iloc[0] = ""
        validator.validate_adata()
        assert validator.errors == ["ERROR: Column 'donor_id' in dataframe 'obs' " "must not contain empty values."]

    def test_donor_id_must_not_be_nan(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.obs["donor_id"][0] = numpy.nan
        validator.validate_adata()
        assert validator.errors == ["ERROR: Column 'donor_id' in dataframe 'obs' " "must not contain NaN values."]

    @pytest.mark.parametrize(
        "assay,suspension_types",
        {
            "EFO:0010010": ["cell", "nucleus"],
            "EFO:0008720": ["nucleus"],
            "EFO:0008722": ["cell", "nucleus"],
            "EFO:0030002": ["cell"],
            "EFO:0008853": ["cell"],
            "EFO:0030026": ["nucleus"],
            "EFO:0010550": ["cell", "nucleus"],
            "EFO:0008796": ["cell"],
            "EFO:0700003": ["cell"],
            "EFO:0700004": ["cell"],
            "EFO:0008780": ["cell", "nucleus"],
            "EFO:0008953": ["cell"],
            "EFO:0700010": ["cell", "nucleus"],
            "EFO:0700011": ["cell", "nucleus"],
            "EFO:0009919": ["cell", "nucleus"],
            "EFO:0030060": ["cell", "nucleus"],
            "EFO:0022490": ["cell", "nucleus"],
            "EFO:0030028": ["cell", "nucleus"],
            "EFO:0008992": ["na"],
        }.items(),
    )
    def test_suspension_type(self, validator, assay, suspension_types):
        """
        suspension_id categorical with str categories. This field MUST be "cell", "nucleus", or "na". The allowed
        values depend on the assay_ontology_term_id. MUST  fail if the corresponding assay is present in the table, but
        the value of the suspension_type does not match the required value(s) in the table.
        """
        # Resetting validator
        validator.adata = examples.adata.copy()
        validator.errors = []
        validator.warnings = []

        invalid_suspension_type = "na"
        if "na" in suspension_types:
            invalid_suspension_type = "nucleus" if "nucleus" not in suspension_types else "cell"
        obs = validator.adata.obs
        obs["suspension_type"] = invalid_suspension_type
        obs["assay_ontology_term_id"] = assay
        obs["suspension_type"] = obs["suspension_type"].astype("category")
        obs["assay_ontology_term_id"] = obs["assay_ontology_term_id"].astype("category")
        validator.validate_adata()
        assert validator.errors == [
            f"ERROR: Column 'suspension_type' in dataframe 'obs' contains invalid values "
            f"'['{invalid_suspension_type}']'. Values must be one of {suspension_types} when "
            f"'assay_ontology_term_id' is in ['{assay}']"
        ]

    @pytest.mark.parametrize(
        "assay,suspension_types",
        {
            "EFO:0030080": ["cell", "nucleus"],
            "EFO:0007045": ["nucleus"],
            "EFO:0010184": ["cell", "nucleus"],
            "EFO:0008994": ["na"],
            "EFO:0008919": ["cell"],
            "EFO:0002761": ["nucleus"],
        }.items(),
    )
    def test_suspension_type_ancestors_inclusive(self, validator_with_adata, assay, suspension_types):
        """
        suspension_id categorical with str categories. This field MUST be "cell", "nucleus", or "na". The allowed
        values depend on the assay_ontology_term_id. MUST  fail if the corresponding assay is present in the table, but
        the value of the suspension_type does not match the required value(s) in the table.
        """
        validator = validator_with_adata
        obs = validator.adata.obs

        invalid_suspension_type = "na"
        if "na" in suspension_types:
            invalid_suspension_type = "nucleus" if "nucleus" not in suspension_types else "cell"
            obs["suspension_type"] = obs["suspension_type"].cat.remove_unused_categories()
        obs["suspension_type"] = invalid_suspension_type
        obs["assay_ontology_term_id"] = assay
        obs["suspension_type"] = obs["suspension_type"].astype("category")
        obs["assay_ontology_term_id"] = obs["assay_ontology_term_id"].astype("category")
        validator.validate_adata()
        assert validator.errors == [
            f"ERROR: Column 'suspension_type' in dataframe 'obs' contains invalid values "
            f"'['{invalid_suspension_type}']'. Values must be one of {suspension_types} when "
            f"'assay_ontology_term_id' is in ['{assay}']"
        ]

    def test_suspension_type_with_descendant_term_id_failure(self, validator_with_adata):
        """
        suspension_id categorical with str categories. This field MUST be "cell", "nucleus", or "na". The allowed
        values depend on the assay_ontology_term_id. MUST support matching against ancestor term rules if specified.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs["suspension_type"] = "nucleus"
        obs["assay_ontology_term_id"] = "EFO:0022615"  # descendant of EFO:0008994
        obs["suspension_type"] = obs["suspension_type"].astype("category")
        obs["assay_ontology_term_id"] = obs["assay_ontology_term_id"].astype("category")
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Column 'suspension_type' in dataframe 'obs' contains invalid values "
            "'['nucleus']'. Values must be one of ['na'] when "
            "'assay_ontology_term_id' is in ['EFO:0022615']"
        ]

    def test_suspension_type_with_descendant_term_id_success(self, validator_with_adata):
        """
        suspension_id categorical with str categories. This field MUST be "cell", "nucleus", or "na". The allowed
        values depend on the assay_ontology_term_id. MUST support matching against ancestor term rules if specified.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "assay_ontology_term_id"] = "EFO:0008904"  # descendant of EFO:0007045
        obs["suspension_type"][0] = "nucleus"

        validator.validate_adata()
        assert validator.errors == []

    def test_suspension_type_unrecognized_assay(self, validator_with_adata):
        """
        suspension_id categorical with str categories. This field MUST be "cell", "nucleus", or "na". The allowed
        values depend on the assay_ontology_term_id. MUST warn if the corresponding assay is not recognized.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[1], "assay_ontology_term_id"] = "EFO:0700005"
        validator.validate_adata()
        assert validator.errors == []
        assert (
            "WARNING: Data contains assay(s) that are not represented in the 'suspension_type' schema definition table. Ensure you have selected the most appropriate value for the assay(s) between 'cell', 'nucleus', and 'na'. Please contact cellxgene@chanzuckerberg.com during submission so that the assay(s) can be added to the schema definition document."
            in validator.warnings
        )

    def test_categories_with_zero_values_warn(self, validator_with_adata):
        validator = validator_with_adata
        obs = validator.adata.obs
        modified_donor_id = obs["donor_id"].cat.add_categories("donor_3")
        obs["donor_id"] = modified_donor_id
        validator.validate_adata()
        assert (
            "WARNING: Column 'donor_id' in dataframe 'obs' contains a category 'donor_3' with zero observations. These categories will be removed when `--add-labels` flag is present."
            in validator.warnings
        )

    def test_deprecated_fields(self, validator_with_adata):
        validator = validator_with_adata
        obs = validator.adata.obs
        obs["ethnicity"] = "test"
        obs["ethnicity_ontology_term_id"] = "unknown"

        validator.validate_adata()
        assert validator.errors == [
            "ERROR: The field 'ethnicity' is present in 'obs', but it is deprecated.",
            "ERROR: The field 'ethnicity_ontology_term_id' is present in 'obs', but it is deprecated.",
        ]

    def test_fields_with_double_underscore_fail(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.obs["__test_field"] = "test"

        validator.validate_adata()
        assert validator.errors == [
            "ERROR: The field '__test_field' in 'obs' is invalid. Fields that start with '__' are reserved.",
        ]

    def test_obs_column_name_uniqueness(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.obs = pd.concat([validator.adata.obs, validator.adata.obs["suspension_type"]], axis=1)

        validator.validate_adata()
        assert len(validator.errors) > 0

    def test_nan_values_must_be_rejected(self, validator_with_adata):
        """
        NaN values should not be allowed in dataframes
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = numpy.nan
        validator.validate_adata()
        assert (
            "ERROR: Column 'tissue_ontology_term_id' in dataframe 'obs' must not contain NaN values."
            in validator.errors
        )


class TestVar:
    """
    Fail cases in adata.var and adata.raw.var
    """

    def test_var_and_raw_var_same_index(self, validator_with_adata):
        """
        var.index MUST contain unique identifiers for features. raw.var.index MUST be identical to var.index.
        """
        validator = validator_with_adata

        # Swap first row for second one
        var = getattr_anndata(validator.adata, "var")

        # First swap the index
        new_index = list(var.index)
        tmp = new_index[0]
        new_index[0] = new_index[1]
        new_index[1] = tmp
        var.set_index(pd.Index(new_index), inplace=True)

        # Then swap the actual rows
        tmp = var.iloc[0, :].copy()
        var.iloc[0, :] = var.iloc[1, :].copy()
        var.iloc[1, :] = tmp

        validator.validate_adata()
        print("FOO", validator.errors)
        assert validator.errors == ["ERROR: Index of 'raw.var' is not identical to index of 'var'."]

    @pytest.mark.parametrize("component_name", ["var", "raw.var"])
    def test_check_unique_var(self, validator_with_adata, component_name):
        """
        var.index MUST contain unique ENSEMBL gene identifiers for features.
        """
        validator = validator_with_adata
        # Duplicate 1st row in var and assign it to 2nd
        component = getattr_anndata(validator.adata, component_name)
        new_index = list(component.index)
        new_index[1] = new_index[0]
        component.set_index(pd.Index(new_index), inplace=True)
        component.iloc[1, :] = component.iloc[0, :]

        validator.validate_adata()
        assert validator.errors == [f"ERROR: Column 'index' in dataframe '{component_name}' is not unique."]

    def test_column_presence(self, validator_with_adata):
        """
        var is a pandas.DataFrame. Curators MUST annotate the following columns in the obs dataframe.
        feature_is_filtered must not be in raw.var, and it's only checked in var
        """
        validator = validator_with_adata

        column = "feature_is_filtered"
        component_name = "var"
        component = getattr_anndata(validator.adata, component_name)
        component.drop(column, axis=1, inplace=True)

        validator.validate_adata()
        assert validator.errors == [f"ERROR: Dataframe '{component_name}' is missing " f"column " f"'{column}'."]

    def test_feature_is_filtered(self, validator_with_adata):
        """
        feature_is_filtered bool. This MUST be True if the feature was filtered out in the final matrix (X)
        but is present in the raw matrix (raw.X). The value for all cells of the given feature in the
        final matrix MUST be 0.

        Otherwise, this MUST be False.
        """
        validator = validator_with_adata
        var = validator.adata.var
        X = validator.adata.X
        # Duplicate 1st row in var and assigned to 2nd
        var["feature_is_filtered"][0] = True
        for i in range(X.shape[0]):
            X[i, 0] = 0
        X[0, 0] = 1
        validator.adata.X = X.map_blocks(lambda x: (x.eliminate_zeros() or x), dtype=X.dtype, meta=X._meta)
        validator.reset(None, 2)
        validator.validate_adata()
        assert validator.warnings[0] == (
            "WARNING: Some features are 'True' in 'feature_is_filtered' of dataframe 'var', "
            "but there are 1 non-zero values in the corresponding columns of the matrix 'X'. "
            "All values for these features must be 0."
        )

        # Test that feature_is_filtered is a bool and not a string
        var["feature_is_filtered"] = "string"
        validator.reset(None, 2)
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Column 'feature_is_filtered' in dataframe 'var' must be boolean, not 'object'."
        ]

    def test_feature_is_filtered_should_be_true(self, validator_with_adata):
        validator = validator_with_adata
        var = validator.adata.var
        X = validator.adata.X

        # zero out first column of X, but set feature_is_filtered to False
        var["feature_is_filtered"][0] = False
        for i in range(X.shape[0]):
            X[i, 0] = 0
        validator.adata.X = X.map_blocks(lambda x: (x.eliminate_zeros() or x), dtype=X.dtype, meta=X._meta)

        validator.reset(None, 2)
        validator.validate_adata()
        assert (
            validator.warnings[0]
            == "WARNING: Genes 'ENSG00000127603' have all-zero values in adata.X. Either feature_is_filtered should be set to True or adata.raw.X should be set to all-zero values."
        )

    def test_feature_is_filtered_var_mishapen(self, validator_with_adata):
        validator = validator_with_adata
        gene_index = 0
        gene_name = validator_with_adata.adata.raw.var_names[gene_index]
        var_to_keep = [v for v in validator_with_adata.adata.raw.var_names if v != gene_name]
        raw_adata = anndata.AnnData(
            X=validator_with_adata.adata.raw.X,
            obs=validator_with_adata.adata.obs,
            var=validator_with_adata.adata.raw.var,
        )
        new_raw_adata = raw_adata[:, var_to_keep]
        validator_with_adata.adata.raw = new_raw_adata
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            "ERROR: Could not complete full validation of feature_is_filtered because of size differences between var and raw.var."
            in validator.errors
        )

    def test_columns_not_in_raw_var(self, validator_with_adata):
        """
        Curators MUST annotate the following column only in the var dataframe.
        This column MUST NOT be present in raw.var:
            feature_is_filtered
        """
        validator = validator_with_adata

        validator.adata.raw = validator.adata
        validator.validate_adata()
        assert validator.errors == ["ERROR: Column 'feature_is_filtered' must not be present in 'raw.var'."]

    @pytest.mark.parametrize("component_name", ["var", "raw.var"])
    def test_feature_id_wrong_format(self, validator_with_adata, component_name):
        """
        feature_id (var.index) str.

        This tests the case of an ID with an incorrect format "ENSEBML_NOGENE"
        """
        validator = validator_with_adata
        component = getattr_anndata(validator.adata, component_name)

        new_index = list(component.index)
        new_index[0] = "ENSEBML_NOGENE"
        component.set_index(pd.Index(new_index), inplace=True)

        validator.validate_adata()
        assert validator.errors == [
            f"ERROR: Could not infer organism from feature ID 'ENSEBML_NOGENE' "
            f"in '{component_name}', make sure it is a valid ID."
        ]

    @pytest.mark.parametrize("component_name", ["var", "raw.var"])
    def test_feature_id_non_existent_ensembl(self, validator_with_adata, component_name):
        """
        feature_id (var.index) str.

        This tests the case of an ENSEMBL ID that has the right format but doesn't exist
        """
        validator = validator_with_adata
        component = getattr_anndata(validator.adata, component_name)

        new_index = list(component.index)
        new_index[0] = "ENSG000"
        component.set_index(pd.Index(new_index), inplace=True)

        validator.validate_adata()
        assert len(validator.errors) > 0

    @pytest.mark.parametrize("component_name", ["var", "raw.var"])
    def test_feature_id_non_existent_ercc(self, validator_with_adata, component_name):
        """
        feature_id (var.index) str.

        This tests the case of an ERCC ID that has the right format but doesn't exist
        """
        validator = validator_with_adata
        component = getattr_anndata(validator.adata, component_name)

        new_index = list(component.index)
        new_index[0] = "ERCC-000000"
        component.set_index(pd.Index(new_index), inplace=True)

        validator.validate_adata()
        assert len(validator.errors) > 0

    def test_should_warn_for_low_gene_count(self, validator_with_adata):
        """
        Raise a warning if there are too few genes
        """
        validator = validator_with_adata
        # NOTE:[EM] changing the schema def here is stateful and results in unpredictable test results.
        #  Reset after mutating.
        _old_schema = deepcopy(validator.schema_def.copy())

        validator.schema_def["components"]["var"]["warn_if_less_than_rows"] = 100
        validator.validate_adata()
        assert (
            f"WARNING: Dataframe 'var' only has {NUMBER_OF_GENES} rows. Features SHOULD NOT be filtered from expression matrix."
            in validator.warnings
        )
        validator.schema_def = _old_schema

    @pytest.mark.parametrize(
        "df,column",
        [
            (df, label["to_column"])
            for df in ["var", "raw.var"]
            for label in schema_def["components"][df]["index"]["add_labels"]
        ],
    )
    def test_add_label_fields_are_reserved(self, validator_with_adata, df, column):
        """
        Raise an error if column names flagged as 'add_label' -> 'to_column' in the schema definition are not available.
        """
        validator = validator_with_adata
        component = getattr_anndata(validator.adata, df)
        component[column] = "dummy_value"

        validator.validate_adata()
        assert validator.errors == [
            f"ERROR: Add labels error: Column '{column}' is a reserved column name "
            f"of '{df}'. Remove it from h5ad and try again."
        ]

    def test_var_column_name_uniqueness(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.var = pd.concat([validator.adata.var, validator.adata.var["feature_is_filtered"]], axis=1)

        validator.validate_adata()
        assert len(validator.errors) > 0

    def test_raw_var_column_name_uniqueness(self, validator_with_adata):
        validator = validator_with_adata
        original_var = validator.adata.var.copy()
        validator.adata.var = pd.concat([validator.adata.var, validator.adata.var["feature_is_filtered"]], axis=1)
        validator.adata.raw = validator.adata
        validator.adata.var = original_var  # Ensure only the raw.var has duplicate columns

        validator.validate_adata()
        assert len(validator.errors) > 0


class TestUns:
    """
    Fail cases in adata.uns
    """

    @pytest.mark.parametrize("reserved_column", schema_def["components"]["uns"]["reserved_columns"])
    def test_reserved_columns_presence(self, validator_with_adata, reserved_column):
        """
        Reserved columns must NOT be used in uns
        """
        validator = validator_with_adata
        validator.adata.uns[reserved_column] = "dummy_value"
        validator.validate_adata()
        assert validator.errors == [
            f"ERROR: Column '{reserved_column}' is a reserved column name "
            f"of 'uns'. Remove it from h5ad and try again."
        ]

    def test_required_fields_title(self, validator_with_adata):
        """
        Curators MUST annotate `schema_version` and values in uns (title)
        """
        validator = validator_with_adata

        del validator.adata.uns["title"]
        validator.validate_adata()
        assert validator.errors == ["ERROR: 'title' in 'uns' is not present."]

    def test_leading_trailing_double_spaces_in_strings(self, validator_with_adata):
        """
        The following sequences MUST NOT appear in str types documented in the schema:
            Leading control or space separators - ”     This is an example”
            Trailing control or space separators - “This is an example     ”
            Multiple (internal) control or space separators - "This     is an example"
        """
        validator = validator_with_adata
        uns = validator.adata.uns
        uns["title"] = " There is a leading space"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: ' There is a leading space' in 'uns['title']' is not valid, it contains leading spaces."
        ]

        uns["title"] = "There is a trailing space "
        validator.errors = []
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'There is a trailing space ' in 'uns['title']' is not valid, it contains trailing spaces."
        ]

        uns["title"] = "There are   double   spaces"
        validator.errors = []
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'There are   double   spaces' in 'uns['title']' is not valid, it contains double spaces."
        ]

    def test_title(self, validator_with_adata):
        """
        Title MUST be a string
        """
        validator = validator_with_adata

        # list instead of string
        validator.adata.uns["title"] = ["title"]
        validator.validate_adata()
        assert validator.errors == ["ERROR: '['title']' in 'uns['title']' is not valid, " "it must be a string."]

    def test_batch_condition_is_list(self, validator_with_adata):
        """
        batch_condition list[str]
        """
        validator = validator_with_adata
        uns = validator.adata.uns

        # Check valid case of numpy array which is interchangeable with lists
        uns["batch_condition"] = numpy.array(uns["batch_condition"])
        validator.validate_adata()
        assert validator.errors == []

        # Check fail case: not a list nor numpy array
        uns["batch_condition"] = "cell_type_ontology_term_id"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'cell_type_ontology_term_id' in 'uns['batch_condition']' "
            "is not valid, it must be a list or numpy array."
        ]

    def test_batch_condition_is_column_from_obs(self, validator_with_adata):
        """
        batch_condition list[str]. str values MUST refer to cell metadata keys in obs.
        """
        validator = validator_with_adata

        validator.adata.uns["batch_condition"] = ["NO_COLUMN"]
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Value 'NO_COLUMN' of list 'batch_condition' is not a " "column in 'adata.obs'."
        ]

    def test_default_embedding_is_str(self, validator_with_adata):
        """
        Default_embedding str.
        """
        validator = validator_with_adata

        validator.adata.uns["default_embedding"] = ["X_umap"]
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: '['X_umap']' in 'uns['default_embedding']' is not valid, " "it must be a string."
        ]

    def test_default_embedding_is_key_from_obsm(self, validator_with_adata):
        """
        Default_embedding str. The value MUST match a key to an embedding in obsm
        """
        validator = validator_with_adata

        validator.adata.uns["default_embedding"] = "X_other"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'X_other' in 'uns['default_embedding']' is not valid, " "it must be a key of 'adata.obsm'."
        ]

    def test_X_approximate_distribution_is_str(self, validator_with_adata):
        """
        X_approximate_distribution str. The value MUST be "count" [...] or "normal".
        Note that `normal` is tested in the happy path test case using `good_uns`.
        """
        validator = validator_with_adata
        uns = validator.adata.uns
        # Check valid case of "count" which is not included in valid object
        uns["X_approximate_distribution"] = "count"
        validator.validate_adata()
        assert validator.errors == []

        # Invalid type: list
        uns["X_approximate_distribution"] = ["count"]
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: '['count']' in 'uns['X_approximate_distribution']' " "is not valid, it must be a string."
        ]

    def test_X_approximate_distribution_is_valid(self, validator_with_adata):
        """
        X_approximate_distribution str. The value MUST be "count" [...] or "normal"
        """
        validator = validator_with_adata

        validator.adata.uns["X_approximate_distribution"] = "COUNT"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'COUNT' in 'uns['X_approximate_distribution']' is " "not valid. Allowed terms: ['count', 'normal']."
        ]

    def test_deprecated_fields(self, validator_with_adata):
        validator = validator_with_adata
        uns = validator.adata.uns
        uns["X_normalization"] = "test_value"
        uns["default_field"] = "test_value"
        uns["layer_descriptions"] = "test_value"
        uns["tags"] = "test_value"
        uns["version"] = "test_value"
        uns["contributors"] = "test_value"
        uns["preprint_doi"] = "test_value"
        uns["project_description"] = "test_value"
        uns["project_links"] = "test_value"
        uns["project_name"] = "test_value"
        uns["publication_doi"] = "test_value"

        validator.validate_adata()
        assert validator.errors == [
            "ERROR: The field 'X_normalization' is present in 'uns', but it is deprecated.",
            "ERROR: The field 'default_field' is present in 'uns', but it is deprecated.",
            "ERROR: The field 'layer_descriptions' is present in 'uns', but it is deprecated.",
            "ERROR: The field 'tags' is present in 'uns', but it is deprecated.",
            "ERROR: The field 'version' is present in 'uns', but it is deprecated.",
            "ERROR: The field 'contributors' is present in 'uns', but it is deprecated.",
            "ERROR: The field 'preprint_doi' is present in 'uns', but it is deprecated.",
            "ERROR: The field 'project_description' is present in 'uns', but it is deprecated.",
            "ERROR: The field 'project_links' is present in 'uns', but it is deprecated.",
            "ERROR: The field 'project_name' is present in 'uns', but it is deprecated.",
            "ERROR: The field 'publication_doi' is present in 'uns', but it is deprecated.",
        ]

    def test_colors_happy_path(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata = examples.adata_with_colors.copy()
        assert validator.validate_adata()

    def test_colors_happy_path_no_column_def(self, validator_with_adata):
        """
        Creates a new column in `obs[test_column]` where the dtype is categorical, but its column definition
        is not in the schema. Having valid colors in `uns[test_column_colors]` should not raise an error, since
        we should be able to find the corresponding `obs` column with the `category` dtype.
        """
        validator = validator_with_adata
        validator.adata.obs["test_column"] = ["one", "two"]
        validator.adata.obs["test_column"] = validator.adata.obs["test_column"].astype("category")
        validator.adata.uns["test_column_colors"] = numpy.array(["#000000", "#ffffff"])
        assert validator.validate_adata()

    def test_uns_true_value(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["log1p"] = True
        assert validator.validate_adata()

    def test_uns_false_value(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["log1p"] = False
        assert validator.validate_adata()

    def test_uns_none_value(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["log1p"] = None
        assert validator.validate_adata()

    def test_uns_empty_numpy_array(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["log1p"] = numpy.array([])
        validator.validate_adata()
        assert validator.errors == ["ERROR: uns['log1p'] cannot be an empty value."]

    def test_uns_empty_python_array(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["log1p"] = []
        validator.validate_adata()
        assert validator.errors == ["ERROR: uns['log1p'] cannot be an empty value."]

    def test_uns_empty_string(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["log1p"] = ""
        validator.validate_adata()
        assert validator.errors == ["ERROR: uns['log1p'] cannot be an empty value."]

    def test_uns_empty_dictionary(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["log1p"] = {}
        validator.validate_adata()
        assert validator.errors == ["ERROR: uns['log1p'] cannot be an empty value."]

    def test_uns_bool_allowed(self, validator_with_adata):
        validator = validator_with_adata

        # Regular bool value
        validator.adata.uns["log1p"] = True
        validator.validate_adata()
        assert validator.errors == []

        # Numpy bool_ value
        validator.adata.uns["log1p"] = numpy.bool_(True)
        validator.validate_adata()
        assert validator.errors == []

        # Numpy bool value
        validator.adata.uns["log1p"] = numpy.bool_(True)
        validator.validate_adata()
        assert validator.errors == []

    @pytest.mark.parametrize("matrix_factory", [scipy.sparse.csr_matrix, scipy.sparse.csc_matrix])
    def test_uns_scipy_matrices_cannot_be_empty(self, validator_with_adata, matrix_factory):
        validator: Validator = validator_with_adata

        validator.adata.uns["test"] = from_array(matrix_factory([[1]], dtype=int))
        validator.validate_adata()
        assert validator.errors == []
        validator.reset()

        validator.adata.uns["test"] = from_array(matrix_factory([[]], dtype=int))
        validator.validate_adata()
        assert validator.errors == ["ERROR: uns['test'] cannot be an empty value."]

    def test_uns_numbers_are_allowed(self, validator_with_adata):
        validator = validator_with_adata

        validator.adata.uns["numpy.int64"] = numpy.int64(10)
        validator.adata.uns["int"] = 1
        validator.adata.uns["int_zero"] = 0
        validator.adata.uns["float"] = float(4)
        validator.adata.uns["numpy.float32"] = numpy.float32(3)
        validator.validate_adata()
        assert validator.errors == []

    def test_colors_happy_path_duplicates(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["suspension_type_colors"] = numpy.array(["lightgrey", "lightgrey"])
        assert validator.validate_adata()

    def test_colors_not_numpy_array(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["suspension_type_colors"] = ["green", "purple"]
        validator.adata.uns["donor_id_colors"] = {"donor_1": "#000000", "donor_2": "#ffffff"}
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Colors field uns['suspension_type_colors'] must be of 'numpy.ndarray' type, it is <class 'list'>",
            "ERROR: Colors field uns['donor_id_colors'] must be of 'numpy.ndarray' type, it is <class 'dict'>",
        ]

    def test_colors_bool_obs_counterpart(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["is_primary_data_colors"] = numpy.array(["green", "purple"])
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Colors field uns[is_primary_data_colors] does not have a corresponding categorical field in obs. is_primary_data is present but is dtype bool"
        ]

    def test_colors_portal_obs_counterpart(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["cell_type_colors"] = numpy.array(["green", "purple"])
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Colors field uns[cell_type_colors] does not have a corresponding categorical field in obs. Annotate cell_type_ontology_term_id_colors instead"
        ]

    def test_colors_without_obs_counterpart(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["fake_field_colors"] = numpy.array(["green", "purple"])
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Colors field uns[fake_field_colors] does not have a corresponding categorical field in obs"
        ]

    def test_empty_color_array(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["suspension_type_colors"] = numpy.array([])
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: uns['suspension_type_colors'] cannot be an empty value.",
            "ERROR: Annotated categorical field suspension_type must have at least 1 color options in uns[suspension_type_colors]. Found: []",
        ]

    def test_different_color_types(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["suspension_type_colors"] = numpy.array(["#000000", "pink"])
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Colors in uns[suspension_type_colors] must be either all hex colors or all CSS4 named colors. Found: ['#000000' 'pink']"
        ]

    def test_invalid_hex_color_option(sel, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["suspension_type_colors"] = numpy.array(["#000", "#ffffff"])
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Colors in uns[suspension_type_colors] must be either all hex colors or all CSS4 named colors. Found: ['#000' '#ffffff']"
        ]

    def test_invalid_named_color_option(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["suspension_type_colors"] = numpy.array(["green", "pynk"])
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Colors in uns[suspension_type_colors] must be either all hex colors or all CSS4 named colors. Found: ['green' 'pynk']"
        ]

    def test_invalid_named_color_option_empty_strings(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["suspension_type_colors"] = numpy.array(["", "green"])
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Colors in uns[suspension_type_colors] must be either all hex colors or all CSS4 named colors. "
            "Found: ['' 'green']"
        ]

    def test_invalid_named_color_option_integer(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["suspension_type_colors"] = numpy.array([3, 4])
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Colors in uns[suspension_type_colors] must be strings. Found: [3 4] which are int64"
        ]

    def test_invalid_named_color_option_none(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["suspension_type_colors"] = numpy.array(["green", None])
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Colors in uns[suspension_type_colors] must be strings. Found: ['green' None] which are object"
        ]

    def test_invalid_named_color_option_nan(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["suspension_type_colors"] = numpy.array([numpy.nan, numpy.nan])
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Colors in uns[suspension_type_colors] must be strings. Found: [nan nan] which are float64"
        ]


class TestObsm:
    """
    Fail cases for adata.obsm
    """

    @pytest.mark.parametrize("key", ["X_tsne", "spatial"])
    def test_obsm_values_ara_numpy(self, validator_with_visium_assay, key):
        """
        values in obsm MUST be a numpy.ndarray
        """
        validator = validator_with_visium_assay
        obsm = validator.adata.obsm
        obsm[key] = pd.DataFrame(obsm["X_umap"], index=validator.adata.obs_names)
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: All embeddings have to be of 'numpy.ndarray' type, "
            f"'adata.obsm['{key}']' is <class 'pandas.core.frame.DataFrame'>')."
        ]

    @pytest.mark.parametrize("key", ["X_umap", "spatial"])
    def test_obsm_values_infinity(self, validator_with_visium_assay, key):
        """
        values in obsm cannot have any infinity values
        """
        validator = validator_with_visium_assay
        validator.adata.obsm[key][0:100, 1] = numpy.inf
        validator.validate_adata()
        assert validator.errors == [
            f"ERROR: adata.obsm['{key}'] contains positive infinity or negative infinity values."
        ]

    @pytest.mark.parametrize("key", ["X_umap", "spatial"])
    def test_obsm_values_str(self, validator_with_visium_assay, key):
        """
        values in obsm must be numerical types, strings are not valid
        """
        validator = validator_with_visium_assay
        obsm = validator.adata.obsm
        all_string = numpy.full(obsm[key].shape, "test")
        obsm[key] = all_string
        validator.validate_adata()
        assert validator.errors == [
            f"ERROR: adata.obsm['{key}'] has an invalid data type. It should be float, integer, or unsigned "
            "integer of any precision (8, 16, 32, or 64 bits)."
        ]

    @pytest.mark.parametrize("key", ["X_umap", "spatial"])
    def test_obsm_values_nan(self, validator_with_visium_assay, key):
        """
        test obsm NaN restrictions for different embedding types.
        feature embeddings: X_* cannot be all NaN
        spatial emeddings: 'spatial' cannot have any NaNs
        """
        validator = validator_with_visium_assay
        obsm = validator.adata.obsm

        # Check embedding has any NaN
        obsm[key][0:100, 1] = numpy.nan
        validator.reset(None, 2)
        validator.validate_adata()

        if key != "spatial":
            assert validator.errors == []
        else:
            assert validator.errors == ["ERROR: adata.obsm['spatial'] contains at least one NaN value."]

        # Check embedding has all NaNs
        all_nan = numpy.full(obsm[key].shape, numpy.nan)
        obsm[key] = all_nan
        validator.reset(None, 2)
        validator.validate_adata()
        if key != "spatial":
            assert validator.errors == [f"ERROR: adata.obsm['{key}'] contains all NaN values."]
        else:
            assert validator.errors == ["ERROR: adata.obsm['spatial'] contains at least one NaN value."]

    def test_obsm_values_no_X_embedding__non_spatial_dataset(self, validator_with_adata):
        """
        X_{suffix} embeddings MUST exist for non-spatial datasets
        """
        validator: Validator = validator_with_adata
        validator.adata.obsm["harmony"] = validator.adata.obsm["X_umap"]
        validator.adata.uns["default_embedding"] = "harmony"
        del validator.adata.obsm["X_umap"]
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: At least one embedding in 'obsm' has to have a key with an 'X_' prefix.",
        ]
        assert validator.is_spatial is False
        assert (
            "WARNING: Embedding key in 'adata.obsm' harmony is not 'spatial' nor does it start with 'X_'. Thus, it will not be available in Explorer"
            in validator.warnings
        )

    @pytest.mark.parametrize("assay_ontology_term_id", ["EFO:0022859", "EFO:0030062", "EFO:0022860"])
    def test_obsm_values_no_X_embedding__visium_dataset(self, validator_with_visium_assay, assay_ontology_term_id):
        """
        X_{suffix} embeddings MAY exist for spatial datasets
        """
        validator: Validator = validator_with_visium_assay
        validator.adata.uns["default_embedding"] = "spatial"
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id

        # may have X_{suffix} embedding
        validator._validate_obsm()
        assert validator.is_spatial is True
        assert validator.errors == []
        validator.reset()

        # may also have no X_{suffix} embedding
        del validator.adata.obsm["X_umap"]
        validator._validate_obsm()
        assert validator.is_spatial is True
        assert validator.errors == []
        validator.reset()

    def test_obsm_values_no_X_embedding__slide_seq_v2_dataset(self, validator_with_slide_seq_v2_assay):
        validator = validator_with_slide_seq_v2_assay
        validator.adata.uns["default_embedding"] = "spatial"
        del validator.adata.obsm["X_umap"]
        validator.validate_adata()
        assert validator.errors == []
        assert validator.is_spatial is True

    def test_obsm_values_spatial_embedding_missing__is_single_true(self, validator_with_visium_assay):
        validator = validator_with_visium_assay
        del validator.adata.obsm["spatial"]
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'spatial' embedding is required in 'adata.obsm' if adata.uns['spatial']['is_single'] is True."
        ]

    def test_obsm_values_spatial_embedding_missing__is_single_false(self, validator_with_spatial_and_is_single_false):
        validator = validator_with_spatial_and_is_single_false
        del validator.adata.obsm["spatial"]
        validator.validate_adata()
        assert validator.errors == []

    def test_obsm_values_spatial_embedding_present__is_single_none(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.obsm["spatial"] = validator.adata.obsm["X_umap"]
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'spatial' embedding is forbidden in 'adata.obsm' if "
            "adata.uns['spatial']['is_single'] is not set."
        ]

    def test_obsm_values_warn_start_with_X(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.obsm["harmony"] = pd.DataFrame(validator.adata.obsm["X_umap"], index=validator.adata.obs_names)
        validator.validate_adata()
        assert (
            "WARNING: Embedding key in 'adata.obsm' harmony is not 'spatial' nor does it start with 'X_'. Thus, it will not be available in Explorer"
            in validator.warnings
        )
        assert validator.errors == [
            "ERROR: All embeddings have to be of 'numpy.ndarray' type, 'adata.obsm['harmony']' is <class 'pandas.core.frame.DataFrame'>')."
        ]

    def test_obsm_values_suffix_is_forbidden(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.obsm["X_spatial"] = validator.adata.obsm["X_umap"]
        validator.validate_adata()
        assert validator.errors == ["ERROR: Embedding key in 'adata.obsm' X_spatial cannot be used."]

    def test_obsm_values_suffix_start_with_number(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.obsm["X_3D"] = pd.DataFrame(validator.adata.obsm["X_umap"], index=validator.adata.obs_names)
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Suffix for embedding key in 'adata.obsm' X_3D does not match the regex pattern ^[a-zA-Z][a-zA-Z0-9_.-]*$.",
            "ERROR: All embeddings have to be of 'numpy.ndarray' type, 'adata.obsm['X_3D']' is <class 'pandas.core.frame.DataFrame'>').",
        ]

    def test_obsm_values_key_start_with_number(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.obsm["3D"] = pd.DataFrame(validator.adata.obsm["X_umap"], index=validator.adata.obs_names)
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Embedding key in 'adata.obsm' 3D does not match the regex pattern ^[a-zA-Z][a-zA-Z0-9_.-]*$.",
            "ERROR: All embeddings have to be of 'numpy.ndarray' type, 'adata.obsm['3D']' is <class "
            "'pandas.core.frame.DataFrame'>').",
        ]
        assert (
            "WARNING: Embedding key in 'adata.obsm' 3D is not 'spatial' nor does it start with 'X_'. Thus, it will not be available in Explorer"
            in validator.warnings
        )

    def test_obsm_suffix_name_valid(self, validator_with_adata):
        """
        Suffix after X_ must be at least 1 character long
        """
        validator = validator_with_adata
        validator.adata.obsm["X_"] = validator.adata.obsm["X_umap"]
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Suffix for embedding key in 'adata.obsm' X_ does not match the regex pattern ^[a-zA-Z][a-zA-Z0-9_.-]*$."
        ]

    def test_obsm_key_name_whitespace(self, validator_with_adata):
        """
        Embedding keys beginning with X_ with whitespace are not valid
        """
        validator = validator_with_adata
        obsm = validator.adata.obsm
        obsm["X_ umap"] = obsm["X_umap"]
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Suffix for embedding key in 'adata.obsm' X_ umap does not match the regex pattern ^[a-zA-Z][a-zA-Z0-9_.-]*$.",
        ]

        del obsm["X_ umap"]
        obsm["u m a p"] = obsm["X_umap"]
        validator.reset(None, 2)
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Embedding key in 'adata.obsm' u m a p does not match the regex pattern ^[a-zA-Z][a-zA-Z0-9_.-]*$."
        ]

    def test_obsm_suffix_has_special_characters_valid(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.obsm["X_umap_MinDist_0.2_N_Neighbors-15"] = validator.adata.obsm["X_umap"]
        validator.validate_adata()
        assert validator.errors == []

    @pytest.mark.parametrize("key", ["X_umap", "spatial"])
    def test_obsm_shape_one_column(self, validator_with_visium_assay, key):
        """
        Curators MUST annotate one or more two-dimensional (m >= 2) embeddings
        """
        # Makes 1 column array
        validator = validator_with_visium_assay
        validator.adata.obsm[key] = numpy.delete(validator.adata.obsm[key], 0, 1)
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: All 'X_' and 'spatial' embeddings must have at least two columns. "
            f"'adata.obsm['{key}']' has columns='1'."
        ]

    def test_obsm_shape_zero_column_with_unknown_key(self, validator_with_adata):
        """
        embeddings that are not 'X_' or 'spatial' that are ndarrays must have at least one column
        """
        # Makes 0 column array
        validator = validator_with_adata
        n_obs = validator_with_adata.adata.n_obs
        validator.adata.obsm["unknown"] = numpy.zeros((n_obs, 0))
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: The size of the ndarray stored for a 'adata.obsm['unknown']' MUST NOT " "be zero.",
            "ERROR: All unspecified embeddings must have at least one column. "
            "'adata.obsm['unknown']' has columns='0'.",
        ]

    def test_obsm_shape_same_rows_and_columns(self, validator_with_adata):
        """
        The number of rows must be equal to the number of columns
        """
        validator = validator_with_adata
        obsm = validator.adata.obsm
        # Create a 3 row array
        arr1 = numpy.array([0, 0])
        arr2 = numpy.array([0, 0])
        arr3 = numpy.array([0, 0])
        three_row_array = numpy.vstack((arr1, arr2, arr3))
        with pytest.raises(ValueError):
            obsm["X_umap"] = three_row_array
            validator.validate_adata()

    def test_obsm_size_zero(self, validator_with_adata):
        """
        size of obsm an key cannot be zero.
        """
        validator = validator_with_adata
        adata = validator.adata
        adata.obsm["badsize"] = numpy.empty((2, 0))
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: The size of the ndarray stored for a 'adata.obsm['badsize']' MUST NOT " "be zero.",
            "ERROR: All unspecified embeddings must have at least one column. "
            "'adata.obsm['badsize']' has columns='0'.",
        ]


class TestObsp:
    def test_obsp_size_zero(self, validator_with_adata):
        """
        size of obsp an key cannot be zero.
        """
        validator = validator_with_adata
        adata = validator.adata
        adata.obsp["badsize"] = numpy.empty((2, 2, 0))
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: The size of the ndarray stored for a 'adata.obsp['badsize']' MUST NOT be zero."
        ]


class TestVarm:
    def test_varm_size_zero(self, validator_with_adata):
        """
        size of varm an key cannot be zero.
        """
        validator = validator_with_adata
        adata = validator.adata
        adata.varm["badsize"] = numpy.empty((NUMBER_OF_GENES, 0))
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: The size of the ndarray stored for a 'adata.varm['badsize']' MUST NOT be " "zero."
        ]


class TestVarp:
    def test_varp_size_zero(self, validator_with_adata):
        """
        size of varp an key cannot be zero.
        """
        validator = validator_with_adata
        adata = validator.adata
        adata.varp["badsize"] = numpy.empty((NUMBER_OF_GENES, NUMBER_OF_GENES, 0))
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: The size of the ndarray stored for a 'adata.varp['badsize']' MUST NOT be zero."
        ]


class TestAddingLabels:
    """
    Tests the addition of labels from IDs based on schema specification. The test is done by comparing manually
    created dataframes (positive control) against the ones produced by the validator
    """

    @pytest.mark.parametrize(
        "column", [label["to_column"] for label in schema_def["components"]["var"]["index"]["add_labels"]]
    )
    def test_var_added_labels(self, label_writer, adata_with_labels, column):
        """
        When a dataset is uploaded, cellxgene Data Portal MUST automatically add the matching human-readable
        name for the corresponding feature identifier and the inferred NCBITaxon term for the reference organism
        to the var dataframe. Curators MUST NOT annotate the columns below:
        """
        expected_column = adata_with_labels.var[column]
        obtained_column = label_writer.adata.var[column]

        for i, j in zip(expected_column.tolist(), obtained_column.tolist()):
            assert i == j

    @pytest.mark.parametrize(
        "column",
        [
            "assay",
            "cell_type",
            "development_stage",
            "disease",
            "self_reported_ethnicity",
            "sex",
            "tissue",
        ],
    )
    def test_obs_added_labels(self, label_writer, adata_with_labels, column):
        """
        When a dataset is uploaded, the cellxgene Data Portal MUST automatically add the matching human-readable
        name for the corresponding ontology term to the obs dataframe.
        Curators MUST NOT annotate the following columns.

            - assay. categorical with str categories. This MUST be the human-readable name assigned to the value
            of assay_ontology_term_id.
            - cell_type. categorical with str categories. This MUST be the human-readable name assigned to the value
            of cell_type_ontology_term_id or "unknown" if set in cell_type_ontology_term_id.
            - development_stage. categorical with str categories. This MUST be "unknown" if set in
            development_stage_ontology_term_id; otherwise, this MUST be the human-readable name assigned to
            the value of development_stage_ontology_term_id.
            - disease. categorical with str categories. This MUST be the human-readable name assigned to
            the value of disease_ontology_term_id.
            - self_reported_ethnicity. categorical with str categories. This MUST be "na" or "unknown" if
            set in self_reported_ethnicity_ontology_term_id; otherwise, this MUST be the human-readable
            name assigned to the value of self_reported_ethnicity_ontology_term_id.
            - sex. categorical with str categories. This MUST be "unknown" if set in sex_ontology_term_id;
            otherwise, this MUST be the human-readable name assigned to the value of sex_ontology_term_id.
            - tissue. categorical with str categories. This MUST be the human-readable name assigned to the
            value of tissue_ontology_term_id. " (cell culture)" or " (organoid)" MUST
            be appended if present in tissue_ontology_term_id.
        """
        expected_column = adata_with_labels.obs[column]
        obtained_column = label_writer.adata.obs[column]

        for i, j in zip(expected_column.tolist(), obtained_column.tolist()):
            assert i == j

    def test_uns_added_labels(self, label_writer, adata_with_labels):
        expected_uns = adata_with_labels.uns
        obtained_uns = label_writer.adata.uns

        assert expected_uns == obtained_uns

    def test_obs_added_tissue_type_label__unknown(self):
        adata = examples.adata.copy()
        obs = adata.obs

        # Arrange
        obs.at["Y", "tissue_type"] = "cell culture"  # Already set in example data, just setting explicitly here
        obs.at["Y", "tissue_ontology_term_id"] = "unknown"  # Testing this term case
        labeler = AnnDataLabelAppender(adata)
        labeler._add_labels()  # Annotate

        assert labeler.adata.obs.at["Y", "tissue"] == "unknown"

    def test_obs_added_cell_type_label__unknown(self):
        adata = examples.adata.copy()
        obs = adata.obs

        # Arrange
        obs.at["Y", "cell_type_ontology_term_id"] = "unknown"  # Testing this term case
        labeler = AnnDataLabelAppender(adata)
        labeler._add_labels()  # Annotate

        assert labeler.adata.obs.at["Y", "cell_type"] == "unknown"

    def test_remove_unused_categories(self, label_writer, adata_with_labels):
        modified_donor_id = label_writer.adata.obs["donor_id"].cat.add_categories("donor_2")
        label_writer.adata.obs["donor_id"] = modified_donor_id
        case = unittest.TestCase()
        case.assertCountEqual(label_writer.adata.obs["donor_id"].dtype.categories, ["donor_1", "donor_2"])
        label_writer._remove_categories_with_zero_values()
        case.assertCountEqual(label_writer.adata.obs["donor_id"].dtype.categories, ["donor_1"])


class TestZebrafish:
    """
    Tests for the zebrafish schema
    """

    @pytest.fixture
    def zebrafish_obs(self):
        obs = examples.adata.copy().obs
        for i in range(2):
            obs.loc[obs.index[i], "cell_type_ontology_term_id"] = "ZFA:0000003"
            obs.loc[obs.index[i], "development_stage_ontology_term_id"] = "ZFS:0000016"
            obs.loc[obs.index[i], "self_reported_ethnicity_ontology_term_id"] = "na"
            obs.loc[obs.index[i], "tissue_ontology_term_id"] = "ZFA:0001262"
        return obs

    @pytest.fixture
    def zebrafish_visium_obs(self):
        obs = examples.adata_visium.copy().obs
        for i in range(2):
            obs.loc[obs.index[i], "cell_type_ontology_term_id"] = "unknown"
            obs.loc[obs.index[i], "development_stage_ontology_term_id"] = "ZFS:0000016"
            obs.loc[obs.index[i], "self_reported_ethnicity_ontology_term_id"] = "na"
            obs.loc[obs.index[i], "tissue_ontology_term_id"] = "ZFA:0001262"
        return obs

    @pytest.fixture
    def zebrafish_var(self):
        return pd.DataFrame(
            [[False]],
            index=[
                "ENSDARG00000103202",
                "ENSDARG00000009657",
                "ENSDARG00000096472",
                "ENSDARG00000096156",
                "ENSDARG00000076160",
                "ENSDARG00000117163",
                "ENSDARG00000096187",
            ],
            columns=["feature_is_filtered"],
        )

    @pytest.fixture
    def validator_with_zebrafish_adata(self, validator_with_adata, zebrafish_obs, zebrafish_var):
        validator_with_adata.adata.obs = zebrafish_obs
        validator_with_adata.adata.uns["organism_ontology_term_id"] = "NCBITaxon:7955"
        validator_with_adata.adata.var = zebrafish_var.copy()
        new_raw = anndata.AnnData(
            X=validator_with_adata.adata.raw.X.copy(), var=zebrafish_var.copy(), obs=zebrafish_obs.copy()
        )
        new_raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        validator_with_adata.adata.raw = new_raw
        return validator_with_adata

    @pytest.fixture
    def validator_with_visium_zebrafish_adata(self, validator_with_visium_assay, zebrafish_visium_obs, zebrafish_var):
        validator_with_visium_assay.adata.obs = zebrafish_visium_obs
        validator_with_visium_assay.adata.uns["organism_ontology_term_id"] = "NCBITaxon:7955"
        validator_with_visium_assay.adata.var = zebrafish_var.copy()
        new_raw = anndata.AnnData(
            X=validator_with_visium_assay.adata.raw.X.copy(), var=zebrafish_var.copy(), obs=zebrafish_visium_obs.copy()
        )
        new_raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        validator_with_visium_assay.adata.raw = new_raw
        return validator_with_visium_assay

    @pytest.mark.parametrize(
        "development_stage_ontology_term_id",
        ["ZFS:0000016", "unknown"],
    )
    def test_development_stage_ontology_term_id_zebrafish(
        self, validator_with_zebrafish_adata, development_stage_ontology_term_id
    ):
        """
        If organism_ontology_term_id is "NCBITaxon:7955" for Danio rerio,
        this MUST be the most accurate ZFS:0100000 descendant or "unknown" and MUST NOT be ZFS:0000000.
        """
        validator = validator_with_zebrafish_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = development_stage_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    @pytest.mark.parametrize(
        "development_stage_ontology_term_id,error",
        [
            (
                "HsapDv:0000001",  # Wrong ontology
                "ERROR: 'HsapDv:0000001' in 'development_stage_ontology_term_id' is not a valid ontology term id of "
                "'ZFA'.",
            ),
            (
                "ZFA:0000001",  # Same ontology, not a descendant of ZFS:0100000
                "ERROR: 'ZFA:0000001' in 'development_stage_ontology_term_id' is not an allowed term id.",
            ),
            (
                "ZFS:0100000",  # Do not accept ZFS:0100000 itself, must be a descendant
                "ERROR: 'ZFS:0100000' in 'development_stage_ontology_term_id' is not an allowed term id.",
            ),
            (
                "ZFS:0000000",  # Descendant of ZFS:0100000 but explicitly forbidden term
                "ERROR: 'ZFS:0000000' in 'development_stage_ontology_term_id' is not allowed.",
            ),
        ],
    )
    def test_development_stage_ontology_term_id_zebrafish__invalid(
        self, validator_with_zebrafish_adata, development_stage_ontology_term_id, error
    ):
        """
        If organism_ontology_term_id is "NCBITaxon:7955" for Danio rerio,
        this MUST be the most accurate ZFS:0100000 descendant or "unknown" and MUST NOT be ZFS:0000000.
        """
        zebrafish_error_message_suffix = (
            "When 'organism_ontology_term_id' is 'NCBITaxon:7955' (Danio rerio), "
            "'development_stage_ontology_term_id' MUST be the most accurate descendant of 'ZFS:0100000' and it "
            "MUST NOT be 'ZFS:0000000' for Unknown. The str 'unknown' is acceptable."
        )
        validator = validator_with_zebrafish_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = development_stage_ontology_term_id
        validator.validate_adata()
        assert validator.errors == [error + " " + zebrafish_error_message_suffix]

    @pytest.mark.parametrize(
        "cell_type_ontology_term_id",
        ["ZFA:0000003", "CL:4023077", "unknown"],
    )
    def test_cell_type_ontology_term_id(self, validator_with_zebrafish_adata, cell_type_ontology_term_id):
        """
        If organism_ontology_term_id is "NCBITaxon:7955" for Danio rerio,
        MUST be a descendant term id of 'ZFA:0009000' (cell) or 'unknown'
        """
        validator = validator_with_zebrafish_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "cell_type_ontology_term_id"] = cell_type_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    @pytest.mark.parametrize(
        "cell_type_ontology_term_id",
        [
            "UBERON:0000001",  # Wrong ontology
            "ZFA:0001094",  # Same ontology, not a descendant of ZFA:0009000
            "ZFA:0009000",  # Do not accept ZFA:0009000 itself, must be a descendant
            "na",  # Allowed for other organisms, not allowed if organism is zebrafih
        ],
    )
    def test_cell_type_ontology_term_id__invalid(self, validator_with_zebrafish_adata, cell_type_ontology_term_id):
        validator = validator_with_zebrafish_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "cell_type_ontology_term_id"] = cell_type_ontology_term_id
        validator.validate_adata()
        assert len(validator.errors) > 0

    def test_organism_cell_type_ontology_term_id__visium_in_tissue_0(self, validator_with_visium_zebrafish_adata):
        validator: Validator = validator_with_visium_zebrafish_adata
        obs = validator.adata.obs
        for i in range(2):
            obs.loc[obs.index[i], "in_tissue"] = i
            obs.loc[obs.index[i], "cell_type_ontology_term_id"] = "unknown"
        validator.reset(None, 2)
        validator.validate_adata()
        assert not validator.errors

    def test_organism_cell_type_ontology_term_id__visium_in_tissue_0_invalid(
        self, validator_with_visium_zebrafish_adata
    ):
        validator = validator_with_visium_zebrafish_adata
        obs = validator.adata.obs
        for i in range(2):
            obs.loc[obs.index[i], "in_tissue"] = 0
            obs.loc[obs.index[i], "cell_type_ontology_term_id"] = "ZFA:0000003"
        validator.validate_adata()
        assert (
            f"obs['cell_type_ontology_term_id'] must be 'unknown' when {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_IN_TISSUE_0}"
            in validator.errors[0]
        )

    @pytest.mark.parametrize(
        "tissue_ontology_term_id",
        [
            "ZFA:0001262",  # valid descendant of ZFA:0100000
            "UBERON:0002048",  # valid UBERON term
        ],
    )
    def test_organism_tissue_type_ontology_term_id(self, validator_with_zebrafish_adata, tissue_ontology_term_id):
        validator = validator_with_zebrafish_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = tissue_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    @pytest.mark.parametrize(
        "tissue_ontology_term_id",
        [
            "CL:0000001",  # Wrong ontology
            "ZFS:0000016",  # Same ontology, not a descendant of ZFA:0100000
            "ZFA:0100000",  # Must be descendant of ZFA:0100000, not itself
            "ZFA:0009000",  # ZFA:0009000 is an explicitly forbidden term
            "ZFA:0000003",  # ZFA:0009000 descendant, an explicitly forbidden ancestor
            "na",
            "unknown",
        ],
    )
    def test_tissue_ontology_term_id__invalid(self, validator_with_zebrafish_adata, tissue_ontology_term_id):
        validator = validator_with_zebrafish_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = tissue_ontology_term_id
        validator.validate_adata()
        assert len(validator.errors) > 0

    @pytest.mark.parametrize(
        "tissue_type",
        ["tissue", "cell culture", "organoid"],
    )
    def test_organism_tissue_type_valid(self, validator_with_zebrafish_adata, tissue_type):
        validator = validator_with_zebrafish_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = tissue_type
        assert not validator.errors

    @pytest.mark.parametrize(
        "tissue_ontology_term_id",
        [
            "CL:4023077",  # valid CL term for cell culture
            "ZFA:0000003",  # valid ZFA term
        ],
    )
    def test_cell_culture_tissue_ontology_term_id(self, validator_with_zebrafish_adata, tissue_ontology_term_id):
        validator = validator_with_zebrafish_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = "cell culture"
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = tissue_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    def test_cell_culture_tissue_ontology_term_id_invalid(self, validator_with_zebrafish_adata):
        validator = validator_with_zebrafish_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = "cell culture"
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = (
            "UBERON:0002048"  # only valid UBERON term if not cell culture
        )
        validator.validate_adata()
        assert len(validator.errors) > 0


class TestFruitFly:
    """
    Tests for the fruit fly schema
    """

    @pytest.fixture
    def fruitfly_obs(self):
        obs = examples.adata.copy().obs
        for i in range(2):
            obs.loc[obs.index[i], "cell_type_ontology_term_id"] = "FBbt:00049192"
            obs.loc[obs.index[i], "self_reported_ethnicity_ontology_term_id"] = "na"
            obs.loc[obs.index[i], "development_stage_ontology_term_id"] = "FBdv:00005370"
            obs.loc[obs.index[i], "tissue_ontology_term_id"] = "FBbt:00007337"
        return obs

    @pytest.fixture
    def fruitfly_visium_obs(self):
        obs = examples.adata_visium.copy().obs
        for i in range(2):
            obs.loc[obs.index[i], "cell_type_ontology_term_id"] = "unknown"
            obs.loc[obs.index[i], "self_reported_ethnicity_ontology_term_id"] = "na"
            obs.loc[obs.index[i], "development_stage_ontology_term_id"] = "FBdv:00005370"
            obs.loc[obs.index[i], "tissue_ontology_term_id"] = "FBbt:00007337"
        return obs

    @pytest.fixture
    def fruitfly_var(self):
        return pd.DataFrame(
            [[False]],
            index=[
                "FBgn0038542",
                "RR42085_transposable_element",
                "FBgn0039592",
                "FBgn0038067",
                "FBgn0053534",
                "FBgn0039602",
                "FBgn0264897",
            ],
            columns=["feature_is_filtered"],
        )

    @pytest.fixture
    def validator_with_fruitfly_adata(self, validator_with_adata, fruitfly_obs, fruitfly_var):
        validator_with_adata.adata.obs = fruitfly_obs
        validator_with_adata.adata.uns["organism_ontology_term_id"] = "NCBITaxon:7227"
        validator_with_adata.adata.var = fruitfly_var.copy()
        new_raw = anndata.AnnData(
            X=validator_with_adata.adata.raw.X.copy(), var=fruitfly_var.copy(), obs=fruitfly_obs.copy()
        )
        new_raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        validator_with_adata.adata.raw = new_raw
        return validator_with_adata

    @pytest.fixture
    def validator_with_visium_fruitfly_adata(self, validator_with_visium_assay, fruitfly_visium_obs, fruitfly_var):
        validator_with_visium_assay.adata.obs = fruitfly_visium_obs
        validator_with_visium_assay.adata.uns["organism_ontology_term_id"] = "NCBITaxon:7227"
        validator_with_visium_assay.adata.var = fruitfly_var.copy()
        new_raw = anndata.AnnData(
            X=validator_with_visium_assay.adata.raw.X.copy(), var=fruitfly_var.copy(), obs=fruitfly_visium_obs.copy()
        )
        new_raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        validator_with_visium_assay.adata.raw = new_raw
        return validator_with_visium_assay

    @pytest.mark.parametrize(
        "development_stage_ontology_term_id",
        [
            "FBdv:00007117",  # descendant of FBdv:00007014 for adult age in days
            "FBdv:00005370",  # descendant of FBdv:00005259 for developmental stage
            "unknown",
        ],
    )
    def test_development_stage_ontology_term_id_fruitfly(
        self, validator_with_fruitfly_adata, development_stage_ontology_term_id
    ):
        validator = validator_with_fruitfly_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = development_stage_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    @pytest.mark.parametrize(
        "development_stage_ontology_term_id",
        [
            "HsapDv:0000001",  # Wrong ontology
            "FBdv:00007012",  # Explicitly forbidden term, life stage
        ],
    )
    def test_development_stage_ontology_term_id_fruitfly__invalid(
        self, validator_with_fruitfly_adata, development_stage_ontology_term_id
    ):
        validator = validator_with_fruitfly_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = development_stage_ontology_term_id
        validator.validate_adata()
        assert len(validator.errors) > 0

    @pytest.mark.parametrize(
        "cell_type_ontology_term_id",
        ["FBbt:00049192", "CL:4023077", "unknown"],
    )
    def test_cell_type_ontology_term_id(self, validator_with_fruitfly_adata, cell_type_ontology_term_id):
        """
        If organism_ontology_term_id is "NCBITaxon:7227" for Drosophila melanogaster,
        MUST be a descendant term id of 'FBbt:0007002' (cell) or 'unknown'
        """
        validator = validator_with_fruitfly_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "cell_type_ontology_term_id"] = cell_type_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    @pytest.mark.parametrize(
        "cell_type_ontology_term_id",
        [
            "UBERON:0000001",  # Wrong ontology
            "FBbt:00007001",  # Same ontology, not a descendant of FBbt:00007002
            "FBbt:00007002",  # Do not accept FBbt:00007002 itself, must be a descendant
            "na",  # Allowed for other organisms, not allowed if organism is fruit fly
        ],
    )
    def test_cell_type_ontology_term_id__invalid(self, validator_with_fruitfly_adata, cell_type_ontology_term_id):
        validator = validator_with_fruitfly_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "cell_type_ontology_term_id"] = cell_type_ontology_term_id
        validator.validate_adata()
        assert len(validator.errors) > 0

    def test_organism_cell_type_ontology_term_id__visium_in_tissue_0(self, validator_with_visium_fruitfly_adata):
        validator = validator_with_visium_fruitfly_adata
        obs = validator.adata.obs
        for i in range(2):
            obs.loc[obs.index[i], "in_tissue"] = i
            obs.loc[obs.index[i], "cell_type_ontology_term_id"] = "unknown"
        validator.reset(None, 2)
        validator.validate_adata()
        assert not validator.errors

    def test_organism_cell_type_ontology_term_id__visium_in_tissue_0_invalid(
        self, validator_with_visium_fruitfly_adata
    ):
        validator: Validator = validator_with_visium_fruitfly_adata
        obs = validator.adata.obs
        for i in range(2):
            obs.loc[obs.index[i], "in_tissue"] = i
            obs.loc[obs.index[i], "cell_type_ontology_term_id"] = "FBbt:00049192"
        validator.reset(None, 2)
        validator.validate_adata()
        assert (
            f"obs['cell_type_ontology_term_id'] must be 'unknown' when {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_IN_TISSUE_0}."
            in validator.errors[0]
        )

    @pytest.mark.parametrize(
        "tissue_ontology_term_id",
        [
            "FBbt:00007337",  # valid descendant of FBbt:10000000
            "UBERON:0002048",  # valid UBERON term
        ],
    )
    def test_organism_tissue_type_ontology_term_id(self, validator_with_fruitfly_adata, tissue_ontology_term_id):
        validator = validator_with_fruitfly_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = tissue_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    @pytest.mark.parametrize(
        "tissue_ontology_term_id",
        [
            "CL:0000001",  # Wrong ontology
            "FBbt:10000000",  # Must be descendant of FBbt:10000000, not itself
            "FBbt:00007002",  # FBbt:00007002 is an explicitly forbidden term
            "FBbt:00007294",  # FBbt:00007002 descendant, an explicitly forbidden ancestor
            "na",
            "unknown",
        ],
    )
    def test_tissue_ontology_term_id__invalid(self, validator_with_fruitfly_adata, tissue_ontology_term_id):
        validator = validator_with_fruitfly_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = tissue_ontology_term_id
        validator.validate_adata()
        assert len(validator.errors) > 0

    @pytest.mark.parametrize(
        "tissue_type",
        ["tissue", "cell culture", "organoid"],
    )
    def test_organism_tissue_type_valid(self, validator_with_fruitfly_adata, tissue_type):
        validator = validator_with_fruitfly_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = tissue_type
        assert not validator.errors

    @pytest.mark.parametrize(
        "tissue_ontology_term_id",
        [
            "CL:4023077",  # valid CL term for cell culture
            "FBbt:00049192",  # valid FBbt term
        ],
    )
    def test_cell_culture_tissue_ontology_term_id(self, validator_with_fruitfly_adata, tissue_ontology_term_id):
        validator = validator_with_fruitfly_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = "cell culture"
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = tissue_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    def test_cell_culture_tissue_ontology_term_id_invalid(self, validator_with_fruitfly_adata):
        validator = validator_with_fruitfly_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = "cell culture"
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = (
            "UBERON:0002048"  # only valid UBERON term if not cell culture
        )
        validator.validate_adata()
        assert len(validator.errors) > 0


class TestRoundworm:
    """
    Tests for the roundworm / c. elegans schema
    """

    @pytest.fixture
    def roundworm_obs(self):
        obs = examples.adata.copy().obs
        for i in range(2):
            obs.loc[obs.index[i], "cell_type_ontology_term_id"] = "WBbt:0008611"
            obs.loc[obs.index[i], "self_reported_ethnicity_ontology_term_id"] = "na"
            obs.loc[obs.index[i], "development_stage_ontology_term_id"] = "WBls:0000532"
            obs.loc[obs.index[i], "tissue_ontology_term_id"] = "WBbt:0006749"
            obs.loc[obs.index[i], "sex_ontology_term_id"] = "PATO:0000384"
        return obs

    @pytest.fixture
    def roundworm_visium_obs(self):
        obs = examples.adata_visium.copy().obs
        for i in range(2):
            obs.loc[obs.index[i], "self_reported_ethnicity_ontology_term_id"] = "na"
            obs.loc[obs.index[i], "development_stage_ontology_term_id"] = "WBls:0000532"
            obs.loc[obs.index[i], "tissue_ontology_term_id"] = "WBbt:0006749"
            obs.loc[obs.index[i], "sex_ontology_term_id"] = "PATO:0000384"
            obs.loc[obs.index[i], "cell_type_ontology_term_id"] = "unknown"
        return obs

    @pytest.fixture
    def roundworm_var(self):
        return pd.DataFrame(
            [[False]],
            index=[
                "WBGene00000003",
                "WBGene00000007",
                "WBGene00000014",
                "WBGene00000015",
                "WBGene00000022",
                "WBGene00000024",
                "WBGene00000027",
            ],
            columns=["feature_is_filtered"],
        )

    @pytest.fixture
    def validator_with_roundworm_adata(self, validator_with_adata, roundworm_obs, roundworm_var):
        validator_with_adata.adata.obs = roundworm_obs
        validator_with_adata.adata.uns["organism_ontology_term_id"] = "NCBITaxon:6239"
        validator_with_adata.adata.var = roundworm_var.copy()
        new_raw = anndata.AnnData(
            X=validator_with_adata.adata.raw.X.copy(), var=roundworm_var.copy(), obs=roundworm_obs.copy()
        )
        new_raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        validator_with_adata.adata.raw = new_raw
        return validator_with_adata

    @pytest.fixture
    def validator_with_visium_roundworm_adata(self, validator_with_visium_assay, roundworm_visium_obs, roundworm_var):
        validator_with_visium_assay.adata.obs = roundworm_visium_obs
        validator_with_visium_assay.adata.uns["organism_ontology_term_id"] = "NCBITaxon:6239"
        validator_with_visium_assay.adata.var = roundworm_var.copy()
        new_raw = anndata.AnnData(
            X=validator_with_visium_assay.adata.raw.X.copy(), var=roundworm_var.copy(), obs=roundworm_visium_obs.copy()
        )
        new_raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        validator_with_visium_assay.adata.raw = new_raw
        return validator_with_visium_assay

    @pytest.mark.parametrize(
        "development_stage_ontology_term_id",
        [
            "WBls:0000669",  # unfertilized egg Ce
            "WBls:0000805",  # descendant of WBls:0000803
            "WBls:0000816",  # descendant of WBls:0000804
            "unknown",
        ],
    )
    def test_development_stage_ontology_term_id_roundworm(
        self, validator_with_roundworm_adata, development_stage_ontology_term_id
    ):
        """
        If organism_ontology_term_id is "NCBITaxon:6239" for C. elegans,
        this MUST be the most accurate WBls term or 'unknown'
        """
        validator = validator_with_roundworm_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = development_stage_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    @pytest.mark.parametrize(
        "development_stage_ontology_term_id",
        [
            "HsapDv:0000001",  # Wrong ontology
            "WBls:0000825",  # Not a descendant of WBls:0000803 or WBls:0000804
        ],
    )
    def test_development_stage_ontology_term_id_roundworm__invalid(
        self, validator_with_roundworm_adata, development_stage_ontology_term_id
    ):
        validator = validator_with_roundworm_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = development_stage_ontology_term_id
        validator.validate_adata()
        assert len(validator.errors) > 0

    @pytest.mark.parametrize(
        "cell_type_ontology_term_id",
        ["WBbt:0005762", "CL:4023077", "unknown"],
    )
    def test_cell_type_ontology_term_id(self, validator_with_roundworm_adata, cell_type_ontology_term_id):
        """
        If organism_ontology_term_id is "NCBITaxon:6239" for C. elegans,
        MUST be a descendant term id of 'WBbt:0004017' (cell) or 'unknown'
        """
        validator = validator_with_roundworm_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "cell_type_ontology_term_id"] = cell_type_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    @pytest.mark.parametrize(
        "cell_type_ontology_term_id",
        [
            "UBERON:0000001",  # Wrong ontology
            "WBbt:0000100",  # Same ontology, not a descendant of WBbt:0004017
            "WBbt:0004017",  # Do not accept WBbt:0004017 itself, must be a descendant
            "na",  # Allowed for other organisms, not allowed if organism is fruit fly
        ],
    )
    def test_cell_type_ontology_term_id__invalid(self, validator_with_roundworm_adata, cell_type_ontology_term_id):
        validator = validator_with_roundworm_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "cell_type_ontology_term_id"] = cell_type_ontology_term_id
        validator.validate_adata()
        assert len(validator.errors) > 0

    def test_organism_cell_type_ontology_term_id__visium_in_tissue_0(self, validator_with_visium_roundworm_adata):
        validator: Validator = validator_with_visium_roundworm_adata
        obs = validator.adata.obs
        for i in range(2):
            obs.loc[obs.index[i], "in_tissue"] = i
            obs.loc[obs.index[i], "cell_type_ontology_term_id"] = "unknown"
        validator.reset(None, 2)
        validator.validate_adata()
        assert not validator.errors

    def test_organism_cell_type_ontology_term_id__visium_in_tissue_0_invalid(
        self, validator_with_visium_roundworm_adata
    ):
        validator: Validator = validator_with_visium_roundworm_adata
        obs = validator.adata.obs
        for i in range(2):
            obs.loc[obs.index[i], "in_tissue"] = i
            obs.loc[obs.index[i], "cell_type_ontology_term_id"] = "WBbt:0005739"
        validator.reset(None, 2)
        validator.validate_adata()
        assert (
            f"obs['cell_type_ontology_term_id'] must be 'unknown' when {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_IN_TISSUE_0}"
            in validator.errors[0]
        )

    @pytest.mark.parametrize(
        "tissue_ontology_term_id",
        [
            "WBbt:0006750",  # valid descendant of WBbt:0005766
            "UBERON:0002048",  # valid UBERON term
        ],
    )
    def test_organism_tissue_type_ontology_term_id(self, validator_with_roundworm_adata, tissue_ontology_term_id):
        validator = validator_with_roundworm_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = tissue_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    @pytest.mark.parametrize(
        "tissue_ontology_term_id",
        [
            "CL:0000001",  # Wrong ontology
            "WBbt:0005766",  # Anatomy, explicitly forbidden term - must be a descendant of this term
            "WBbt:0007849",  # hermaphrodite, explicitly forbidden term
            "WBbt:0007850",  # male, explicitly forbidden term
            "WBbt:0008595",  # female, explicitly forbidden term
            "WBbt:0004017",  # cell, explicitly forbidden term
            "WBbt:0008611",  # descendant of WBbt:0004017 (cell)
            "WBbt:00006803",  # nucleus, explicitly forbidden term
            "WBbt:0002702",  # descendant of WBbt:00006803 (nucleus)
            "na",
            "unknown",
        ],
    )
    def test_tissue_ontology_term_id__invalid(self, validator_with_roundworm_adata, tissue_ontology_term_id):
        validator = validator_with_roundworm_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = tissue_ontology_term_id
        validator.validate_adata()
        assert len(validator.errors) > 0

    @pytest.mark.parametrize(
        "tissue_type",
        ["tissue", "cell culture", "organoid"],
    )
    def test_organism_tissue_type_valid(self, validator_with_roundworm_adata, tissue_type):
        validator = validator_with_roundworm_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = tissue_type
        assert not validator.errors

    @pytest.mark.parametrize(
        "sex_ontology_term_id",
        ["unknown", "PATO:0000384", "PATO:0001340"],
    )
    def test_sex_ontology_term_id_valid(self, validator_with_roundworm_adata, sex_ontology_term_id):
        validator = validator_with_roundworm_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "sex_ontology_term_id"] = sex_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    def test_sex_ontology_term_id__invalid(self, validator_with_roundworm_adata):
        validator = validator_with_roundworm_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "sex_ontology_term_id"] = "PATO:0000383"  # allowed for other organisms, not c. elegans
        validator.validate_adata()
        error_message = (
            "ERROR: 'PATO:0000383' in 'sex_ontology_term_id' is not an allowed term id. When "
            "'organism_ontology_term_id' is 'NCBITaxon:6239' (Caenorhabditis elegans), "
            "'sex_ontology_term_id' MUST be 'PATO:0000384' for male, 'PATO:0001340' for hermaphrodite, or 'unknown'."
        )
        assert error_message in validator.errors

    @pytest.mark.parametrize(
        "tissue_ontology_term_id",
        [
            "CL:4023077",  # valid CL term for cell culture
            "WBbt:0005762",  # valid WBbt term
        ],
    )
    def test_cell_culture_tissue_ontology_term_id(self, validator_with_roundworm_adata, tissue_ontology_term_id):
        validator = validator_with_roundworm_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = "cell culture"
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = tissue_ontology_term_id
        validator.validate_adata()
        assert not validator.errors

    def test_cell_culture_tissue_ontology_term_id_invalid(self, validator_with_roundworm_adata):
        validator = validator_with_roundworm_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = "cell culture"
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = "UBERON:0002048"
        validator.validate_adata()
        assert len(validator.errors) > 0


class TestMultiSpecies:
    """
    Tests to verify our support for human / mouse is not impacted by support for additional species
    """

    @pytest.mark.parametrize(
        "cell_type_ontology_term_id",
        [
            "UBERON:0000001",  # Wrong ontology
            "ZFA:0000003",  # Valid for zebrafish, not valid for human or mouse data
            "FBbt:00049192",  # Valid for fruit fly, not valid for human or mouse data
            "WBbt:0008611",  # Valid for roundworm, not valid for human or mouse data
            "na",  # Allowed for other organisms, not allowed if organism is fruit fly
        ],
    )
    def test_cell_type_ontology_term_id__invalid(self, validator_with_adata, cell_type_ontology_term_id):
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "cell_type_ontology_term_id"] = cell_type_ontology_term_id
        validator.validate_adata()
        assert len(validator.errors) > 0

    @pytest.mark.parametrize(
        "tissue_ontology_term_id",
        [
            "CL:0000001",  # Wrong ontology
            "ZFA:0001262",  # Valid for zebrafish, not valid for human or mouse data
            "FBbt:00007337",  # Valid for fruit fly, not valid for human or mouse data
            "WBbt:0006749",  # Valid for roundworm, not valid for human or mouse data
            "na",
            "unknown",
        ],
    )
    def test_tissue_ontology_term_id__invalid(self, validator_with_adata, tissue_ontology_term_id):
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = tissue_ontology_term_id
        validator.validate_adata()
        assert len(validator.errors) > 0
