"""
Tests for schema compliance of an AnnData object
"""

import tempfile
import unittest

import anndata
import fixtures.examples_validate as examples
import numpy
import pandas as pd
import pytest
import scipy.sparse
from cellxgene_schema.schema import get_schema_definition
from cellxgene_schema.utils import getattr_anndata
from cellxgene_schema.validate import (
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE,
    VISIUM_AND_IS_SINGLE_TRUE_MATRIX_SIZE,
    Validator,
)
from cellxgene_schema.write_labels import AnnDataLabelAppender

schema_def = get_schema_definition()


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
    validator.visium_and_is_single_true_matrix_size = 2
    return validator


@pytest.fixture
def validator_with_slide_seq_v2_assay(validator) -> Validator:
    validator.adata = examples.adata_slide_seqv2.copy()
    return validator


@pytest.fixture
def label_writer(validator_with_validated_adata) -> AnnDataLabelAppender:
    """
    Fixture that returns an AnnDataLabelAppender object
    """
    label_writer = AnnDataLabelAppender(validator_with_validated_adata)
    label_writer._add_labels()
    return label_writer


def save_and_read_adata(adata: anndata.AnnData) -> anndata.AnnData:
    """
    Saves adata to a temporary file and reads it back. Used to test read/write errors.
    :param adata: AnnData object
    :return: AnnData object
    """
    with tempfile.NamedTemporaryFile(suffix=".h5ad") as f:
        adata.write_h5ad(f.name)
        return anndata.read_h5ad(f.name)


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
        assert "ERROR: Number of genes in X (3) is different than raw.X (4)." in validator.errors

    def test_sparsity(self, validator_with_adata):
        """
        In any layer, if a matrix has 50% or more values that are zeros, it is STRONGLY RECOMMENDED that
        the matrix be encoded as a scipy.sparse.csr_matrix
        """
        validator = validator_with_adata
        sparse_X = numpy.zeros([validator.adata.obs.shape[0], validator.adata.var.shape[0]], dtype=numpy.float32)
        sparse_X[0, 1] = 1
        sparse_X[1, 1] = 1
        validator.adata.X = sparse_X
        validator.validate_adata()
        assert validator.warnings == [
            "WARNING: Sparsity of 'X' is 0.75 which is greater than 0.5, "
            "and it is not a 'scipy.sparse.csr_matrix'. It is "
            "STRONGLY RECOMMENDED to use this type of matrix for "
            "the given sparsity."
        ]

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
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: All non-zero values in raw matrix must be positive integers of type numpy.float32.",
            "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements.",
        ]

    @pytest.mark.parametrize("datatype", [int, "float16", "float64"])
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

        validator = validator_with_visium_assay
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

        validator = validator_with_visium_assay
        validator.adata.X[1] = numpy.zeros(validator.adata.var.shape[0], dtype=numpy.float32)
        validator.adata.raw.X[1] = numpy.zeros(validator.adata.var.shape[0], dtype=numpy.float32)
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
        validator.adata.X = numpy.zeros(
            [validator.adata.obs.shape[0], validator.adata.var.shape[0]], dtype=numpy.float32
        )
        validator.adata.raw = validator.adata.copy()
        validator.adata.raw.var.drop("feature_is_filtered", axis=1, inplace=True)
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
        validator.validate_adata()
        assert validator.errors == []

    def test_raw_values__invalid_visium_and_is_single_true_row_length(self, validator_with_visium_assay):
        """
        Dataset is visium and uns['is_single'] is True, but raw.X is the wrong length.
        """
        validator = validator_with_visium_assay
        validator.visium_and_is_single_true_matrix_size = VISIUM_AND_IS_SINGLE_TRUE_MATRIX_SIZE

        validator.validate_adata()
        assert validator.errors == [
            f"ERROR: When {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE}, the raw matrix must be the "
            "unfiltered feature-barcode matrix 'raw_feature_bc_matrix'. It must have exactly "
            f"{validator.visium_and_is_single_true_matrix_size} rows. Raw matrix row count is 2.",
            "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements.",
        ]

    def test_raw_values__multiple_invalid_in_tissue_errors(self, validator_with_visium_assay):
        """
        Dataset is visium and uns['is_single'] is True, in_tissue has both 0 and 1 values and there
        are issues validating rows of both in the matrix.
        """

        validator = validator_with_visium_assay
        validator.visium_and_is_single_true_matrix_size = VISIUM_AND_IS_SINGLE_TRUE_MATRIX_SIZE
        validator.adata.X = numpy.zeros(
            [validator.adata.obs.shape[0], validator.adata.var.shape[0]], dtype=numpy.float32
        )
        validator.adata.raw = validator.adata.copy()
        validator.adata.raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        validator.validate_adata()
        assert validator.errors == [
            f"ERROR: When {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE}, the raw matrix must be the "
            "unfiltered feature-barcode matrix 'raw_feature_bc_matrix'. It must have exactly "
            f"{validator.visium_and_is_single_true_matrix_size} rows. Raw matrix row count is 2.",
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
        with unittest.mock.patch.object(validator_with_adata._chunk_matrix, "__defaults__", (1,)):
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
        assert validator.warnings == [
            "WARNING: Only raw data was found, i.e. there is no 'raw.X'. "
            "It is STRONGLY RECOMMENDED that 'final' (normalized) data is provided."
        ]


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

    def test_column_presence_organism(self, validator_with_adata):
        """
        obs is a pandas.DataFrame. Curators MUST annotate the following columns in the obs dataframe.

        A separate check is need for organism_ontology_term_id because removing from anndata results in multiple
        errors given that other columns depend on its presence
        """
        validator = validator_with_adata
        validator.adata.obs.drop("organism_ontology_term_id", axis=1, inplace=True)
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Dataframe 'obs' is missing column " "'organism_ontology_term_id'.",
            "ERROR: Checking values with dependencies failed for "
            "adata.obs['self_reported_ethnicity_ontology_term_id'], this is likely due "
            "to missing dependent column in adata.obs.",
            "ERROR: Checking values with dependencies failed for "
            "adata.obs['development_stage_ontology_term_id'], this is likely due "
            "to missing dependent column in adata.obs.",
        ]

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
            "ERROR: Checking values with dependencies failed for "
            "adata.obs['suspension_type'], this is likely due "
            "to missing dependent column in adata.obs.",
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

    def test_cell_type_ontology_term_id_invalid_term(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.obs.loc[validator.adata.obs.index[0], "cell_type_ontology_term_id"] = "EFO:0000001"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'EFO:0000001' in 'cell_type_ontology_term_id' is not a valid ontology term id of 'CL'.",
        ]

    @pytest.mark.parametrize(
        "term",
        schema_def["components"]["obs"]["columns"]["cell_type_ontology_term_id"]["curie_constraints"]["forbidden"][
            "terms"
        ],
    )
    def test_cell_type_ontology_term_id(self, validator_with_adata, term):
        """
        cell_type_ontology_term_id categorical with str categories. This MUST be a CL term, and must NOT match forbidden
        columns defined in schema
        """
        validator = validator_with_adata
        validator.adata.obs.loc[validator.adata.obs.index[0], "cell_type_ontology_term_id"] = term
        validator.validate_adata()
        # Forbidden columns may be marked as either "not allowed" or "deprecated"
        assert validator.errors == [
            f"ERROR: '{term}' in 'cell_type_ontology_term_id' is not allowed."
        ] or validator.errors == [f"ERROR: '{term}' in 'cell_type_ontology_term_id' is a deprecated term id of 'CL'."]

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
        If organism_ontolology_term_id is "NCBITaxon:9606" for Homo sapiens,
        this MUST be the most accurate HsapDv:0000001 descendant.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "organism_ontology_term_id"] = "NCBITaxon:9606"
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
        self, validator_with_adata, development_stage_ontology_term_id, error
    ):
        """
        If organism_ontolology_term_id is "NCBITaxon:10090" for Mus musculus,
        this MUST be the most accurate MmusDv:0000001 descendant.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "organism_ontology_term_id"] = "NCBITaxon:10090"
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = development_stage_ontology_term_id
        obs.loc[
            obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "na"
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
        obs.loc[obs.index[0], "organism_ontology_term_id"] = "NCBITaxon:10114"
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = "EFO:0000001"
        obs.loc[
            obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "na"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'EFO:0000001' in 'development_stage_ontology_term_id' is "
            "not a valid ontology term id of 'UBERON'. When 'organism_ontology_term_id' is not 'NCBITaxon:10090' "
            "nor 'NCBITaxon:9606', 'development_stage_ontology_term_id' MUST be a descendant term id of "
            "'UBERON:0000105' excluding 'UBERON:0000071', or unknown."
        ]

        # All other it MUST be descendants of UBERON:0000105 and not UBERON:0000071
        # Fail case UBERON:0000071
        validator.errors = []
        obs.loc[obs.index[0], "organism_ontology_term_id"] = "NCBITaxon:10114"
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = "UBERON:0000071"
        obs.loc[
            obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "na"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'UBERON:0000071' in 'development_stage_ontology_term_id' is not allowed. When "
            "'organism_ontology_term_id' is not 'NCBITaxon:10090' "
            "nor 'NCBITaxon:9606', 'development_stage_ontology_term_id' MUST be a descendant term id of "
            "'UBERON:0000105' excluding 'UBERON:0000071', or unknown.",
        ]

    def test_disease_ontology_term_id(self, validator_with_adata):
        """
        disease_ontology_term_id categorical with str categories. This MUST be one of:
        - PATO:0000461 for normal or healthy
        - descendant of MONDO:0000001 for disease
        - self or descendant of MONDO:0021178 for injury
        """
        validator = validator_with_adata
        obs = validator.adata.obs

        # Invalid ontology
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "EFO:0000001"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'EFO:0000001' in 'disease_ontology_term_id' is not a valid ontology term id of 'MONDO, PATO'. "
            "Only 'PATO:0000461' (normal), 'MONDO:0021178' (injury) or descendant terms thereof, or descendant terms of 'MONDO:0000001' (disease) are allowed"
        ]

        # Invalid PATO term id
        validator.errors = []
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "PATO:0001894"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'PATO:0001894' in 'disease_ontology_term_id' is not an allowed term id. "
            "Only 'PATO:0000461' (normal), 'MONDO:0021178' (injury) or descendant terms thereof, or descendant terms of 'MONDO:0000001' (disease) are allowed"
        ]

        # Invalid MONDO term id - disease characteristic
        validator.errors = []
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "MONDO:0021125"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'MONDO:0021125' in 'disease_ontology_term_id' is not an allowed term id. "
            "Only 'PATO:0000461' (normal), 'MONDO:0021178' (injury) or descendant terms thereof, or descendant terms of 'MONDO:0000001' (disease) are allowed"
        ]

        # Invalid MONDO term id - disease parent term
        validator.errors = []
        obs.loc[obs.index[0], "disease_ontology_term_id"] = "MONDO:0000001"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'MONDO:0000001' in 'disease_ontology_term_id' is not an allowed term id. "
            "Only 'PATO:0000461' (normal), 'MONDO:0021178' (injury) or descendant terms thereof, or descendant terms of 'MONDO:0000001' (disease) are allowed"
        ]

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
        assert validator.errors == [
            "ERROR: 'unknown' in 'tissue_ontology_term_id' is not a valid ontology term id of 'UBERON'. "
            "When 'tissue_type' is 'tissue' or 'organoid', 'tissue_ontology_term_id' MUST be a descendant "
            "term id of 'UBERON:0001062' (anatomical entity)."
        ]

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
        ] = "HANCESTRO:0005,HANCESTRO:0014,unknown"
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
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'HANCESTRO:0306' in 'self_reported_ethnicity_ontology_term_id' is not allowed. Descendant terms "
                "of 'HANCESTRO:0304' are not allowed.",
            )
        ]

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
        ] = "HANCESTRO:0005,HANCESTRO:0014,HANCESTRO:0018"
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
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "self_reported_ethnicity_ontology_term_id"
        ]["error_message_suffix"]
        # Mouse organism ID
        validator.adata.obs.loc[validator.adata.obs.index[0], "organism_ontology_term_id"] = "NCBITaxon:10090"
        # Required to set to avoid development_stage_ontology_term_id errors
        validator.adata.obs.loc[validator.adata.obs.index[0], "development_stage_ontology_term_id"] = "MmusDv:0000003"
        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "HANCESTRO:0005"
        validator.validate_adata()
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'HANCESTRO:0005' in 'self_reported_ethnicity_ontology_term_id' is not a "
                "valid value of 'self_reported_ethnicity_ontology_term_id'.",
            )
        ]

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
        ] = "HANCESTRO:0014,HANCESTRO:0005"
        validator.validate_adata()
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'HANCESTRO:0014,HANCESTRO:0005' in 'self_reported_ethnicity_ontology_term_id' is not in "
                "ascending lexical order.",
            )
        ]

    def test_self_reported_ethnicity_ontology_term_id__invalid_delimiters(self, validator_with_adata):
        """
        Test error message for self_reported_ethnicity_ontology_term_id involving
        delimiters that are not specified in the schema definition yaml, such as whitespace
        """
        validator = validator_with_adata
        error_message_suffix = validator.schema_def["components"]["obs"]["columns"][
            "self_reported_ethnicity_ontology_term_id"
        ]["dependencies"][0]["error_message_suffix"]

        validator.adata.obs.loc[
            validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "HANCESTRO:0005, HANCESTRO:0014"
        validator.validate_adata()
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: ' HANCESTRO:0014' in 'self_reported_ethnicity_ontology_term_id' is not a valid ontology "
                "term id of 'HANCESTRO'.",
            )
        ]

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
        ] = "EFO:0000001,HANCESTRO:0004,HANCESTRO:0014,HANCESTRO:1"
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
        ] = "HANCESTRO:0014,HANCESTRO:0014"
        validator.validate_adata()
        assert validator.errors == [
            self.get_format_error_message(
                error_message_suffix,
                "ERROR: 'HANCESTRO:0014,HANCESTRO:0014' in 'self_reported_ethnicity_ontology_term_id' contains "
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
        ] = ["HANCESTRO:0005,HANCESTRO:0014"]
        validator.validate_adata()
        assert validator.errors[1] == self.get_format_error_message(
            error_message_suffix,
            "ERROR: '['HANCESTRO:0005,HANCESTRO:0014']' in 'self_reported_ethnicity_ontology_term_id' is not "
            "a valid ontology term value, it must be a string.",
        )

    def test_organism_ontology_term_id(self, validator_with_adata):
        """
        organism_ontology_term_id categorical with str categories. This MUST be a descendant of NCBITaxon:33208.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        # Setting "organism_ontology_term_id" to "EFO:0000001" is the fail case. However since this represents neither
        # human nor mouse, then two other columns that are dependent on it need to be set appropriately to avoid
        # other error messages: "development_stage_ontology_term_id" and "self_reported_ethnicity_ontology_term_id"
        obs.loc[obs.index[0], "organism_ontology_term_id"] = "EFO:0000001"
        obs.loc[obs.index[0], "development_stage_ontology_term_id"] = "unknown"
        obs.loc[
            obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "na"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'EFO:0000001' in 'organism_ontology_term_id' is not a valid "
            "ontology term id of 'NCBITaxon'. Only descendant term ids of 'NCBITaxon:33208' for metazoan are allowed."
        ]

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
        assert validator.errors == [
            "ERROR: 'EFO:0000001' in 'tissue_ontology_term_id' is not a valid ontology term id of "
            "'UBERON'. When 'tissue_type' is 'tissue' or 'organoid', 'tissue_ontology_term_id' MUST be a "
            "descendant term id of 'UBERON:0001062' (anatomical entity)."
        ]

    def test_tissue_ontology_term_id_cell_culture__suffix_in_term_id(self, validator_with_adata):
        """
        Cell Culture - Cannot include suffixes.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = "cell culture"
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = "CL:0000057 (cell culture)"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'CL:0000057 (cell culture)' in 'tissue_ontology_term_id' is not a valid ontology term id "
            "of 'CL'. When 'tissue_type' is 'cell culture', 'tissue_ontology_term_id' MUST be either a CL term "
            "(excluding 'CL:0000255' (eukaryotic cell), 'CL:0000257' (Eumycetozoan cell), "
            "and 'CL:0000548' (animal cell)) or 'unknown'."
        ]

    def test_tissue_ontology_term_id_cell_culture__not_a_CL_term(self, validator_with_adata):
        """
        Cell Culture - MUST be CL term
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_type"] = "cell culture"
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = "EFO:0000001"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'EFO:0000001' in 'tissue_ontology_term_id' is not a valid ontology term id of "
            "'CL'. When 'tissue_type' is 'cell culture', 'tissue_ontology_term_id' MUST be either a CL term "
            "(excluding 'CL:0000255' (eukaryotic cell), 'CL:0000257' (Eumycetozoan cell), "
            "and 'CL:0000548' (animal cell)) or 'unknown'."
        ]

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
        assert validator.errors == [
            f"ERROR: '{term}' in 'tissue_ontology_term_id' is not allowed. When 'tissue_type' is "
            f"'cell culture', 'tissue_ontology_term_id' MUST be either a CL term "
            "(excluding 'CL:0000255' (eukaryotic cell), 'CL:0000257' (Eumycetozoan cell), "
            "and 'CL:0000548' (animal cell)) or 'unknown'."
        ] or validator.errors == [
            f"ERROR: '{term}' in 'tissue_ontology_term_id' is a deprecated term id of 'CL'. When 'tissue_type' is "
            f"'cell culture', 'tissue_ontology_term_id' MUST be either a CL term "
            "(excluding 'CL:0000255' (eukaryotic cell), 'CL:0000257' (Eumycetozoan cell), "
            "and 'CL:0000548' (animal cell)) or 'unknown'."
        ]

    def test_tissue_ontology_term_id_organoid(self, validator_with_adata):
        """
        Organoid - must not accept suffixes like "(organoid)"
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = "UBERON:0000057 (organoid)"
        obs.tissue_type = obs.tissue_type.cat.add_categories(["organoid"])
        obs.loc[obs.index[0], "tissue_type"] = "organoid"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'UBERON:0000057 (organoid)' in 'tissue_ontology_term_id' is not a valid ontology term id of "
            "'UBERON'. When 'tissue_type' is 'tissue' or 'organoid', 'tissue_ontology_term_id' MUST be a "
            "descendant term id of 'UBERON:0001062' (anatomical entity)."
        ]

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
        assert validator.errors == [
            "ERROR: 'UBERON:0001062' in 'tissue_ontology_term_id' is not an allowed term id. "
            "When 'tissue_type' is 'tissue' or 'organoid', 'tissue_ontology_term_id' "
            "MUST be a descendant term id of 'UBERON:0001062' (anatomical entity)."
        ]

    def test_tissue_ontology_term_id_descendant_of_anatomical_entity__organoid(self, validator_with_adata):
        """
        Tissue ontology term ID must be a descendant term of 'UBERON:0001062' (anatomical entity) if tissue_type is
        organoid or tissue.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "tissue_ontology_term_id"] = "UBERON:0001062"
        obs.tissue_type = obs.tissue_type.cat.add_categories(["organoid"])
        obs.loc[obs.index[0], "tissue_type"] = "organoid"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: 'UBERON:0001062' in 'tissue_ontology_term_id' is not an allowed term id. "
            "When 'tissue_type' is 'tissue' or 'organoid', 'tissue_ontology_term_id' "
            "MUST be a descendant term id of 'UBERON:0001062' (anatomical entity)."
        ]

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
        obs.loc[obs.index[1], "suspension_type"] = invalid_suspension_type
        obs.loc[obs.index[1], "assay_ontology_term_id"] = assay
        validator.validate_adata()
        assert validator.errors == [
            f"ERROR: Column 'suspension_type' in dataframe 'obs' contains invalid values "
            f"'['{invalid_suspension_type}']'. Values must be one of {suspension_types} when "
            f"'assay_ontology_term_id' is {assay}"
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
        obs.loc[obs.index[1], "assay_ontology_term_id"] = assay
        obs.loc[obs.index[1], "suspension_type"] = invalid_suspension_type
        validator.validate_adata()
        assert validator.errors == [
            f"ERROR: Column 'suspension_type' in dataframe 'obs' contains invalid values "
            f"'['{invalid_suspension_type}']'. Values must be one of {suspension_types} when "
            f"'assay_ontology_term_id' is {assay} or its descendants"
        ]

    def test_suspension_type_with_descendant_term_id_failure(self, validator_with_adata):
        """
        suspension_id categorical with str categories. This field MUST be "cell", "nucleus", or "na". The allowed
        values depend on the assay_ontology_term_id. MUST support matching against ancestor term rules if specified.
        """
        validator = validator_with_adata
        obs = validator.adata.obs
        obs.loc[obs.index[0], "assay_ontology_term_id"] = "EFO:0022615"  # descendant of EFO:0008994
        obs.loc[obs.index[0], "suspension_type"] = "nucleus"

        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Column 'suspension_type' in dataframe 'obs' contains invalid values "
            "'['nucleus']'. Values must be one of ['na'] when "
            "'assay_ontology_term_id' is EFO:0008994 or its descendants"
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
        assert validator.warnings == [
            "WARNING: Data contains assay(s) that are not represented in the 'suspension_type' schema "
            "definition table. Ensure you have selected the most appropriate value for the assay(s) between "
            "'cell', 'nucleus', and 'na'. Please contact cellxgene@chanzuckerberg.com "
            "during submission so that the assay(s) can be added to the schema definition document."
        ]

    def test_categories_with_zero_values_warn(self, validator_with_adata):
        validator = validator_with_adata
        obs = validator.adata.obs
        modified_donor_id = obs["donor_id"].cat.add_categories("donor_3")
        obs["donor_id"] = modified_donor_id
        validator.validate_adata()
        assert validator.warnings == [
            "WARNING: Column 'donor_id' in dataframe 'obs' "
            "contains a category 'donor_3' with zero observations. "
            "These categories will be removed when `--add-labels` "
            "flag is present."
        ]

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

        with pytest.raises(ValueError):
            validator.validate_adata()

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

        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Some features are 'True' in 'feature_is_filtered' of dataframe 'var', "
            "but there are 1 non-zero values in the corresponding columns of the matrix 'X'. "
            "All values for these features must be 0."
        ]

        # Test that feature_is_filtered is a bool and not a string
        var["feature_is_filtered"] = "string"
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Column 'feature_is_filtered' in dataframe 'var' must be boolean, not 'object'."
        ]

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
        assert validator.errors == [f"ERROR: 'ENSG000' is not a valid feature ID in '{component_name}'."]

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
        assert validator.errors == [f"ERROR: 'ERCC-000000' is not a valid feature ID in '{component_name}'."]

    def test_should_warn_for_low_gene_count(self, validator_with_adata):
        """
        Raise a warning if there are too few genes
        """
        validator = validator_with_adata
        validator.schema_def["components"]["var"]["warn_if_less_than_rows"] = 100
        validator.validate_adata()
        assert validator.warnings == [
            "WARNING: Dataframe 'var' only has 4 rows. Features SHOULD NOT be filtered from expression matrix."
        ]

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

        with pytest.raises(ValueError):
            validator.validate_adata()

    def test_raw_var_column_name_uniqueness(self, validator_with_adata):
        validator = validator_with_adata
        original_var = validator.adata.var.copy()
        validator.adata.var = pd.concat([validator.adata.var, validator.adata.var["feature_is_filtered"]], axis=1)
        validator.adata.raw = validator.adata
        validator.adata.var = original_var  # Ensure only the raw.var has duplicate columns

        with pytest.raises(ValueError):
            validator.validate_adata()


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

    def test_uns_scipy_matrices_cannot_be_empty(self, validator_with_adata):
        validator = validator_with_adata

        validator.adata.uns["test"] = scipy.sparse.csr_matrix([[1]], dtype=int)
        validator.validate_adata()
        assert validator.errors == []

        validator.adata.uns["test"] = scipy.sparse.csr_matrix([[]], dtype=int)
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
            "ERROR: Annotated categorical field suspension_type must have at least 2 color options in uns[suspension_type_colors]. Found: []",
        ]

    def test_not_enough_color_options(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.uns["suspension_type_colors"] = numpy.array(["green"])
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: Annotated categorical field suspension_type must have at least 2 color options in uns[suspension_type_colors]. Found: ['green']"
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
            "ERROR: Colors in uns[suspension_type_colors] must be either all hex colors or all CSS4 named colors. Found: ['' 'green']"
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
        values in obsm cannot all be NaN
        """
        validator = validator_with_visium_assay
        obsm = validator.adata.obsm
        # It's okay if only one value is NaN
        obsm[key][0:100, 1] = numpy.nan
        validator.validate_adata()
        assert validator.errors == []

        # It's not okay if all values are NaN
        all_nan = numpy.full(obsm[key].shape, numpy.nan)
        obsm[key] = all_nan
        validator.validate_adata()
        assert validator.errors == [f"ERROR: adata.obsm['{key}'] contains all NaN values."]

    def test_obsm_values_no_X_embedding__non_spatial_dataset(self, validator_with_adata):
        validator = validator_with_adata
        validator.adata.obsm["harmony"] = validator.adata.obsm["X_umap"]
        validator.adata.uns["default_embedding"] = "harmony"
        del validator.adata.obsm["X_umap"]
        validator.validate_adata()
        assert validator.errors == [
            "ERROR: At least one embedding in 'obsm' has to have a key with an 'X_' prefix.",
        ]
        assert validator.is_spatial is False
        assert validator.warnings == [
            "WARNING: Dataframe 'var' only has 4 rows. Features SHOULD NOT be filtered from expression matrix.",
            "WARNING: Embedding key in 'adata.obsm' harmony is not 'spatial' nor does it start with 'X_'. "
            "Thus, it will not be available in Explorer",
            "WARNING: Validation of raw layer was not performed due to current errors, try again after fixing current errors.",
        ]

    def test_obsm_values_no_X_embedding__visium_dataset(self, validator_with_visium_assay):
        validator = validator_with_visium_assay
        validator.adata.uns["default_embedding"] = "spatial"
        del validator.adata.obsm["X_umap"]
        validator.validate_adata()
        assert validator.errors == []
        assert validator.is_spatial is True

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
        assert validator.warnings == [
            "WARNING: Dataframe 'var' only has 4 rows. Features SHOULD NOT be filtered from expression matrix.",
            "WARNING: Embedding key in 'adata.obsm' harmony is not 'spatial' nor does it start with 'X_'. "
            "Thus, it will not be available in Explorer",
            "WARNING: Validation of raw layer was not performed due to current errors, try again after fixing current errors.",
        ]
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
        assert validator.warnings == [
            "WARNING: Dataframe 'var' only has 4 rows. Features SHOULD NOT be filtered from expression matrix.",
            "WARNING: Embedding key in 'adata.obsm' 3D is not 'spatial' nor does it start with 'X_'. "
            "Thus, it will not be available in Explorer",
            "WARNING: Validation of raw layer was not performed due to current errors, try again after fixing current errors.",
        ]

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
        validator.adata = save_and_read_adata(adata)
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
        validator.adata = save_and_read_adata(adata)
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
        adata.varm["badsize"] = numpy.empty((4, 0))
        validator.adata = save_and_read_adata(adata)
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
        adata.varp["badsize"] = numpy.empty((4, 4, 0))
        validator.adata = save_and_read_adata(adata)
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
            "organism",
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
            - organism. categorical with str categories. This MUST be the human-readable name assigned
            to the value of organism_ontology_term_id.
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

    def test_obs_added_tissue_type_label__unknown(self, validator_with_adata):
        obs = validator_with_adata.adata.obs

        # Arrange
        obs.at["Y", "tissue_type"] = "cell culture"  # Already set in example data, just setting explicitly here
        obs.at["Y", "tissue_ontology_term_id"] = "unknown"  # Testing this term case
        validator_with_adata.validate_adata()  # Validate
        labeler = AnnDataLabelAppender(validator_with_adata)
        labeler._add_labels()  # Annotate

        assert labeler.adata.obs.at["Y", "tissue"] == "unknown"

    def test_obs_added_cell_type_label__unknown(self, validator_with_adata):
        obs = validator_with_adata.adata.obs

        # Arrange
        obs.at["Y", "cell_type_ontology_term_id"] = "unknown"  # Testing this term case
        validator_with_adata.validate_adata()  # Validate
        labeler = AnnDataLabelAppender(validator_with_adata)
        labeler._add_labels()  # Annotate

        assert labeler.adata.obs.at["Y", "cell_type"] == "unknown"

    def test_remove_unused_categories(self, label_writer, adata_with_labels):
        modified_donor_id = label_writer.adata.obs["donor_id"].cat.add_categories("donor_3")
        label_writer.adata.obs["donor_id"] = modified_donor_id
        case = unittest.TestCase()
        case.assertCountEqual(label_writer.adata.obs["donor_id"].dtype.categories, ["donor_1", "donor_2", "donor_3"])
        label_writer._remove_categories_with_zero_values()
        case.assertCountEqual(label_writer.adata.obs["donor_id"].dtype.categories, ["donor_1", "donor_2"])
