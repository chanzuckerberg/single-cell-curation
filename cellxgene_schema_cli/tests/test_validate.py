import hashlib
import os
import tempfile
from typing import Union
from unittest import mock

import anndata
import numpy as np
import pandas as pd
import pytest
from cellxgene_schema.ontology import OntologyChecker
from cellxgene_schema.schema import get_schema_definition
from cellxgene_schema.validate import Validator, validate
from cellxgene_schema.write_labels import AnnDataLabelAppender
from fixtures.examples_validate import (
    adata as adata_valid,
)
from fixtures.examples_validate import (
    adata_minimal,
    adata_with_labels,
    good_obs,
    good_obsm,
    good_uns,
    good_var,
    h5ad_invalid,
    h5ad_valid,
)
from numpy import ndarray
from scipy import sparse
from scipy.sparse import spmatrix

# Tests for internal functions of the Validator and LabelWriter classes.


@pytest.fixture(scope="class")
def schema_def():
    return get_schema_definition()


@pytest.fixture()
def validator_with_minimal_adata():
    validator = Validator()
    validator.adata = adata_minimal.copy()
    return validator


@pytest.fixture
def label_writer():
    validator = Validator()
    validator.adata = adata_valid.copy()
    validator.validate_adata()
    return AnnDataLabelAppender(validator)


@pytest.fixture
def valid_adata():
    return adata_valid.copy()


class TestFieldValidation:
    def test_schema_definition(self, schema_def):
        """
        Tests that the definition of schema is well-defined
        """
        ontology_checker = OntologyChecker()
        assert isinstance(schema_def["components"], dict)
        assert isinstance(schema_def["components"]["obs"], dict)
        assert isinstance(schema_def["components"]["obs"]["columns"], dict)

        # Check that any columns in obs that are "curie" have "curie_constraints" and "ontologies" under the constraints
        for i in schema_def["components"]["obs"]["columns"]:
            assert "type" in schema_def["components"]["obs"]["columns"][i]
            if schema_def["components"]["obs"]["columns"][i]["type"] == "curie":
                if "curie_constraints" in schema_def["components"]["obs"]["columns"][i]:
                    assert isinstance(
                        schema_def["components"]["obs"]["columns"][i]["curie_constraints"],
                        dict,
                    )
                    assert isinstance(
                        schema_def["components"]["obs"]["columns"][i]["curie_constraints"]["ontologies"],
                        list,
                    )

                    # Check that the allowed ontologies are in the ontology checker or 'NA' (special case)
                    for ontology_name in schema_def["components"]["obs"]["columns"][i]["curie_constraints"][
                        "ontologies"
                    ]:
                        if ontology_name != "NA":
                            assert ontology_checker.is_valid_ontology(ontology_name)
                else:
                    # if no curie_constraints in top-level for type curie, assert that 'dependencies' list exists
                    assert isinstance(
                        schema_def["components"]["obs"]["columns"][i]["dependencies"],
                        list,
                    )

    def test_validate_ontology_good(self, validator_with_minimal_adata, schema_def):
        validator = validator_with_minimal_adata
        column_name = "cell_type_ontology_term_id"
        curie_constraints = schema_def["components"]["obs"]["columns"][column_name]["curie_constraints"]
        validator._validate_curie("CL:0000066", column_name, curie_constraints)
        validator._validate_curie("CL:0000192", column_name, curie_constraints)
        assert not validator.errors

    def test_validate_ontology_wrong_ontology(self, validator_with_minimal_adata, schema_def):
        validator = validator_with_minimal_adata
        column_name = "cell_type_ontology_term_id"
        curie_constraints = schema_def["components"]["obs"]["columns"][column_name]["curie_constraints"]
        validator._validate_curie("EFO:0009899", column_name, curie_constraints)
        assert validator.errors

    def test_validate_ontology_wrong_term(self, validator_with_minimal_adata, schema_def):
        validator = validator_with_minimal_adata
        column_name = "cell_type_ontology_term_id"
        curie_constraints = schema_def["components"]["obs"]["columns"][column_name]["curie_constraints"]
        validator._validate_curie("NO_TERM2", column_name, curie_constraints)
        assert validator.errors


class TestAddLabelFunctions:
    def test_get_dictionary_mapping_feature_id(self, label_writer):
        # Good
        ids = [
            "ERCC-00002",
            "ENSG00000127603",
            "ENSMUSG00000059552",
            "ENSSASG00005000004",
        ]
        labels = [
            "ERCC-00002 (spike-in control)",
            "MACF1",
            "Trp53",
            "S",
        ]
        expected_dict = dict(zip(ids, labels))
        assert label_writer._get_mapping_dict_feature_id(ids), expected_dict

        # Bad
        ids = ["NO_GENE"]
        with pytest.raises(KeyError):
            label_writer._get_mapping_dict_feature_id(ids)

    def test_get_dictionary_mapping_feature_reference(self, label_writer):
        # Good
        ids = [
            "ERCC-00002",
            "ENSG00000127603",
            "ENSMUSG00000059552",
            "ENSSASG00005000004",
        ]
        labels = [
            "NCBITaxon:32630",
            "NCBITaxon:9606",
            "NCBITaxon:10090",
            "NCBITaxon:2697049",
        ]
        expected_dict = dict(zip(ids, labels))
        assert label_writer._get_mapping_dict_feature_reference(ids) == expected_dict

        # Bad
        ids = ["NO_GENE"]
        with pytest.raises(KeyError):
            label_writer._get_mapping_dict_feature_id(ids)

    def test_get_dictionary_mapping_feature_length(self, label_writer):
        # Good
        ids = [
            "ERCC-00002",
            "ENSG00000127603",
            "ENSMUSG00000059552",
            "ENSSASG00005000004",
        ]
        # values derived from csv
        gene_lengths = [
            0,  # non-gene feature, so set to 0 regardless of csv value
            42738,
            4045,
            3822,
        ]
        expected_dict = dict(zip(ids, gene_lengths))
        assert label_writer._get_mapping_dict_feature_length(ids) == expected_dict

        # Bad
        ids = ["NO_GENE"]
        with pytest.raises(KeyError):
            label_writer._get_mapping_dict_feature_id(ids)

    @pytest.mark.parametrize(
        "ids,labels,curie_constraints",
        [
            (["CL:0000066", "CL:0000192"], ["epithelial cell", "smooth muscle cell"], "cell_type_ontology_term_id"),
            (["EFO:0009899", "EFO:0009922"], ["10x 3' v2", "10x 3' v3"], "assay_ontology_term_id"),
            (["MONDO:0100096"], ["COVID-19"], "disease_ontology_term_id"),
        ],
    )
    def test_get_dictionary_mapping_curie__good(self, schema_def, label_writer, ids, labels, curie_constraints):
        # Good
        curie_constraints = schema_def["components"]["obs"]["columns"][curie_constraints]["curie_constraints"]
        expected_dict = dict(zip(ids, labels))
        assert label_writer._get_mapping_dict_curie(ids, curie_constraints) == expected_dict

    def test_get_dictionary_mapping_curie__self_reported_ethnicity_ontology_term_id(self, schema_def, label_writer):
        ids = ["HANCESTRO:0005", "HANCESTRO:0014", "HANCESTRO:0005,HANCESTRO:0014", "unknown"]
        labels = ["European", "Hispanic or Latin American", "European,Hispanic or Latin American", "unknown"]
        curie_constraints = schema_def["components"]["obs"]["columns"]["self_reported_ethnicity_ontology_term_id"][
            "dependencies"
        ][0]["curie_constraints"]
        expected_dict = dict(zip(ids, labels))
        assert label_writer._get_mapping_dict_curie(ids, curie_constraints) == expected_dict

    @pytest.mark.parametrize(
        "ids", [["CL:0000066", "CL:0000192_FOO"], ["CL:0000066", "CL:0000192", "UBERON:0002048"], ["CL:NO_TERM"]]
    )
    def test_get_dictionary_mapping_curie__bad(self, schema_def, label_writer, ids):
        curie_constraints = schema_def["components"]["obs"]["columns"]["cell_type_ontology_term_id"][
            "curie_constraints"
        ]
        with pytest.raises(ValueError):
            label_writer._get_mapping_dict_curie(ids, curie_constraints)

    def test__write__Success(self, label_writer):
        with tempfile.TemporaryDirectory() as temp_dir:
            labels_path = "/".join([temp_dir, "labels.h5ad"])
            label_writer.write_labels(labels_path)
        assert label_writer.was_writing_successful
        assert not label_writer.errors

    def test__write__Fail(self, label_writer):
        label_writer.adata.write_h5ad = mock.Mock(side_effect=Exception("Test Fail"))
        with tempfile.TemporaryDirectory() as temp_dir:
            labels_path = "/".join([temp_dir, "labels.h5ad"])
            label_writer.write_labels(labels_path)
        assert not label_writer.was_writing_successful
        assert label_writer.errors


class TestIgnoreLabelFunctions:
    def test_validating_labeled_h5ad_should_fail_if_no_flag_set(self):
        validator = Validator()
        validator.adata = adata_with_labels.copy()
        is_valid = validator.validate_adata()

        assert not is_valid

    def test_validating_labeled_h5ad_should_pass_if_flag_set(self):
        validator = Validator(ignore_labels=True)
        validator.adata = adata_with_labels.copy()
        is_valid = validator.validate_adata()

        assert is_valid


class TestValidate:
    @staticmethod
    def hash_file(file_name: str) -> str:
        with open(file_name, "rb") as f:
            # Read the contents of the file in chunks
            chunk_size = 1024
            hasher = hashlib.sha256()
            while chunk := f.read(chunk_size):
                hasher.update(chunk)
        return hasher.hexdigest()

    def test__validate_with_h5ad_valid_and_labels(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            labels_path = "/".join([temp_dir, "labels.h5ad"])

            success, errors, is_seurat_convertible = validate(h5ad_valid, labels_path)

            import anndata as ad

            adata = ad.read_h5ad(labels_path)
            assert adata.X.has_canonical_format
            assert adata.raw.X.has_canonical_format
            assert success
            assert not errors
            assert is_seurat_convertible
            assert os.path.exists(labels_path)
            expected_hash = "55fbc095218a01cad33390f534d6690af0ecd6593f27d7cd4d26e91072ea8835"
            original_hash = self.hash_file(h5ad_valid)
            assert original_hash != expected_hash, "Writing labels did not change the dataset from the original."

    def test__validate_with_h5ad_valid_and_without_labels(self):
        success, errors, is_seurat_convertible = validate(h5ad_valid)

        assert success
        assert not errors
        assert is_seurat_convertible

    def test__validate_with_h5ad_invalid_and_with_labels(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            labels_path = "/".join([temp_dir, "labels.h5ad"])

            success, errors, is_seurat_convertible = validate(h5ad_invalid, labels_path)

            assert not success
            assert errors
            assert is_seurat_convertible
            assert not os.path.exists(labels_path)

    def test__validate_with_h5ad_invalid_and_without_labels(self):
        success, errors, is_seurat_convertible = validate(h5ad_invalid)

        assert not success
        assert errors
        assert is_seurat_convertible


class TestSeuratConvertibility:
    def validation_helper(self, matrix, raw=None):
        data = anndata.AnnData(X=matrix, obs=good_obs, uns=good_uns, obsm=good_obsm, var=good_var)
        if raw:
            data.raw = raw
        self.validator: Validator = Validator()
        self.validator._set_schema_def()
        self.validator.schema_def["max_size_for_seurat"] = 2**3 - 1  # Reduce size required to fail (faster tests)
        self.validator.adata = data

    def test_determine_seurat_convertibility(self):
        # Sparse matrix with too many nonzero values is not Seurat-convertible
        sparse_matrix_too_large = sparse.csr_matrix(np.ones((good_obs.shape[0], good_var.shape[0]), dtype=np.float32))
        self.validation_helper(sparse_matrix_too_large)
        self.validator._validate_seurat_convertibility()
        assert len(self.validator.warnings) == 1
        assert not self.validator.is_seurat_convertible

        # Reducing nonzero count by 1, to within limit, makes it Seurat-convertible
        sparse_matrix_with_zero = sparse.csr_matrix(np.ones((good_obs.shape[0], good_var.shape[0]), dtype=np.float32))
        sparse_matrix_with_zero[0, 0] = 0
        self.validation_helper(sparse_matrix_with_zero)
        self.validator._validate_seurat_convertibility()
        assert len(self.validator.warnings) == 0
        assert self.validator.is_seurat_convertible

        # Dense matrices with a dimension that exceeds limit will fail -- zeros are irrelevant
        dense_matrix_with_zero = np.zeros((good_obs.shape[0], good_var.shape[0]), dtype=np.float32)
        self.validation_helper(dense_matrix_with_zero)
        self.validator.schema_def["max_size_for_seurat"] = 2**2 - 1
        self.validator._validate_seurat_convertibility()
        assert len(self.validator.warnings) == 1
        assert not self.validator.is_seurat_convertible

        # Dense matrices with dimensions in bounds but total count over will succeed
        dense_matrix = np.ones((good_obs.shape[0], good_var.shape[0]), dtype=np.float32)
        self.validation_helper(dense_matrix)
        self.validator.schema_def["max_size_for_seurat"] = 2**3 - 1
        self.validator._validate_seurat_convertibility()
        assert len(self.validator.warnings) == 0
        assert self.validator.is_seurat_convertible

        # h5ad where raw matrix variable count != length of raw var variables array is not Seurat-convertible
        matrix = sparse.csr_matrix(np.zeros([good_obs.shape[0], good_var.shape[0]], dtype=np.float32))
        raw = anndata.AnnData(X=matrix, var=good_var)
        raw.var.drop("ENSSASG00005000004", axis=0, inplace=True)
        self.validation_helper(matrix, raw)
        self.validator._validate_seurat_convertibility()
        assert len(self.validator.errors) == 1
        assert not self.validator.is_seurat_convertible
        assert not self.validator.is_valid


class TestValidatorValidateDataFrame:
    @pytest.mark.parametrize("_type", [np.int64, np.int32, int, np.float64, np.float32, float, str])
    def test_succeed_categorical_types(self, tmp_path, _type, valid_adata):
        # Arrange
        categories = [*map(_type, range(adata_valid.n_obs))]
        self._add_catagorical_obs(valid_adata, categories)
        validator = self._create_validator(valid_adata)

        # Act
        validator._validate_dataframe("obs")

        # Assert
        assert not validator.errors
        valid_adata.write_h5ad(f"{tmp_path}/test.h5ad")  # Succeed write

    def test_fail_categorical_mixed_types(self, tmp_path, valid_adata):
        # Arrange
        categories = ["hello", 123]
        self._add_catagorical_obs(valid_adata, categories)
        validator = self._create_validator(valid_adata)

        # Act
        validator._validate_dataframe("obs")

        # Assert
        assert "in dataframe 'obs' containes 2 categorical types. Only one type is allowed." in validator.errors[0]
        self._fail_write_h5ad(tmp_path, valid_adata)

    def test_fail_categorical_bool(self, tmp_path, valid_adata):
        # Arrange
        categories = [True, False]
        self._add_catagorical_obs(valid_adata, categories)
        validator = self._create_validator(valid_adata)

        # Act
        validator._validate_dataframe("obs")

        # Assert
        assert "in dataframe 'obs' containes illegal_categorical_types={<class 'bool'>}." in validator.errors[0]
        self._fail_write_h5ad(tmp_path, valid_adata)

    def _add_catagorical_obs(self, adata, categories):
        t = pd.CategoricalDtype(categories=categories)
        adata.obs["test_cat"] = pd.Series(data=categories, index=["X", "Y"], dtype=t)

    def _create_validator(self, adata):
        validator = Validator()
        validator._set_schema_def()
        validator.adata = adata
        return validator

    def _fail_write_h5ad(self, tmp_path, adata):
        with pytest.raises(TypeError):
            # If this tests starts to fail here it means the anndata version has be upgraded and this check is no
            # longer needed
            adata.write_h5ad(f"{tmp_path}/test.h5ad")

    def test_fail_mixed_column_types(self, tmp_path, valid_adata):
        # Arrange
        valid_adata.obs["mixed"] = pd.Series(data=["1234", 0], index=["X", "Y"])
        validator = self._create_validator(valid_adata)

        # Act
        validator._validate_dataframe("obs")

        # Assert
        assert "in dataframe 'obs' cannot contain mixed types." in validator.errors[0]
        self._fail_write_h5ad(tmp_path, valid_adata)


class TestIsRaw:
    @staticmethod
    def create_validator(data: Union[ndarray, spmatrix], matrix_format: str) -> Validator:
        """
        Create a sample AnnData instance with the given data and format.

        :param data: The data matrix.
        :param matrix_format: The format of the data matrix (e.g., "dense", "csr", "csc").

        :return anndata.AnnData: An AnnData instance with the specified data and format.
        """
        validator = Validator()

        adata = anndata.AnnData(X=data)
        adata.obsm["X_" + matrix_format] = data

        validator.adata = adata
        return validator

    @pytest.mark.parametrize(
        "data, matrix_format, expected_result",
        [
            # Test case with integer values in a dense matrix
            (np.array([[1, 2, 3], [4, 5, 6]], dtype=int), "dense", True),
            # Test case with float values in a dense matrix
            (np.array([[1.1, 2.2, 3.3], [4.4, 5.5, 6.6]]), "dense", False),
            # Test case with integer values in a sparse matrix (CSR format)
            (sparse.csr_matrix([[1, 0, 3], [0, 5, 0]], dtype=int), "csr", True),
            # Test case with float values in a sparse matrix (CSC format)
            (sparse.csc_matrix([[1.1, 0, 3.3], [0, 5.5, 0]]), "csc", False),
            # Test case with mixed integer and float values in a dense matrix
            (np.array([[1, 2.2, 3], [4.4, 5, 6.6]]), "dense", False),
        ],
    )
    def test_has_valid_raw(self, data, matrix_format, expected_result):
        validator = self.create_validator(data, matrix_format)
        assert validator._has_valid_raw() == expected_result

    @mock.patch("cellxgene_schema.validate.get_matrix_format", return_value="unknown")
    def test_has_valid_raw_with_unknown_format(self, mock_get_matrix_format):
        data = np.array([[1, 2, 3], [4, 5, 6]], dtype=int)
        validator = self.create_validator(data, "unknown")
        with pytest.raises(AssertionError):
            validator._has_valid_raw()
