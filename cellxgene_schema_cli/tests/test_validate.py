import copy
import hashlib
import os
import tempfile
import unittest
from unittest import mock

import anndata
import numpy as np
from cellxgene_schema.ontology import OntologyChecker
from cellxgene_schema.schema import get_schema_definition
from cellxgene_schema.validate import Validator, validate
from cellxgene_schema.write_labels import AnnDataLabelAppender
from fixtures.examples_validate import (
    adata,
    adata_minimal,
    adata_with_labels,
    good_obs,
    good_obsm,
    good_uns,
    good_var,
    h5ad_invalid,
    h5ad_valid,
)
from scipy import sparse

# Tests for internal functions of the Validator and LabelWriter classes.


class TestFieldValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema_def = get_schema_definition()
        cls.OntologyChecker = OntologyChecker()

    def setUp(self):
        self.validator = Validator()
        self.validator.adata = adata_minimal
        self.column_name = "cell_type_ontology_term_id"
        self.column_schema = self.schema_def["components"]["obs"]["columns"][self.column_name]
        self.curie_constraints = self.schema_def["components"]["obs"]["columns"][self.column_name]["curie_constraints"]

    def test_schema_definition(self):
        """
        Tests that the definition of schema is well-defined
        """

        self.assertIsInstance(self.schema_def["components"], dict)
        self.assertIsInstance(self.schema_def["components"]["obs"], dict)
        self.assertIsInstance(self.schema_def["components"]["obs"]["columns"], dict)

        # Check that any columns in obs that are "curie" have "curie_constraints" and "ontologies" under the constraints
        for i in self.schema_def["components"]["obs"]["columns"]:
            self.assertTrue("type" in self.schema_def["components"]["obs"]["columns"][i])
            if i == "curie":
                self.assertIsInstance(
                    self.schema_def["components"]["obs"]["columns"][i]["curie_constraints"],
                    dict,
                )
                self.assertIsInstance(
                    self.schema_def["components"]["obs"]["columns"][i]["curie_constraints"]["ontolgies"],
                    list,
                )

                # Check that the allowed ontologies are in the ontology checker
                for ontology_name in self.schema_def["components"]["obs"]["columns"][i]["curie_constraints"][
                    "ontolgies"
                ]:
                    self.assertTrue(self.OntologyChecker.is_valid_ontology(ontology_name))

    def test_validate_ontology_good(self):
        self.validator._validate_curie("CL:0000066", self.column_name, self.curie_constraints)
        self.validator._validate_curie("CL:0000192", self.column_name, self.curie_constraints)
        self.assertFalse(self.validator.errors)

    def test_validate_ontology_wrong_ontology(self):
        self.validator._validate_curie("EFO:0009899", self.column_name, self.curie_constraints)
        self.assertTrue(self.validator.errors)

    def test_validate_ontology_wrong_term(self):
        self.validator._validate_curie("NO_TERM2", self.column_name, self.curie_constraints)
        self.assertTrue(self.validator.errors)


class TestAddLabelFunctions(unittest.TestCase):
    def setUp(self):
        # Set up test data
        self.test_adata = adata.copy()
        self.test_adata_with_labels = adata_with_labels
        self.schema_def = get_schema_definition()

        validator = Validator()
        validator.adata = self.test_adata
        validator.validate_adata()
        self.writer = AnnDataLabelAppender(validator)

    def test_get_dictionary_mapping_feature_id(self):
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
        self.assertEqual(self.writer._get_mapping_dict_feature_id(ids), expected_dict)

        # Bad
        ids = ["NO_GENE"]
        expected_dict = dict(zip(ids, labels))
        with self.assertRaises(KeyError):
            self.writer._get_mapping_dict_feature_id(ids)

    def test_get_dictionary_mapping_feature_reference(self):
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
        self.assertEqual(self.writer._get_mapping_dict_feature_reference(ids), expected_dict)

        # Bad
        ids = ["NO_GENE"]
        expected_dict = dict(zip(ids, labels))
        with self.assertRaises(KeyError):
            self.writer._get_mapping_dict_feature_id(ids)

    def test_get_dictionary_mapping_curie(self):
        # Good
        ids = ["CL:0000066", "CL:0000192"]
        labels = ["epithelial cell", "smooth muscle cell"]
        curie_constraints = self.schema_def["components"]["obs"]["columns"]["cell_type_ontology_term_id"][
            "curie_constraints"
        ]
        expected_dict = dict(zip(ids, labels))
        self.assertEqual(self.writer._get_mapping_dict_curie(ids, curie_constraints), expected_dict)

        ids = ["EFO:0009899", "EFO:0009922"]
        labels = ["10x 3' v2", "10x 3' v3"]
        curie_constraints = self.schema_def["components"]["obs"]["columns"]["assay_ontology_term_id"][
            "curie_constraints"
        ]
        expected_dict = dict(zip(ids, labels))
        self.assertEqual(self.writer._get_mapping_dict_curie(ids, curie_constraints), expected_dict)

        ids = ["MONDO:0100096"]
        labels = ["COVID-19"]
        curie_constraints = self.schema_def["components"]["obs"]["columns"]["disease_ontology_term_id"][
            "curie_constraints"
        ]
        expected_dict = dict(zip(ids, labels))
        self.assertEqual(self.writer._get_mapping_dict_curie(ids, curie_constraints), expected_dict)

        # Bad
        curie_constraints = self.schema_def["components"]["obs"]["columns"]["cell_type_ontology_term_id"][
            "curie_constraints"
        ]

        ids = ["CL:0000066", "CL:0000192_FOO"]
        with self.assertRaises(ValueError):
            self.writer._get_mapping_dict_curie(ids, curie_constraints)

        ids = ["CL:0000066", "CL:0000192", "UBERON:0002048"]
        with self.assertRaises(ValueError):
            self.writer._get_mapping_dict_curie(ids, curie_constraints)

        ids = ["CL:NO_TERM"]
        with self.assertRaises(ValueError):
            self.writer._get_mapping_dict_curie(ids, curie_constraints)

    def test__write__Success(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            labels_path = "/".join([temp_dir, "labels.h5ad"])
            self.writer.write_labels(labels_path)
        self.assertTrue(self.writer.was_writing_successful)
        self.assertFalse(self.writer.errors)

    def test__write__Fail(self):
        self.writer.adata.write_h5ad = mock.Mock(side_effect=Exception("Test Fail"))
        with tempfile.TemporaryDirectory() as temp_dir:
            labels_path = "/".join([temp_dir, "labels.h5ad"])
            self.writer.write_labels(labels_path)
        self.assertFalse(self.writer.was_writing_successful)
        self.assertTrue(self.writer.errors)


class TestIgnoreLabelFunctions(unittest.TestCase):
    def setUp(self):
        # Set up test data
        self.test_adata = adata
        self.test_adata_with_labels = adata_with_labels

    def test_validating_labeled_h5ad_should_fail_if_no_flag_set(self):
        validator = Validator()
        validator.adata = self.test_adata_with_labels
        is_valid = validator.validate_adata()

        self.assertFalse(is_valid)

    def test_validating_labeled_h5ad_should_pass_if_flag_set(self):
        validator = Validator(ignore_labels=True)
        validator.adata = copy.deepcopy(self.test_adata_with_labels)
        is_valid = validator.validate_adata()

        self.assertTrue(is_valid)


class TestValidate(unittest.TestCase):
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
            self.assertTrue(adata.X.has_canonical_format)
            self.assertTrue(adata.raw.X.has_canonical_format)
            self.assertTrue(success)
            self.assertListEqual(errors, [])
            self.assertTrue(is_seurat_convertible)
            self.assertTrue(os.path.exists(labels_path))
            expected_hash = "55fbc095218a01cad33390f534d6690af0ecd6593f27d7cd4d26e91072ea8835"
            original_hash = self.hash_file(h5ad_valid)
            self.assertNotEqual(
                original_hash,
                expected_hash,
                "Writing labels did not change the dataset from the original.",
            )

    def test__validate_with_h5ad_valid_and_without_labels(self):
        success, errors, is_seurat_convertible = validate(h5ad_valid)

        self.assertTrue(success)
        self.assertListEqual(errors, [])
        self.assertTrue(is_seurat_convertible)

    def test__validate_with_h5ad_invalid_and_with_labels(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            labels_path = "/".join([temp_dir, "labels.h5ad"])

            success, errors, is_seurat_convertible = validate(h5ad_invalid, labels_path)

            self.assertFalse(success)
            self.assertTrue(errors)
            self.assertTrue(is_seurat_convertible)
            self.assertFalse(os.path.exists(labels_path))

    def test__validate_with_h5ad_invalid_and_without_labels(self):
        success, errors, is_seurat_convertible = validate(h5ad_invalid)

        self.assertFalse(success)
        self.assertTrue(errors)
        self.assertTrue(is_seurat_convertible)


class TestSeuratConvertibility(unittest.TestCase):
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
        self.assertTrue(len(self.validator.warnings) == 1)
        self.assertFalse(self.validator.is_seurat_convertible)

        # Reducing nonzero count by 1, to within limit, makes it Seurat-convertible
        sparse_matrix_with_zero = sparse.csr_matrix(np.ones((good_obs.shape[0], good_var.shape[0]), dtype=np.float32))
        sparse_matrix_with_zero[0, 0] = 0
        self.validation_helper(sparse_matrix_with_zero)
        self.validator._validate_seurat_convertibility()
        self.assertTrue(len(self.validator.warnings) == 0)
        self.assertTrue(self.validator.is_seurat_convertible)

        # Dense matrices with a dimension that exceeds limit will fail -- zeros are irrelevant
        dense_matrix_with_zero = np.zeros((good_obs.shape[0], good_var.shape[0]), dtype=np.float32)
        self.validation_helper(dense_matrix_with_zero)
        self.validator.schema_def["max_size_for_seurat"] = 2**2 - 1
        self.validator._validate_seurat_convertibility()
        self.assertTrue(len(self.validator.warnings) == 1)
        self.assertFalse(self.validator.is_seurat_convertible)

        # Dense matrices with dimensions in bounds but total count over will succeed
        dense_matrix = np.ones((good_obs.shape[0], good_var.shape[0]), dtype=np.float32)
        self.validation_helper(dense_matrix)
        self.validator.schema_def["max_size_for_seurat"] = 2**3 - 1
        self.validator._validate_seurat_convertibility()
        self.assertTrue(len(self.validator.warnings) == 0)
        self.assertTrue(self.validator.is_seurat_convertible)

        # h5ad where raw matrix variable count != length of raw var variables array is not Seurat-convertible
        matrix = sparse.csr_matrix(np.zeros([good_obs.shape[0], good_var.shape[0]], dtype=np.float32))
        raw = anndata.AnnData(X=matrix, var=good_var)
        raw.var.drop("ENSSASG00005000004", axis=0, inplace=True)
        self.validation_helper(matrix, raw)
        self.validator._validate_seurat_convertibility()
        self.assertTrue(len(self.validator.warnings) == 1)
        self.assertFalse(self.validator.is_seurat_convertible)
