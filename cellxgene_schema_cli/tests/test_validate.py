import hashlib
import os
import re
import tempfile
from typing import Union
from unittest import mock

import anndata
import numpy as np
import pandas as pd
import pytest
from cellxgene_ontology_guide.entities import Ontology
from cellxgene_schema.schema import get_schema_definition
from cellxgene_schema.validate import (
    ERROR_SUFFIX_SPARSE_FORMAT,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_FORBIDDEN,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_IN_TISSUE_0,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_NOTNULL,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_REQUIRED,
    SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE,
    SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE_VISIUM_11MM,
    Validator,
    validate,
)
from cellxgene_schema.write_labels import AnnDataLabelAppender
from dask.array import from_array
from fixtures.examples_validate import adata as adata_valid
from fixtures.examples_validate import (
    adata_minimal,
    adata_slide_seqv2,
    adata_spatial_is_single_false,
    adata_visium,
    adata_with_labels,
    good_obs,
    good_uns,
    good_uns_with_visium_spatial,
    h5ad_invalid,
    h5ad_valid,
    visium_library_id,
)
from numpy import ndarray, zeros
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
def label_writer(valid_adata):
    return AnnDataLabelAppender(adata_valid)


@pytest.fixture
def valid_adata():
    return adata_valid.copy()


class TestFieldValidation:
    def test_schema_definition(self, schema_def):
        """
        Tests that the definition of schema is well-defined
        """
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
                            assert getattr(Ontology, ontology_name, None) is not None
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
            "FBtr0472816_df_nrg",
            "ENSDARG00000009657",
            "WBGene00000003",
        ]
        labels = ["ERCC-00002 (spike-in control)", "MACF1", "Trp53", "S", "FBtr0472816_df_nrg", "fgfr1op2", "aat-2"]
        expected_dict = dict(zip(ids, labels))
        assert label_writer._get_mapping_dict_feature_id(ids), expected_dict

        # Bad
        ids = ["NO_GENE"]
        with pytest.raises(ValueError):
            label_writer._get_mapping_dict_feature_id(ids)

    def test_get_dictionary_mapping_feature_reference(self, label_writer):
        # Good
        ids = [
            "ERCC-00002",
            "ENSG00000127603",
            "ENSMUSG00000059552",
            "ENSSASG00005000004",
            "FBtr0472816_df_nrg",
            "ENSDARG00000009657",
            "WBGene00000003",
        ]
        labels = [
            "NCBITaxon:32630",
            "NCBITaxon:9606",
            "NCBITaxon:10090",
            "NCBITaxon:2697049",
            "NCBITaxon:7227",
            "NCBITaxon:7955",
            "NCBITaxon:6239",
        ]
        expected_dict = dict(zip(ids, labels))
        assert label_writer._get_mapping_dict_feature_reference(ids) == expected_dict

        # Bad
        ids = ["NO_GENE"]
        with pytest.raises(ValueError):
            label_writer._get_mapping_dict_feature_id(ids)

    def test_get_dictionary_mapping_feature_length(self, label_writer):
        # Good
        ids = [
            "ERCC-00002",
            "ENSG00000127603",
            "ENSMUSG00000059552",
            "ENSSASG00005000004",
            "FBtr0472816_df_nrg",
            "ENSDARG00000009657",
            "WBGene00000003",
        ]
        # values derived from csv
        gene_lengths = [1061, 2821, 1797, 3822, 22, 1088, 1738]
        expected_dict = dict(zip(ids, gene_lengths))
        assert label_writer._get_mapping_dict_feature_length(ids) == expected_dict

        # Bad
        ids = ["NO_GENE"]
        with pytest.raises(ValueError):
            label_writer._get_mapping_dict_feature_id(ids)

    def test_get_dictionary_mapping_feature_type(self, label_writer):
        # Good
        ids = [
            "ERCC-00002",
            "ENSG00000127603",
            "ENSMUSG00000059552",
            "ENSSASG00005000004",
            "FBtr0472816_df_nrg",
            "ENSDARG00000009657",
            "WBGene00000003",
        ]
        # values derived from csv
        gene_types = [
            "synthetic",
            "protein_coding",
            "protein_coding",
            "protein_coding",
            "ncRNA",
            "protein_coding",
            "protein_coding",
        ]
        expected_dict = dict(zip(ids, gene_types))
        assert label_writer._get_mapping_dict_feature_type(ids) == expected_dict

        # Bad
        ids = ["NO_GENE_BAD"]
        with pytest.raises(ValueError):
            label_writer._get_mapping_dict_feature_type(ids)

    def test_get_dictionary_mapping_feature_biotype(self, label_writer):
        # Good
        ids = [
            "ERCC-00002",
            "ENSG00000127603",
            "ENSMUSG00000059552",
            "ENSSASG00005000004",
            "FBtr0472816_df_nrg",
            "ENSDARG00000009657",
            "WBGene00000003",
        ]
        # Values derived from csv
        biotypes = [
            "spike-in",
            "gene",
            "gene",
            "gene",
            "gene",
            "gene",
            "gene",
        ]
        expected_dict = dict(zip(ids, biotypes))
        assert label_writer._get_mapping_dict_feature_biotype(ids) == expected_dict

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
        ids = ["HANCESTRO:0005", "HANCESTRO:0014", "HANCESTRO:0005 || HANCESTRO:0014", "unknown"]
        labels = ["European", "Hispanic or Latin American", "European || Hispanic or Latin American", "unknown"]
        curie_constraints = schema_def["components"]["obs"]["columns"]["self_reported_ethnicity_ontology_term_id"][
            "dependencies"
        ][0]["curie_constraints"]
        expected_dict = dict(zip(ids, labels))
        assert label_writer._get_mapping_dict_curie(ids, curie_constraints) == expected_dict

    def test_get_dictionary_mapping_curie__disease_ontology_term_id(self, schema_def, label_writer):
        ids = ["PATO:0000461", "MONDO:1030008", "MONDO:0004604 || MONDO:0043004 || MONDO:0800349 || MONDO:1030008"]
        labels = [
            "normal",
            "mitral valve insufficiency",
            "Hodgkin's lymphoma, lymphocytic-histiocytic predominance || Weil's disease || atrial fibrillation, familial, 16 || mitral valve insufficiency",
        ]
        curie_constraints = schema_def["components"]["obs"]["columns"]["disease_ontology_term_id"]["curie_constraints"]
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
            assert label_writer.write_labels(labels_path)
        assert not label_writer.errors

    def test__write__Fail(self, label_writer):
        label_writer.adata.write_h5ad = mock.Mock(side_effect=Exception("Test Fail"))
        with tempfile.TemporaryDirectory() as temp_dir:
            labels_path = "/".join([temp_dir, "labels.h5ad"])
            assert not label_writer.write_labels(labels_path)
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

            success, errors, _ = validate(h5ad_valid, labels_path)

            import anndata as ad

            adata = ad.read_h5ad(labels_path)
            assert adata.X.has_canonical_format
            assert adata.raw.X.has_canonical_format
            assert success
            assert not errors
            assert os.path.exists(labels_path)
            expected_hash = "55fbc095218a01cad33390f534d6690af0ecd6593f27d7cd4d26e91072ea8835"
            original_hash = self.hash_file(h5ad_valid)
            assert original_hash != expected_hash, "Writing labels did not change the dataset from the original."

    def test__validate_with_h5ad_valid_and_without_labels(self):
        success, errors, _ = validate(h5ad_valid)

        assert success
        assert not errors

    def test__validate_with_h5ad_invalid_and_with_labels(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            labels_path = "/".join([temp_dir, "labels.h5ad"])

            success, errors, _ = validate(h5ad_invalid, labels_path)

            assert not success
            assert errors
            assert not os.path.exists(labels_path)

    def test__validate_with_h5ad_invalid_and_without_labels(self):
        success, errors, _ = validate(h5ad_invalid)

        assert not success
        assert errors


class TestCheckSpatial:
    @pytest.mark.parametrize(
        "assay_ontology_term_id, expected_is_visium",
        [
            # Parent term for Visium Spatial Gene Expression. This term and all its descendants are Visium
            ("EFO:0022858", True),
            # Visium Spatial Gene Expression V1
            ("EFO:0022857", True),
            # Visium CytAssist Spatial Gene Expression V2
            ("EFO:0022858", True),
            # Visium CytAssist Spatial Gene Expression, 11mm
            ("EFO:0022860", True),
            # Visium CytAssist Spatial Gene Expression, 6.5mm
            ("EFO:0022859", True),
            # Random other EFO term
            ("EFO:0003740", False),
        ],
    )
    def test__is_visium_descendant(self, assay_ontology_term_id, expected_is_visium):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id

        assert validator._is_visium_including_descendants() == expected_is_visium

    def test__validate_spatial_visium_ok(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator._visium_and_is_single_true_matrix_size = 2
        # Confirm spatial is valid.
        validator.validate_adata()
        assert not validator.errors

    def test__forbid_generic_visium(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator._visium_and_is_single_true_matrix_size = 2

        # set assay to the generic visium term
        validator.adata.obs["assay_ontology_term_id"] = "EFO:0010961"

        # Confirm this triggers FORBIDDIN ERROR and downstream errors due to invalid spatial term.
        validator.validate_adata()
        EXPECTED_FORBIDDEN_ERROR = "ERROR: Invalid spatial assay. obs['assay_ontology_term_id'] must be a descendant of EFO:0010961 but NOT EFO:0010961 itself. "
        assert len(validator.errors) == 5
        assert EXPECTED_FORBIDDEN_ERROR in validator.errors

    @mock.patch("cellxgene_schema.validate.VISIUM_AND_IS_SINGLE_TRUE_MATRIX_SIZE", 2)
    def test__validate_from_file(self):
        """Testing compatibility with SparseDataset types in Anndata"""
        validator: Validator = Validator()
        validator._set_schema_def()
        validator._visium_and_is_single_true_matrix_size = 2
        with tempfile.TemporaryDirectory() as temp_dir:
            file_path = os.path.join(temp_dir, "visium.h5ad")
            adata_visium.write_h5ad(file_path)
            # Confirm spatial is valid.
            validator.validate_adata(file_path)
        assert not validator.errors

    def test__validate_spatial_visium_dense_matrix_ok(self):
        """
        Test visium specific requirements on a dense X matrix
        """
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator._visium_and_is_single_true_matrix_size = 2
        _Xdense = validator.adata.X.compute()
        _Xdense[_Xdense == 0] = 1  # ensure the matrix doesn't trigger sparsity error
        validator.adata.X = from_array(_Xdense.toarray())  # daskify
        validator.adata.raw = validator.adata.copy()
        validator.adata.raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        # Confirm spatial is valid.
        validator.validate_adata()
        assert not validator.errors

    @pytest.mark.parametrize(
        "is_single, is_primary_data, validation_result",
        [
            (False, False, True),
            (False, True, False),
            # Same as above cases but with numpy bools
            (np.bool_(False), np.bool_(False), True),
            (np.bool_(False), np.bool_(True), False),
        ],
    )
    def test__validate_spatial_is_primary_data(self, is_single, is_primary_data, validation_result):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"] = {"is_single": is_single}
        del validator.adata.obsm["spatial"]
        validator.adata.obs.pop("array_col")
        validator.adata.obs.pop("array_row")
        validator.adata.obs.pop("in_tissue")
        validator.adata.obs["is_primary_data"] = is_primary_data
        assert validator.validate_adata() == validation_result

    def test__validate_spatial_slide_seqV2_ok(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_slide_seqv2.copy()

        # Confirm spatial is valid.
        validator.validate_adata()
        assert not validator.errors

    @pytest.mark.parametrize("spatial", [None, "invalid", 1, 1.0, True])
    def test__validate_spatial_type_error(self, spatial):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"] = spatial

        # Confirm key type dict is required.
        validator.validate_adata()
        assert (
            "ERROR: A dict in uns['spatial'] is required when obs['assay_ontology_term_id'] is either a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2)."
            in validator.errors
        )

    def test__validate_spatial_is_single_false_ok(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_spatial_is_single_false.copy()

        # Confirm spatial is valid.
        validator.validate_adata()
        assert not validator.errors

    def test__validate_spatial_forbidden_if_not_visium_or_slide_seqv2(self):
        validator: Validator = Validator()
        validator._set_schema_def()

        # Create Visium adata (with spatial) with 10x 3' v2 obs.
        validator.adata = adata_visium.copy()
        validator.adata.obs = good_obs.copy()

        # Confirm spatial is not allowed for 10x 3' v2.
        validator._check_spatial_uns()
        assert validator.errors == [
            "uns['spatial'] is only allowed when obs['assay_ontology_term_id'] is either "
            "a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2)"
        ]

    @pytest.mark.parametrize(
        "assay_ontology_term_id, is_descendant",
        [("EFO:0022859", True), ("EFO:0022858", True), ("EFO:0030029", False), ("EFO:0002697", False)],
    )
    def test__validate_spatial_required_if_visium(self, assay_ontology_term_id, is_descendant):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id

        if is_descendant:
            # check pass if 'spatial' included
            validator.adata.uns = good_uns_with_visium_spatial.copy()
            validator._check_spatial_uns()
            assert len(validator.errors) == 0
            validator.reset()

            # check fail if 'spatial' not included
            validator.adata.uns = good_uns.copy()
            validator._check_spatial_uns()
            assert validator.errors == [
                "A dict in uns['spatial'] is required when obs['assay_ontology_term_id'] is "
                "either a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2)."
            ]
            validator.reset()
        else:
            # check fail if 'spatial' included
            validator.adata.uns = good_uns_with_visium_spatial.copy()
            validator._check_spatial_uns()
            assert len(validator.errors) == 1
            validator.reset()

    def test__validate_spatial_required_if_slide_seqV2(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_slide_seqv2.copy()
        validator.adata.uns = good_uns.copy()

        # Confirm spatial is required for Slide-seqV2.
        validator._check_spatial_uns()
        assert validator.errors == [
            "A dict in uns['spatial'] is required when obs['assay_ontology_term_id'] is either a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2)."
        ]

    def test__validate_spatial_allowed_keys_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"]["invalid_key"] = True

        # Confirm additional key is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            "uns['spatial'] must contain only two top-level keys: 'is_single' and a library_id. "
            "More than two top-level keys detected:" in validator.errors[0]
        )

    @pytest.mark.parametrize(
        "assay_ontology_term_id, is_descendant",
        [("EFO:0022859", True), ("EFO:0022858", True), ("EFO:0030029", False), ("EFO:0002697", False)],
    )
    def test__validate_is_single_required_visium_error(self, assay_ontology_term_id, is_descendant):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        validator.adata.uns["spatial"].pop("is_single")
        validator._check_spatial_uns()

        if is_descendant:
            # if spatial, MUST specify `is_single`
            assert "uns['spatial'] must contain the key 'is_single'." in validator.errors[0]
        else:
            # if not spatial, MUST NOT speciffy `is_single`
            assert validator.errors == [
                "uns['spatial'] is only allowed when obs['assay_ontology_term_id'] is either a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2)"
            ]

    def test__validate_is_single_required_slide_seqV2_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()

        # Remove is_single key to trigger error, but add dummy library_id to pass earlier truthy check
        # on spatial existing.
        validator.adata = adata_slide_seqv2.copy()
        validator.adata.uns["spatial"] = {
            "library_id": "test_library_id",
        }

        # Confirm is_single is identified as required.
        validator._check_spatial_uns()
        assert validator.errors
        assert "uns['spatial'] must contain the key 'is_single'." in validator.errors[0]

    def test__validate_is_single_boolean_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"]["is_single"] = "invalid"

        # Confirm is_single is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors
        assert "uns['spatial']['is_single'] must be of boolean type" in validator.errors[0]

    def test__validate_library_id_forbidden_if_slide_seqV2(self):
        validator: Validator = Validator()
        validator._set_schema_def()

        # Add library_id to Slide-seqV2 uns to trigger error.
        validator.adata = adata_slide_seqv2.copy()
        validator.adata.uns = good_uns_with_visium_spatial

        # Confirm library_id is not allowed for Slide-seqV2.
        validator._check_spatial_uns()
        assert len(validator.errors) == 1
        assert f"uns['spatial'][library_id] {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_FORBIDDEN}." in validator.errors[0]

    def test__validate_library_id_forbidden_if_visium_or_is_single_false(self):
        validator: Validator = Validator()
        validator._set_schema_def()

        # Set is_single to False to trigger error.
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"]["is_single"] = np.bool_(False)

        # Confirm library_id is not allowed for Visium if is_single is False.
        validator._check_spatial_uns()
        assert len(validator.errors) == 1
        assert f"uns['spatial'][library_id] {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_FORBIDDEN}." in validator.errors[0]

    @pytest.mark.parametrize(
        "assay_ontology_term_id, is_descendant",
        [("EFO:0022859", True), ("EFO:0022858", True), ("EFO:0030029", False), ("EFO:0002697", False)],
    )
    def test__validate_library_id_required_if_visium(self, assay_ontology_term_id, is_descendant):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()

        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        if is_descendant:
            # if spatial, `library_id` must exist
            validator._check_spatial_uns()
            assert len(validator.errors) == 0
            validator.reset()

            # if spatial, but missing from `uns`
            validator.adata.uns["spatial"].pop(visium_library_id)
            validator._check_spatial_uns()
            assert validator.errors == [
                f"uns['spatial'] must contain at least one key representing the library_id when {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE}."
            ]
        else:
            # if not spatial, MUST NOT define `library_id`
            validator.adata.uns["spatial"][visium_library_id] = {"images": []}
            validator._check_spatial_uns()
            # Report the most general top level error
            assert validator.errors == [
                "uns['spatial'] is only allowed when obs['assay_ontology_term_id'] is either a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2)"
            ]

    @pytest.mark.parametrize("library_id", [None, "invalid", 1, 1.0, True])
    def test__validate_library_id_type_error(self, library_id):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id] = library_id

        # Confirm library_id is identified as invalid.
        validator.validate_adata()
        assert validator.errors
        assert "uns['spatial'][library_id] must be a dictionary." in validator.errors[0]

    def test__validate_library_id_allowed_keys_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["invalid_key"] = True

        # Confirm additional key is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            "uns['spatial'][library_id] can only contain the keys 'images' and 'scalefactors'." in validator.errors[0]
        )

    def test__validate_images_required_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id].pop("images")

        # Confirm images is required.
        validator._check_spatial_uns()
        assert validator.errors
        assert "uns['spatial'][library_id] must contain the key 'images'." in validator.errors[0]

    @pytest.mark.parametrize(
        "assay_ontology_term_id, is_descendant",
        [("EFO:0010961", True), ("EFO:0022858", True), ("EFO:0030029", False), ("EFO:0002697", False)],
    )
    def test__validate_images_allowed_keys_error(self, assay_ontology_term_id, is_descendant):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["images"]["invalid_key"] = True

        # Confirm additional key is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            "uns['spatial'][library_id]['images'] can only contain the keys 'fullres' and 'hires'."
            in validator.errors[0]
        )

    def test__validate_images_hires_required_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["images"].pop("hires")

        # Confirm hires is required.
        validator._check_spatial_uns()
        assert validator.errors
        assert "uns['spatial'][library_id]['images'] must contain the key 'hires'." in validator.errors[0]

    def test__validate_images_fullres_optional(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["images"].pop("fullres")

        # Confirm fullres is optional.
        validator._check_spatial_uns()
        assert not validator.errors
        assert validator.warnings
        assert (
            "No uns['spatial'][library_id]['images']['fullres'] was found. "
            "It is STRONGLY RECOMMENDED that uns['spatial'][library_id]['images']['fullres'] is provided."
            in validator.warnings[0]
        )

    @pytest.mark.parametrize(
        "image_name, image_shape",
        [
            ("hires", (1, 2000, 4)),
            ("fullres", (1, 1, 4)),
        ],
    )
    def test__validate_images_image_last_dimension_4_ok(self, image_name, image_shape):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["images"][image_name] = np.zeros(image_shape, dtype=np.uint8)

        # Confirm image is valid.
        validator._check_spatial_uns()
        assert not validator.errors

    @pytest.mark.parametrize(
        "image_name, image_shape",
        [
            ("hires", (1, 2000, 3)),
            ("fullres", (1, 1, 3)),
        ],
    )
    def test__validate_images_image_ndarray_type_error(self, image_name, image_shape):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        # Defaults to float64.
        validator.adata.uns["spatial"][visium_library_id]["images"][image_name] = np.zeros(image_shape)

        # Confirm image is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            f"uns['spatial'][library_id]['images']['{image_name}'] must be of type numpy.uint8" in validator.errors[0]
        )

    @pytest.mark.parametrize(
        "image_name",
        [
            "hires",
            "fullres",
        ],
    )
    def test__validate_images_image_is_ndarray_error(self, image_name):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["images"][image_name] = "invalid"

        # Confirm image is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            f"uns['spatial'][library_id]['images']['{image_name}'] must be of numpy.ndarray type" in validator.errors[0]
        )

    @pytest.mark.parametrize(
        "image_name",
        [
            "hires",
            "fullres",
        ],
    )
    def test__validate_images_image_is_shape_error(self, image_name):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["images"][image_name] = np.zeros((1, 1), dtype=np.uint8)

        # Confirm image is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            f"uns['spatial'][library_id]['images']['{image_name}'] must have a length of 3 and either 3 (RGB color model "
            "for example) or 4 (RGBA color model for example) for its last dimension" in validator.errors[0]
        )

    @pytest.mark.parametrize(
        "assay_ontology_term_id, hi_res_size, image_max",
        [
            ("EFO:0022858", 2001, SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE),
            ("EFO:0022860", 4001, SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE_VISIUM_11MM),
        ],
    )
    def test__validate_images_hires_max_dimension_greater_than_error(
        self, assay_ontology_term_id, hi_res_size, image_max
    ):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        validator.adata.uns["spatial"][visium_library_id]["images"]["hires"] = np.zeros(
            (1, hi_res_size, 3), dtype=np.uint8
        )

        # Confirm hires is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors == [
            f"The largest dimension of uns['spatial'][library_id]['images']['hires'] must be {image_max} pixels, it has a largest dimension of {hi_res_size} pixels."
        ]

    @pytest.mark.parametrize(
        "assay_ontology_term_id, hi_res_size, size_requirement",
        [
            ("EFO:0022858", SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE, SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE),
            ("EFO:0022858", SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE_VISIUM_11MM, SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE),
            ("EFO:0022860", SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE, SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE_VISIUM_11MM),
            (
                "EFO:0022860",
                SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE_VISIUM_11MM,
                SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE_VISIUM_11MM,
            ),
        ],
    )
    def test__validate_images_hires_max_dimension(self, assay_ontology_term_id, hi_res_size, size_requirement):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        validator.adata.uns["spatial"][visium_library_id]["images"]["hires"] = np.zeros(
            (1, hi_res_size, 3), dtype=np.uint8
        )

        # Confirm hires is identified as invalid.
        validator.reset()
        validator._check_spatial_uns()
        if hi_res_size == size_requirement:
            assert validator.errors == []
        else:
            assert validator.errors == [
                f"The largest dimension of uns['spatial'][library_id]['images']['hires'] must be {size_requirement} pixels, it has a largest dimension of {hi_res_size} pixels."
            ]

    @pytest.mark.parametrize(
        "assay_ontology_term_id, hi_res_size, image_max",
        [
            ("EFO:0022858", 1999, SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE),
            ("EFO:0022860", 3999, SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE_VISIUM_11MM),
        ],
    )
    def test__validate_images_hires_max_dimension_less_than_error(self, assay_ontology_term_id, hi_res_size, image_max):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        validator.adata.uns["spatial"][visium_library_id]["images"]["hires"] = np.zeros(
            (1, hi_res_size, 3), dtype=np.uint8
        )

        # Confirm hires is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors == [
            f"The largest dimension of uns['spatial'][library_id]['images']['hires'] must be {image_max} pixels, it has a largest dimension of {hi_res_size} pixels."
        ]

    def test__validate_scalefactors_required_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id].pop("scalefactors")

        # Confirm scalefactors is required.
        validator._check_spatial_uns()
        assert validator.errors
        assert "uns['spatial'][library_id] must contain the key 'scalefactors'." in validator.errors[0]

    def test__validate_scalefactors_allowed_keys_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["scalefactors"]["invalid_key"] = True

        # Confirm additional key is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            "uns['spatial'][library_id]['scalefactors'] can only contain the keys "
            "'spot_diameter_fullres' and 'tissue_hires_scalef'." in validator.errors[0]
        )

    def test__validate_scalefactors_spot_diameter_fullres_required_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["scalefactors"].pop("spot_diameter_fullres")

        # Confirm spot_diameter_fullres is required.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            "uns['spatial'][library_id]['scalefactors'] must contain the key 'spot_diameter_fullres'."
            in validator.errors[0]
        )

    def test__validate_scalefactors_tissue_hires_scalef_required_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["scalefactors"].pop("tissue_hires_scalef")

        # Confirm tissue_hires_scalef is required.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            "uns['spatial'][library_id]['scalefactors'] must contain the key 'tissue_hires_scalef'."
            in validator.errors[0]
        )

    def test__validate_scalefactors_spot_diameter_fullres_is_float_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["scalefactors"]["spot_diameter_fullres"] = "invalid"

        # Confirm spot_diameter_fullres is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            "uns['spatial'][library_id]['scalefactors']['spot_diameter_fullres'] must be of type float"
            in validator.errors[0]
        )

    def test__validate_scalefactors_tissue_hires_scalef_is_float_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["scalefactors"]["tissue_hires_scalef"] = "invalid"

        # Confirm tissue_hires_scalef is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            "uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef'] must be of type float"
            in validator.errors[0]
        )

    @pytest.mark.parametrize("key", ["scalefactors", "images"])
    def test__validate_library_id_key_value_type_error(self, key):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id][key] = "invalid"

        # Confirm key type dict is required.
        validator._check_spatial_uns()
        assert validator.errors
        assert f"uns['spatial'][library_id]['{key}'] must be a dictionary." in validator.errors[0]

    def test__validate_assay_type_ontology_term_id_not_unique_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs.assay_ontology_term_id = ["EFO:0022858", "EFO:0030062"]

        # Confirm assay ontology term id is identified as invalid.
        validator._validate_spatial_assay_ontology_term_id()
        assert validator.errors
        assert (
            "When obs['assay_ontology_term_id'] is either a descendant"
            " of 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2), all observations must contain the same value."
        ) in validator.errors[0]

    def test__validate_assay_type_ontology_term_id_not_unique_ok(self, valid_adata):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = valid_adata  # "EFO:0009899" and "EFO:0009918"

        # Confirm assay ontology term id is considered valid.
        validator._validate_spatial_assay_ontology_term_id()
        assert not validator.errors

    @pytest.mark.parametrize(
        "assay_ontology_term_id, is_single",
        [
            (["EFO:0022858", "EFO:0030062"], True),
            (["EFO:0022858", "EFO:0030062"], False),
            ("EFO:0022858", False),
            ("EFO:0030062", True),
            ("EFO:0030062", False),
            ("EFO:0030062", False),
            ("EFO:0009899", True),  # Non-spatial
            ("EFO:0009899", False),  #  Non-spatial
        ],
    )
    def test__validate_tissue_position_forbidden(self, assay_ontology_term_id, is_single):
        validator: Validator = Validator()
        validator._set_schema_def()

        # Create Visium adata with tissue positions, update assay_ontology_term_id and is_single to trigger error.
        validator.adata = adata_visium.copy()
        validator.adata.obs.assay_ontology_term_id = assay_ontology_term_id
        validator.adata.uns["spatial"]["is_single"] = is_single
        validator.adata.obs["is_primary_data"] = False

        # Confirm tissue positions are not allowed.
        validator._validate_spatial_tissue_positions()
        assert len(validator.errors) == 3
        tissue_position_names = ["array_col", "array_row", "in_tissue"]
        for i, tissue_position_name in enumerate(tissue_position_names):
            assert (
                f"obs['{tissue_position_name}'] {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_FORBIDDEN}."
                in validator.errors[i]
            )

    @pytest.mark.parametrize("tissue_position_name", ["array_col", "array_row", "in_tissue"])
    def test__validate_tissue_position_required(self, tissue_position_name):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs.pop(tissue_position_name)

        # check visium
        validator.adata.obs["assay_ontology_term_id"] = "EFO:0022858"
        validator._check_spatial_obs()
        assert validator.errors
        assert (
            f"obs['{tissue_position_name}'] {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_REQUIRED}." in validator.errors[0]
        )
        validator.reset()

        # check visium descendant
        validator.adata.obs["assay_ontology_term_id"] = "EFO:0022860"
        validator._check_spatial_obs()
        assert validator.errors
        assert (
            f"obs['{tissue_position_name}'] {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_REQUIRED}." in validator.errors[0]
        )
        validator.reset()

    @pytest.mark.parametrize("assay_ontology_term_id", ["EFO:0022858", "EFO:0030062", "EFO:0022860"])
    def test__validate_tissue_position_not_required(self, assay_ontology_term_id):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_slide_seqv2.copy()
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        validator.adata.uns["spatial"]["is_single"] = False  # setting to false removes the requirement
        validator.adata.obs["is_primary_data"] = False
        validator._check_spatial_obs()
        assert not validator.errors

    @pytest.mark.parametrize("tissue_position_name", ["array_col", "array_row", "in_tissue"])
    def test__validate_tissue_position_int_error(self, tissue_position_name):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs[tissue_position_name] = 1.0

        # Confirm tissue_position is identified as invalid.
        validator._check_spatial_obs()
        assert validator.errors
        assert f"obs['{tissue_position_name}'] must be of int type" in validator.errors[0]

    @pytest.mark.parametrize("tissue_position_name", ["array_col", "array_row", "in_tissue"])
    def test__validate_tissue_position_nan_error(self, tissue_position_name):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs[tissue_position_name] = np.nan

        # Confirm tissue_position is identified as invalid.
        validator._check_spatial_obs()
        assert validator.errors[0] == f"obs['{tissue_position_name}'] {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_NOTNULL}."

    @pytest.mark.parametrize("assay_ontology_term_id", ["EFO:0022857", "EFO:0022860", "EFO:0022859"])
    @pytest.mark.parametrize("tissue_position_name, min", [("array_col", 0), ("array_row", 0), ("in_tissue", 0)])
    def test__validate_tissue_position_int_min_error(self, assay_ontology_term_id, tissue_position_name, min):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        validator.adata.obs[tissue_position_name] = min - 1

        # Confirm tissue_position is identified as invalid.
        validator._check_spatial_obs()
        assert (
            re.match(f"^obs\['{tissue_position_name}'\] must be (between )?{min} (and|or) [0-9]+", validator.errors[0])
            is not None
        )

    @pytest.mark.parametrize(
        "assay_ontology_term_id, tissue_position_name, tissue_position_max",
        [
            ("EFO:0022857", "array_col", 127),
            ("EFO:0022857", "array_row", 77),
            ("EFO:0022860", "array_col", 223),
            ("EFO:0022860", "array_row", 127),
            ("EFO:0022859", "array_col", 127),
            ("EFO:0022859", "array_row", 77),
            ("EFO:0022859", "in_tissue", 1),
        ],
    )
    def test__validate_tissue_position_int_max_error(
        self, assay_ontology_term_id, tissue_position_name, tissue_position_max
    ):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        validator.adata.obs[tissue_position_name] = tissue_position_max + 1

        # Confirm tissue_position is identified as invalid.
        validator._check_spatial_obs()
        assert (
            re.match(
                f"^obs\['{tissue_position_name}'\] must be (between )?[0-9]+ (and|or) {tissue_position_max}",
                validator.errors[0],
            )
            is not None
        )

    @pytest.mark.parametrize(
        "cell_type_ontology_term_id, in_tissue, assay_ontology_term_id",
        [
            # MUST be unknown when in_tissue = 0 and assay_ontology_term_id = Visium Spatial Gene Expression v2
            ("unknown", 0, "EFO:0022858"),
            # MUST be unknown when in_tissue = 0 and assay_ontology_term_id = Visium CytAssist Spatial Gene Expression, 11mm
            ("unknown", 0, "EFO:0022860"),
            # MUST be unknown when in_tissue = 0 and assay_ontology_term_id = Visium Spatial Gene Expression V1
            # valid CL term is ok when in_tissue = 1 and assay_ontology_term_id = Visium CytAssist Spatial Gene Expression, 11mm
            (["unknown", "CL:0000066"], [0, 1], ["EFO:0022857", "EFO:0022860"]),
            # normal CL term for in_tissue = 1 and assay_ontology_term_id = 10x 3' v2
            ("CL:0000066", 1, "EFO:0009899"),
        ],
    )
    def test__validate_cell_type_ontology_term_id_ok(
        self, cell_type_ontology_term_id, in_tissue, assay_ontology_term_id
    ):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs.cell_type_ontology_term_id = cell_type_ontology_term_id
        validator.adata.obs.in_tissue = in_tissue
        validator.adata.obs.assay_ontology_term_id = assay_ontology_term_id

        # Confirm cell type is valid.
        validator._validate_spatial_cell_type_ontology_term_id()
        assert not validator.errors

    @pytest.mark.parametrize(
        "cell_type_ontology_term_id, in_tissue, assay_ontology_term_id",
        [
            # MUST be unknown when in_tissue = 0 and assay_ontology_term_id = Visium Spatial Gene Expression
            ("CL:0000066", 0, "EFO:0022858"),
            (["CL:0000066", "unknown"], [0, 1], ["EFO:0022858", "EFO:0022858"]),
            # MUST be unknown when in_tissue = 0 and assay_ontology_term_id = Visium CytAssist Spatial Gene Expression, 11mm
            ("CL:0000066", 0, "EFO:0022860"),
            # MUST be unknown when in_tissue = 0 and assay_ontology_term_id = Visium Spatial Gene Expression V1
            ("CL:0000066", 0, "EFO:0022857"),
        ],
    )
    def test__validate_cell_type_ontology_term_id_error(
        self, cell_type_ontology_term_id, in_tissue, assay_ontology_term_id
    ):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs.cell_type_ontology_term_id = cell_type_ontology_term_id
        validator.adata.obs.in_tissue = in_tissue
        validator.adata.obs.assay_ontology_term_id = assay_ontology_term_id

        # Confirm errors.
        validator._validate_spatial_cell_type_ontology_term_id()
        assert validator.errors
        assert (
            f"obs['cell_type_ontology_term_id'] must be 'unknown' when {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_IN_TISSUE_0}."
            in validator.errors[0]
        )

    def test__validate_embeddings_non_nans(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator._visium_and_is_single_true_matrix_size = 2

        # invalidate spatial embeddings with NaN value
        validator.adata.obsm["spatial"][0, 1] = np.nan
        # Confirm spatial is valid.
        validator.validate_adata()
        assert validator.errors == ["ERROR: adata.obsm['spatial'] contains at least one NaN value."]


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
        assert "in dataframe 'obs' contains 2 categorical types. Only one type is allowed." in validator.errors[0]
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
            (from_array(np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)), "dense", True),
            # Test case with float values in a dense matrix
            (from_array(np.array([[1.1, 2.2, 3.3], [4.4, 5.5, 6.6]])), "dense", False),
            # Test case with integer values in a sparse matrix (CSR format)
            (from_array(sparse.csr_matrix([[1, 0, 3], [0, 5, 0]], dtype=np.float32)), "csr", True),
            # Test case with float values in a sparse matrix (CSC format)
            (from_array(sparse.csc_matrix([[1.1, 0, 3.3], [0, 5.5, 0]])), "csc", False),
            # Test case with mixed integer and float values in a dense matrix
            (from_array(np.array([[1, 2.2, 3], [4.4, 5, 6.6]])), "dense", False),
        ],
    )
    def test_has_valid_raw(self, data, matrix_format, expected_result):
        validator = self.create_validator(data, matrix_format)
        assert validator._has_valid_raw() == expected_result

    def test_has_valid_raw_with_unknown_format(self):
        # a matrix with unknown format should be invalid
        data = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
        validator = self.create_validator(data, "unknown")
        assert validator._has_valid_raw() is False
        assert validator.errors == [f"Unknown encoding for matrix X. {ERROR_SUFFIX_SPARSE_FORMAT}"]


class TestValidateAnndataRawCounts:
    """
    Paired Assays (measures both accessibility and gene expression) require raw counts to be present and validated.
    Unpaired Assays (measures only accessibility) do not require raw counts to be present.

    This is the purview of the AnnData Validator, not the AtacValidator, which is primarily focused on the fragment file.
    """

    def test_paired_requires_raw_validation(self, atac_anndata, tmpdir):
        # 10x multiome (EFO:0030059) is paired (both ATAC and RNA single cell sequencing)
        atac_anndata.obs["assay_ontology_term_id"] = ["EFO:0030059"] * 3

        # use a valid count matrix (as dask array)
        X = atac_anndata.X

        # validate with AnnData Validator
        validator = Validator(ignore_labels=True)
        validator._set_schema_def()

        # do validation with a valid count matrix
        atac_anndata.X = from_array(X.astype("float32"))
        validator.adata = atac_anndata
        validator.reset()
        validator._validate_raw()
        assert validator.errors == []

        # do validation with an invalid count matrix
        atac_anndata.X = from_array(zeros(X.shape).astype("float32"))
        validator.adata = atac_anndata
        validator.reset()
        validator._validate_raw()
        assert len(validator.errors) > 0

    def test_unpaired_skips_raw_validation(self, atac_anndata, tmpdir):
        # scATAC-seq (EFO:0010891) is unpaired paired
        atac_anndata.obs["assay_ontology_term_id"] = ["EFO:0010891"] * 3

        # remove matrix - it shouldn't be required
        del atac_anndata.X

        # check that validation passes even without matrix
        validator = Validator(ignore_labels=True)
        validator._set_schema_def()
        validator.adata = atac_anndata
        validator.reset()
        validator._validate_raw()
        assert validator.errors == []
