import hashlib
import os
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
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_FORBIDDEN,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_IN_TISSUE_0,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_REQUIRED,
    Validator,
    validate,
)
from cellxgene_schema.write_labels import AnnDataLabelAppender
from fixtures.examples_validate import adata as adata_valid
from fixtures.examples_validate import (
    adata_minimal,
    adata_slide_seqv2,
    adata_spatial_is_single_false,
    adata_visium,
    adata_with_labels,
    good_obs,
    good_obsm,
    good_uns,
    good_uns_with_visium_spatial,
    good_var,
    h5ad_invalid,
    h5ad_valid,
    visium_library_id,
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
            1061,
            2821,
            1797,
            3822,
        ]
        expected_dict = dict(zip(ids, gene_lengths))
        assert label_writer._get_mapping_dict_feature_length(ids) == expected_dict

        # Bad
        ids = ["NO_GENE"]
        with pytest.raises(KeyError):
            label_writer._get_mapping_dict_feature_id(ids)

    def test_get_dictionary_mapping_feature_type(self, label_writer):
        # Good
        ids = [
            "ERCC-00002",
            "ENSG00000127603",
            "ENSMUSG00000059552",
            "ENSSASG00005000004",
        ]
        # values derived from csv
        gene_types = [
            "synthetic",
            "protein_coding",
            "protein_coding",
            "protein_coding",
        ]
        expected_dict = dict(zip(ids, gene_types))
        assert label_writer._get_mapping_dict_feature_type(ids) == expected_dict

        # Bad
        ids = ["NO_GENE"]
        with pytest.raises(KeyError):
            label_writer._get_mapping_dict_feature_type(ids)

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


class TestCheckSpatial:
    def test__validate_spatial_visium_ok(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.visium_and_is_single_true_matrix_size = 2
        # Confirm spatial is valid.
        validator.validate_adata()
        assert not validator.errors

    def test__validate_from_file(self):
        """Testing compatibility with SparseDatset types in Anndata"""
        validator: Validator = Validator()
        validator._set_schema_def()
        with tempfile.TemporaryDirectory() as temp_dir:
            file_path = os.path.join(temp_dir, "slide_seqv2.h5ad")
            adata_slide_seqv2.write_h5ad(file_path)
            # Confirm spatial is valid.
            validator.validate_adata(file_path)
        assert not validator.errors

    def test__validate_spatial_visium_dense_matrix_ok(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.visium_and_is_single_true_matrix_size = 2
        validator.adata.X = validator.adata.X.toarray()
        validator.adata.raw = validator.adata.copy()
        validator.adata.raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        # Confirm spatial is valid.
        validator.validate_adata()
        assert not validator.errors

    def test__validate_spatial_visium_and_is_single_false_ok(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"] = {"is_single": False}
        del validator.adata.obsm["spatial"]
        # Format adata.obs into valid shape for Visium and is_single False.
        validator.adata.obs.pop("array_col")
        validator.adata.obs.pop("array_row")
        validator.adata.obs.pop("in_tissue")
        validator.adata.obs["is_primary_data"] = False
        # Confirm spatial is valid.
        validator.validate_adata()
        assert not validator.errors

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
        assert validator.errors
        assert (
            "A dict in uns['spatial'] is required for obs['assay_ontology_term_id'] values 'EFO:0010961' (Visium Spatial Gene Expression) and 'EFO:0030062' (Slide-seqV2)."
            in validator.errors[0]
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
        assert len(validator.errors) == 1
        assert (
            "uns['spatial'] is only allowed for obs['assay_ontology_term_id'] values "
            "'EFO:0010961' (Visium Spatial Gene Expression) and 'EFO:0030062' (Slide-seqV2)." in validator.errors[0]
        )

    def test__validate_spatial_required_if_visium(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns = good_uns.copy()

        # Confirm spatial is required for Visium.
        validator._check_spatial_uns()
        assert len(validator.errors) == 1
        assert (
            "A dict in uns['spatial'] is required for obs['assay_ontology_term_id'] values "
            "'EFO:0010961' (Visium Spatial Gene Expression) and 'EFO:0030062' (Slide-seqV2)." in validator.errors[0]
        )

    def test__validate_spatial_required_if_slide_seqV2(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_slide_seqv2.copy()
        validator.adata.uns = good_uns.copy()

        # Confirm spatial is required for Slide-seqV2.
        validator._check_spatial_uns()
        assert len(validator.errors) == 1
        assert (
            "A dict in uns['spatial'] is required for obs['assay_ontology_term_id'] values "
            "'EFO:0010961' (Visium Spatial Gene Expression) and 'EFO:0030062' (Slide-seqV2)." in validator.errors[0]
        )

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

    def test__validate_is_single_required_visium_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"].pop("is_single")

        # Confirm is_single is identified as required.
        validator._check_spatial_uns()
        assert validator.errors
        assert "uns['spatial'] must contain the key 'is_single'." in validator.errors[0]

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

    def test__validate_library_id_required_if_visium(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"].pop(visium_library_id)

        # Confirm library_id is identified as required.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            f"uns['spatial'] must contain at least one key representing the library_id when {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE}."
            in validator.errors[0]
        )

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

    def test__validate_images_allowed_keys_error(self):
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

    def test__validate_images_hires_max_dimension_greater_than_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["images"]["hires"] = np.zeros((1, 2001, 3), dtype=np.uint8)

        # Confirm hires is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            "The largest dimension of uns['spatial'][library_id]['images']['hires'] must be 2000 pixels"
            in validator.errors[0]
        )

    def test__validate_images_hires_max_dimension_less_than_error(self):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.uns["spatial"][visium_library_id]["images"]["hires"] = np.zeros((1, 1999, 3), dtype=np.uint8)

        # Confirm hires is identified as invalid.
        validator._check_spatial_uns()
        assert validator.errors
        assert (
            "The largest dimension of uns['spatial'][library_id]['images']['hires'] must be 2000 pixels"
            in validator.errors[0]
        )

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
        validator.adata.obs.assay_ontology_term_id = ["EFO:0010961", "EFO:0030062"]

        # Confirm assay ontology term id is identified as invalid.
        validator._validate_spatial_assay_ontology_term_id()
        assert validator.errors
        assert (
            "When obs['assay_ontology_term_id'] is either 'EFO:0010961' (Visium Spatial Gene Expression) or "
            "'EFO:0030062' (Slide-seqV2), all observations must contain the same value."
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
            (["EFO:0010961", "EFO:0030062"], True),
            (["EFO:0010961", "EFO:0030062"], False),
            ("EFO:0010961", False),
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

        validator._check_spatial_obs()
        assert validator.errors
        assert (
            f"obs['{tissue_position_name}'] {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_REQUIRED}." in validator.errors[0]
        )

    @pytest.mark.parametrize("assay_ontology_term_id", ["EFO:0010961", "EFO:0030062"])
    def test__validate_tissue_position_not_required(self, assay_ontology_term_id):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_slide_seqv2.copy()
        validator.adata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        validator.adata.uns["spatial"]["is_single"] = False
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

    @pytest.mark.parametrize(
        "tissue_position_name, min, error_message_token",
        [
            ("array_col", 0, "between 0 and 127"),
            ("array_row", 0, "between 0 and 77"),
            ("in_tissue", 0, "0 or 1"),
        ],
    )
    def test__validate_tissue_position_int_min_error(self, tissue_position_name, min, error_message_token):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs[tissue_position_name] = min - 1

        # Confirm tissue_position is identified as invalid.
        validator._check_spatial_obs()
        assert validator.errors
        assert f"obs['{tissue_position_name}'] must be {error_message_token}" in validator.errors[0]

    @pytest.mark.parametrize(
        "tissue_position_name, max, error_message_token",
        [
            ("array_col", 127, "between 0 and 127"),
            ("array_row", 77, "between 0 and 77"),
            ("in_tissue", 1, "0 or 1"),
        ],
    )
    def test__validate_tissue_position_int_max_error(self, tissue_position_name, max, error_message_token):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs[tissue_position_name] = max + 1

        # Confirm tissue_position is identified as invalid.
        validator._check_spatial_obs()
        assert validator.errors
        assert f"obs['{tissue_position_name}'] must be {error_message_token}" in validator.errors[0]

    @pytest.mark.parametrize(
        "cell_type_ontology_term_id, in_tissue",
        [("unknown", 0), (["unknown", "CL:0000066"], [0, 1]), ("CL:0000066", 1)],
    )
    def test__validate_cell_type_ontology_term_id_ok(self, cell_type_ontology_term_id, in_tissue):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs.cell_type_ontology_term_id = cell_type_ontology_term_id
        validator.adata.obs.in_tissue = in_tissue

        # Confirm cell type is valid.
        validator._validate_spatial_cell_type_ontology_term_id()
        assert not validator.errors

    @pytest.mark.parametrize(
        "cell_type_ontology_term_id, in_tissue",
        [
            ("CL:0000066", 0),
            (["CL:0000066", "unknown"], [0, 1]),
        ],
    )
    def test__validate_cell_type_ontology_term_id_error(self, cell_type_ontology_term_id, in_tissue):
        validator: Validator = Validator()
        validator._set_schema_def()
        validator.adata = adata_visium.copy()
        validator.adata.obs.cell_type_ontology_term_id = cell_type_ontology_term_id
        validator.adata.obs.in_tissue = in_tissue

        # Confirm errors.
        validator._validate_spatial_cell_type_ontology_term_id()
        assert validator.errors
        assert (
            f"obs['cell_type_ontology_term_id'] must be 'unknown' when {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_IN_TISSUE_0}."
            in validator.errors[0]
        )


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

        # Visium datasets are not Seurat-convertible
        self.validation_helper(sparse_matrix_with_zero)
        self.validator.adata.obs = adata_visium.obs.copy()
        self.validator._validate_seurat_convertibility()
        assert len(self.validator.warnings) == 1
        assert not self.validator.is_seurat_convertible


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
            (np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32), "dense", True),
            # Test case with float values in a dense matrix
            (np.array([[1.1, 2.2, 3.3], [4.4, 5.5, 6.6]]), "dense", False),
            # Test case with integer values in a sparse matrix (CSR format)
            (sparse.csr_matrix([[1, 0, 3], [0, 5, 0]], dtype=np.float32), "csr", True),
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
        data = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
        validator = self.create_validator(data, "unknown")
        with pytest.raises(AssertionError):
            validator._has_valid_raw()
