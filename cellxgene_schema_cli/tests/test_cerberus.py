import fixtures.examples_validate as examples
import numpy
import pandas as pd
import pytest
from cellxgene_schema.ex_cerberus import get_validator
from test_schema_compliance import save_and_read_adata


@pytest.fixture
def adata():
    return examples.adata.copy()


@pytest.fixture
def validator():
    return get_validator()


def validate_anndata(adata, validator):
    document = dict(
        adata=adata,
        obs=adata.obs,
        var=adata.var,
        obsm=adata.obsm,
        obsp=adata.obsp,
        varm=adata.varm,
        varp=adata.varp,
        uns=adata.uns,
        raw=adata.raw,
    )
    return validator.validate(document, normalize=False)


class TestObsm:
    """
    Fail cases for adata.obsm
    """

    def test_obsm_values_ara_numpy(self, adata, validator):
        """
        values in obsm MUST be a numpy.ndarray
        """

        obsm = adata.obsm
        obsm["X_tsne"] = pd.DataFrame(obsm["X_umap"], index=adata.obs_names)
        obsm["X_umap"] = numpy.array(obsm["X_umap"])
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == [
            "ERROR: All embeddings have to be of 'numpy.ndarray' type, "
            "'adata.obsm['X_tsne']' is <class 'pandas.core.frame.DataFrame'>')."
        ]

    def test_obsm_values_infinity(self, adata, validator):
        """
        values in obsm cannot have any infinity values
        """

        adata.obsm["X_umap"][0:100, 1] = numpy.inf
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == [
            "ERROR: adata.obsm['X_umap'] contains positive infinity or negative infinity values."
        ]

    def test_obsm_values_str(self, adata, validator):
        """
        values in obsm must be numerical types, strings are not valid
        """

        obsm = adata.obsm
        all_string = numpy.full(obsm["X_umap"].shape, "test")
        obsm["X_umap"] = all_string
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == [
            "ERROR: adata.obsm['X_umap'] has an invalid data type. It should be float, integer, or unsigned "
            "integer of any precision (8, 16, 32, or 64 bits)."
        ]

    def test_obsm_values_nan(self, adata, validator):
        """
        values in obsm cannot all be NaN
        """

        obsm = adata.obsm
        # It's okay if only one value is NaN
        obsm["X_umap"][0:100, 1] = numpy.nan
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == []

        # It's not okay if all values are NaN
        all_nan = numpy.full(obsm["X_umap"].shape, numpy.nan)
        obsm["X_umap"] = all_nan
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == ["ERROR: adata.obsm['X_umap'] contains all NaN values."]

    def test_obsm_values_no_X_embedding__non_spatial_dataset(self, adata, validator):
        adata.obsm["harmony"] = adata.obsm["X_umap"]
        adata.uns["default_embedding"] = "harmony"
        del adata.obsm["X_umap"]
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == [
            "ERROR: At least one embedding in 'obsm' has to have a key with an 'X_' prefix.",
        ]
        assert validator.is_spatial is False
        assert validator.warnings == [
            "WARNING: Dataframe 'var' only has 4 rows. Features SHOULD NOT be filtered from expression matrix.",
            "WARNING: Embedding key in 'adata.obsm' harmony does not start with X_ and thus will not be available in Explorer",
            "WARNING: Validation of raw layer was not performed due to current errors, try again after fixing current errors.",
        ]

    @pytest.mark.parametrize("assay_ontology_term_id", ["EFO:0010961", "EFO:0030062"])
    def test_obsm_values_no_X_embedding__spatial_dataset(self, adata, assay_ontology_term_id):
        adata.obsm["harmony"] = adata.obsm["X_umap"]
        adata.uns["default_embedding"] = "harmony"
        del adata.obsm["X_umap"]
        adata.obs["assay_ontology_term_id"] = assay_ontology_term_id
        adata.obs["suspension_type"] = "na"
        adata.obs.loc[:, ["suspension_type"]] = adata.obs.astype("category")
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == []
        assert validator.is_spatial is True

    def test_obsm_values_warn_start_with_X(self, adata, validator):
        adata.obsm["harmony"] = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names)
        assert not validate_anndata(adata, validator)
        assert validator.warnings == [
            "WARNING: Dataframe 'var' only has 4 rows. Features SHOULD NOT be filtered from expression matrix.",
            "WARNING: Embedding key in 'adata.obsm' harmony does not start with X_ and thus will not be available in Explorer",
            "WARNING: All embeddings have to be of 'numpy.ndarray' type, 'adata.obsm['harmony']' is <class 'pandas.core.frame.DataFrame'>').",
        ]

    def test_obsm_values_suffix_is_forbidden(self, adata, validator):
        adata.obsm["X_spatial"] = adata.obsm["X_umap"]
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == ["ERROR: Embedding key in 'adata.obsm' X_spatial cannot be used."]

    def test_obsm_values_suffix_start_with_number(self, adata, validator):
        adata.obsm["X_3D"] = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names)
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == [
            "ERROR: Suffix for embedding key in 'adata.obsm' X_3D does not match the regex pattern ^[a-zA-Z][a-zA-Z0-9_.-]*$.",
            "ERROR: All embeddings have to be of 'numpy.ndarray' type, 'adata.obsm['X_3D']' is <class 'pandas.core.frame.DataFrame'>').",
        ]

    def test_obsm_values_key_start_with_number(self, adata, validator):
        adata.obsm["3D"] = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names)
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == [
            "ERROR: Embedding key in 'adata.obsm' 3D does not match the regex pattern ^[a-zA-Z][a-zA-Z0-9_.-]*$."
        ]
        assert validator.warnings == [
            "WARNING: Dataframe 'var' only has 4 rows. Features SHOULD NOT be filtered from expression matrix.",
            "WARNING: Embedding key in 'adata.obsm' 3D does not start with X_ and thus will not be available in Explorer",
            "WARNING: All embeddings have to be of 'numpy.ndarray' type, 'adata.obsm['3D']' is <class 'pandas.core.frame.DataFrame'>').",
            "WARNING: Validation of raw layer was not performed due to current errors, try again after fixing current errors.",
        ]

    def test_obsm_suffix_name_valid(self, adata, validator):
        """
        Suffix after X_ must be at least 1 character long
        """

        adata.obsm["X_"] = adata.obsm["X_umap"]
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == [
            "ERROR: Suffix for embedding key in 'adata.obsm' X_ does not match the regex pattern ^[a-zA-Z][a-zA-Z0-9_.-]*$."
        ]

    def test_obsm_key_name_whitespace(self, adata, validator):
        """
        Embedding keys beginning with X_ with whitespace are not valid
        """

        obsm = adata.obsm
        obsm["X_ umap"] = obsm["X_umap"]
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == [
            "ERROR: Suffix for embedding key in 'adata.obsm' X_ umap does not match the regex pattern ^[a-zA-Z][a-zA-Z0-9_.-]*$.",
        ]

        del obsm["X_ umap"]
        obsm["u m a p"] = obsm["X_umap"]
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == [
            "ERROR: Embedding key in 'adata.obsm' u m a p does not match the regex pattern ^[a-zA-Z][a-zA-Z0-9_.-]*$."
        ]

    def test_obsm_suffix_has_special_characters_valid(self, adata, validator):
        adata.obsm["X_umap_MinDist_0.2_N_Neighbors-15"] = adata.obsm["X_umap"]
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == []

    def test_obsm_shape_one_column(self, adata, validator):
        """
        Curators MUST annotate one or more two-dimensional (m >= 2) embeddings
        """
        # Makes 1 column array

        adata.obsm["X_umap"] = numpy.delete(adata.obsm["X_umap"], 0, 1)
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == [
            "ERROR: All embeddings must have as many rows as cells, and "
            "at least two columns. 'adata.obsm['X_umap']' has shape "
            "of '(2, 1)'."
        ]

    def test_obsm_shape_same_rows_and_columns(self, adata, validator):
        """
        The number of rows must be equal to the number of columns
        """

        obsm = adata.obsm
        # Create a 3 row array
        arr1 = numpy.array([0, 0])
        arr2 = numpy.array([0, 0])
        arr3 = numpy.array([0, 0])
        three_row_array = numpy.vstack((arr1, arr2, arr3))
        with pytest.raises(ValueError):
            obsm["X_umap"] = three_row_array
            assert validate_anndata(adata, validator)

    def test_obsm_size_zero(self, adata, validator):
        """
        size of obsm an key cannot be zero.
        """

        adata = adata
        adata.obsm["badsize"] = numpy.empty((2, 0))
        adata = save_and_read_adata(adata)
        assert not validate_anndata(adata, validator)
        assert validator.errors["obsm"] == [
            "ERROR: The size of the ndarray stored for a 'adata.obsm['badsize']' MUST NOT be zero.",
        ]
