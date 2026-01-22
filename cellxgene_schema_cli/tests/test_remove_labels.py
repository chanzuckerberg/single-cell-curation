import anndata as ad
import numpy as np
import pandas as pd
import pytest
from cellxgene_schema.remove_labels import AnnDataLabelRemover
from fixtures.examples_validate import (
    adata,
    adata_with_labels,
)
from pandas.testing import assert_frame_equal
from scipy import sparse


@pytest.fixture
def remove_labels_setup():
    anndata_label_remover = AnnDataLabelRemover(adata_with_labels.copy())
    adata_no_labels = adata
    return anndata_label_remover, adata_no_labels


class TestRemoveLabels:
    def test_remove_labels(self, remove_labels_setup):
        anndata_label_remover, adata_no_labels = remove_labels_setup
        anndata_label_remover.adata.raw = adata_with_labels
        anndata_label_remover.adata.raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        anndata_label_remover.remove_labels()
        adata_labels_removed = anndata_label_remover.adata
        assert_frame_equal(adata_labels_removed.obs, adata_no_labels.obs)
        assert_frame_equal(adata_labels_removed.var, adata_no_labels.var)
        assert_frame_equal(adata_labels_removed.raw.var, adata_no_labels.raw.var)
        assert dict(adata_labels_removed.uns) == dict(adata_no_labels.uns)

    def test_remove_labels_no_raw(self, remove_labels_setup):
        anndata_label_remover, adata_no_labels = remove_labels_setup
        anndata_label_remover.remove_labels()
        adata_labels_removed = anndata_label_remover.adata

        assert adata_labels_removed.raw is None
        assert_frame_equal(adata_labels_removed.obs, adata_no_labels.obs)
        assert_frame_equal(adata_labels_removed.var, adata_no_labels.var)
        assert dict(adata_labels_removed.uns) == dict(adata_no_labels.uns)


class TestRemoveLabelsGeneticPerturbations:
    """Test removal of reserved_keys from genetic_perturbations"""

    def test_remove_derived_genomic_regions_from_genetic_perturbations(self):
        """Test that derived_genomic_regions is removed from genetic_perturbations"""
        # Create adata with genetic_perturbations containing derived_genomic_regions
        X = sparse.csr_matrix(np.zeros((2, 10), dtype=np.float32))
        obs = pd.DataFrame(index=["cell1", "cell2"])
        var = pd.DataFrame(index=[f"gene{i}" for i in range(10)])

        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.uns["genetic_perturbations"] = {
            "guide1": {
                "role": "targeting",
                "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
                "protospacer_adjacent_motif": "3' NGG",
                "derived_genomic_regions": ["16:75647615-75647633(-)"],  # Should be removed
            },
            "guide2": {
                "role": "targeting",
                "protospacer_sequence": "GGGCCCTCCGGGAAGATGG",
                "protospacer_adjacent_motif": "3' NGG",
                "derived_genomic_regions": ["1:100-200(+)", "2:300-400(-)"],  # Should be removed
            },
        }

        remover = AnnDataLabelRemover(adata)
        remover.remove_labels()

        # Check that derived_genomic_regions was removed from both guides
        assert "derived_genomic_regions" not in adata.uns["genetic_perturbations"]["guide1"]
        assert "derived_genomic_regions" not in adata.uns["genetic_perturbations"]["guide2"]
        # Check that other keys are preserved
        assert "role" in adata.uns["genetic_perturbations"]["guide1"]
        assert "protospacer_sequence" in adata.uns["genetic_perturbations"]["guide1"]
        assert "protospacer_adjacent_motif" in adata.uns["genetic_perturbations"]["guide1"]

    def test_remove_derived_features_from_genetic_perturbations(self):
        """Test that derived_features is removed from genetic_perturbations"""
        X = sparse.csr_matrix(np.zeros((2, 10), dtype=np.float32))
        obs = pd.DataFrame(index=["cell1", "cell2"])
        var = pd.DataFrame(index=[f"gene{i}" for i in range(10)])

        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.uns["genetic_perturbations"] = {
            "guide1": {
                "role": "targeting",
                "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
                "protospacer_adjacent_motif": "3' NGG",
                "derived_features": {"ENSG00000123456": "GENE1"},  # Should be removed
            },
        }

        remover = AnnDataLabelRemover(adata)
        remover.remove_labels()

        # Check that derived_features was removed
        assert "derived_features" not in adata.uns["genetic_perturbations"]["guide1"]
        # Check that other keys are preserved
        assert "role" in adata.uns["genetic_perturbations"]["guide1"]
        assert "protospacer_sequence" in adata.uns["genetic_perturbations"]["guide1"]

    def test_remove_both_derived_keys_from_genetic_perturbations(self):
        """Test that both derived_genomic_regions and derived_features are removed"""
        X = sparse.csr_matrix(np.zeros((2, 10), dtype=np.float32))
        obs = pd.DataFrame(index=["cell1", "cell2"])
        var = pd.DataFrame(index=[f"gene{i}" for i in range(10)])

        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.uns["genetic_perturbations"] = {
            "guide1": {
                "role": "targeting",
                "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
                "protospacer_adjacent_motif": "3' NGG",
                "derived_genomic_regions": ["16:75647615-75647633(-)"],
                "derived_features": {"ENSG00000123456": "GENE1"},
            },
        }

        remover = AnnDataLabelRemover(adata)
        remover.remove_labels()

        # Check that both derived keys were removed
        assert "derived_genomic_regions" not in adata.uns["genetic_perturbations"]["guide1"]
        assert "derived_features" not in adata.uns["genetic_perturbations"]["guide1"]
        # Check that required keys are preserved
        assert "role" in adata.uns["genetic_perturbations"]["guide1"]
        assert "protospacer_sequence" in adata.uns["genetic_perturbations"]["guide1"]
        assert "protospacer_adjacent_motif" in adata.uns["genetic_perturbations"]["guide1"]

    def test_remove_reserved_keys_from_multiple_perturbations(self):
        """Test removal from multiple genetic perturbations"""
        X = sparse.csr_matrix(np.zeros((2, 10), dtype=np.float32))
        obs = pd.DataFrame(index=["cell1", "cell2"])
        var = pd.DataFrame(index=[f"gene{i}" for i in range(10)])

        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.uns["genetic_perturbations"] = {
            "guide1": {
                "role": "targeting",
                "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
                "protospacer_adjacent_motif": "3' NGG",
                "derived_genomic_regions": ["16:75647615-75647633(-)"],
            },
            "guide2": {
                "role": "control",
                "protospacer_sequence": "TTTTTTTTTTTTTTTTTTTT",
                "protospacer_adjacent_motif": "3' NGG",
                "derived_features": {"ENSG00000123456": "GENE1"},
            },
            "guide3": {
                "role": "targeting",
                "protospacer_sequence": "GGGCCCTCCGGGAAGATGG",
                "protospacer_adjacent_motif": "3' NGG",
                # No derived keys - should remain unchanged
            },
        }

        remover = AnnDataLabelRemover(adata)
        remover.remove_labels()

        # Check guide1
        assert "derived_genomic_regions" not in adata.uns["genetic_perturbations"]["guide1"]
        assert "role" in adata.uns["genetic_perturbations"]["guide1"]

        # Check guide2
        assert "derived_features" not in adata.uns["genetic_perturbations"]["guide2"]
        assert "role" in adata.uns["genetic_perturbations"]["guide2"]

        # Check guide3 (no derived keys, should be unchanged)
        assert "role" in adata.uns["genetic_perturbations"]["guide3"]
        assert "protospacer_sequence" in adata.uns["genetic_perturbations"]["guide3"]

    def test_genetic_perturbations_without_reserved_keys(self):
        """Test that genetic_perturbations without reserved keys are not affected"""
        X = sparse.csr_matrix(np.zeros((2, 10), dtype=np.float32))
        obs = pd.DataFrame(index=["cell1", "cell2"])
        var = pd.DataFrame(index=[f"gene{i}" for i in range(10)])

        adata = ad.AnnData(X=X, obs=obs, var=var)
        original_perturbations = {
            "guide1": {
                "role": "targeting",
                "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
                "protospacer_adjacent_motif": "3' NGG",
            },
        }
        adata.uns["genetic_perturbations"] = original_perturbations.copy()

        remover = AnnDataLabelRemover(adata)
        remover.remove_labels()

        # Should remain unchanged
        assert adata.uns["genetic_perturbations"] == original_perturbations

    def test_genetic_perturbations_not_present(self):
        """Test that adata without genetic_perturbations is handled gracefully"""
        X = sparse.csr_matrix(np.zeros((2, 10), dtype=np.float32))
        obs = pd.DataFrame(index=["cell1", "cell2"])
        var = pd.DataFrame(index=[f"gene{i}" for i in range(10)])

        adata = ad.AnnData(X=X, obs=obs, var=var)
        # No genetic_perturbations in uns

        remover = AnnDataLabelRemover(adata)
        # Should not raise an error
        remover.remove_labels()

        assert "genetic_perturbations" not in adata.uns


class TestRemoveReservedRecursive:
    """Test recursive removal of reserved_keys with nested structures"""

    def test_remove_reserved_recursive_nested_dict(self):
        """Test recursive removal of reserved_keys in nested dictionaries using dict structure"""
        X = sparse.csr_matrix(np.zeros((2, 10), dtype=np.float32))
        obs = pd.DataFrame(index=["cell1", "cell2"])
        var = pd.DataFrame(index=[f"gene{i}" for i in range(10)])

        adata = ad.AnnData(X=X, obs=obs, var=var)
        # Create a nested structure in uns with reserved keys
        adata.uns["test_dict"] = {
            "level1": {
                "level2": {
                    "reserved_key": "should_be_removed",
                    "normal_key": "should_remain",
                },
                "level2_alt": {
                    "reserved_key2": "should_be_removed",
                    "normal_key2": "should_remain",
                },
            },
            "normal_top_level": "should_remain",
        }

        remover = AnnDataLabelRemover(adata)
        # Use dict structure to specify recursive removal paths
        # This tests that the method can recurse into nested dicts and remove keys from nested levels
        reserved_dict = {
            "level1": {
                "level2": ["reserved_key"],
                "level2_alt": ["reserved_key2"],
            }
        }
        remover._remove_reserved(adata.uns["test_dict"], reserved_dict)

        # Check that reserved keys were removed from nested levels
        assert "reserved_key" not in adata.uns["test_dict"]["level1"]["level2"]
        assert "reserved_key2" not in adata.uns["test_dict"]["level1"]["level2_alt"]
        # Check that normal keys remain
        assert "normal_key" in adata.uns["test_dict"]["level1"]["level2"]
        assert "normal_key2" in adata.uns["test_dict"]["level1"]["level2_alt"]
        assert "normal_top_level" in adata.uns["test_dict"]

    def test_remove_reserved_with_key_pattern(self):
        """Test removal of keys matching a pattern"""
        X = sparse.csr_matrix(np.zeros((2, 10), dtype=np.float32))
        obs = pd.DataFrame(index=["cell1", "cell2"])
        var = pd.DataFrame(index=[f"gene{i}" for i in range(10)])

        adata = ad.AnnData(X=X, obs=obs, var=var)
        # Create a dict with keys that match a pattern
        adata.uns["test_dict"] = {
            "guide1": {"role": "targeting"},
            "guide2": {"role": "control"},
            "other_key": {"data": "value"},
            "guide3": {"role": "targeting"},
        }

        remover = AnnDataLabelRemover(adata)
        # Test key_pattern matching - keys starting with "guide"
        reserved_dict = {
            "key_pattern": "^guide",
        }
        remover._remove_reserved(adata.uns["test_dict"], reserved_dict)

        # Check that pattern-matched keys were removed
        assert "guide1" not in adata.uns["test_dict"]
        assert "guide2" not in adata.uns["test_dict"]
        assert "guide3" not in adata.uns["test_dict"]
        # Check that non-matching key remains
        assert "other_key" in adata.uns["test_dict"]

    def test_remove_reserved_mixed_list_and_dict(self):
        """Test removal with list (top-level only) and dict (recursive) reserved_keys"""
        X = sparse.csr_matrix(np.zeros((2, 10), dtype=np.float32))
        obs = pd.DataFrame(index=["cell1", "cell2"])
        var = pd.DataFrame(index=[f"gene{i}" for i in range(10)])

        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.uns["test_dict"] = {
            "key1": "value1",
            "key2": {
                "nested_key": "nested_value",
                "reserved_nested": "should_be_removed",
            },
            "reserved_top": "should_be_removed",
        }

        remover = AnnDataLabelRemover(adata)
        # Test with a list of reserved keys (removes only from top level)
        reserved_list = ["reserved_top"]
        remover._remove_reserved(adata.uns["test_dict"], reserved_list)

        # Check top-level removal
        assert "reserved_top" not in adata.uns["test_dict"]
        # Check that nested reserved key is still there (list doesn't recurse)
        assert "reserved_nested" in adata.uns["test_dict"]["key2"]
        # Check preservation
        assert "key1" in adata.uns["test_dict"]
        assert "nested_key" in adata.uns["test_dict"]["key2"]

        # Now test recursive removal with dict structure
        reserved_dict = {"key2": ["reserved_nested"]}
        remover._remove_reserved(adata.uns["test_dict"], reserved_dict)

        # Now nested reserved key should be removed
        assert "reserved_nested" not in adata.uns["test_dict"]["key2"]
        assert "nested_key" in adata.uns["test_dict"]["key2"]
