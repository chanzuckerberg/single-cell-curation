import pytest
from cellxgene_schema.remove_labels import AnnDataLabelRemover
from fixtures.examples_validate import (
    adata,
    adata_with_labels,
)
from pandas.testing import assert_frame_equal


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
