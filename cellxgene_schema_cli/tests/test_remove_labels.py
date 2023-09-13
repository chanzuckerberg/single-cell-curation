import unittest

from cellxgene_schema.remove_labels import AnnDataLabelRemover
from fixtures.examples_validate import (
    adata,
    adata_with_labels,
)
from pandas.testing import assert_frame_equal


class TestRemoveLabels(unittest.TestCase):
    def setUp(self):
        self.anndata_label_remover = AnnDataLabelRemover(adata_with_labels.copy())
        self.adata_no_labels = adata

    def test_remove_labels(self):
        self.anndata_label_remover.adata.raw = adata_with_labels
        self.anndata_label_remover.adata.raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        self.anndata_label_remover.remove_labels()
        adata_labels_removed = self.anndata_label_remover.adata
        assert_frame_equal(adata_labels_removed.obs, self.adata_no_labels.obs)
        assert_frame_equal(adata_labels_removed.var, self.adata_no_labels.var)
        assert_frame_equal(adata_labels_removed.raw.var, self.adata_no_labels.raw.var)
        self.assertDictEqual(dict(adata_labels_removed.uns), dict(self.adata_no_labels.uns))

    def test_remove_labels_no_raw(self):
        self.anndata_label_remover.remove_labels()
        adata_labels_removed = self.anndata_label_remover.adata

        self.assertIsNone(adata_labels_removed.raw)
        assert_frame_equal(adata_labels_removed.obs, self.adata_no_labels.obs)
        assert_frame_equal(adata_labels_removed.var, self.adata_no_labels.var)
        self.assertDictEqual(dict(adata_labels_removed.uns), dict(self.adata_no_labels.uns))
