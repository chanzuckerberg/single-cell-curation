import unittest

from cellxgene_schema.remove_labels import remove_labels_inplace
from fixtures.examples_validate import (
    adata_with_labels,
    adata,
)
from pandas.testing import assert_frame_equal


class TestRemoveLabels(unittest.TestCase):

    def setUp(self):
        self.adata_with_labels = adata_with_labels.copy()
        self.adata_no_labels = adata

    def test_remove_labels_inplace(self):
        self.adata_with_labels.raw = adata_with_labels
        self.adata_with_labels.raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        remove_labels_inplace(self.adata_with_labels)

        assert_frame_equal(self.adata_with_labels.obs, self.adata_no_labels.obs)
        assert_frame_equal(self.adata_with_labels.var, self.adata_no_labels.var)
        assert_frame_equal(self.adata_with_labels.raw.var, self.adata_no_labels.raw.var)

    def test_remove_labels_inplace_no_raw(self):
        remove_labels_inplace(self.adata_with_labels)

        self.assertIsNone(self.adata_with_labels.raw)
        assert_frame_equal(self.adata_with_labels.obs, self.adata_no_labels.obs)
        assert_frame_equal(self.adata_with_labels.var, self.adata_no_labels.var)
