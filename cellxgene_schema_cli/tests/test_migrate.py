import unittest
from tempfile import TemporaryDirectory
from unittest.mock import patch

import anndata
from cellxgene_schema.migrate import migrate
from fixtures.examples_validate import adata_with_lables_unmigrated


class TestMigrate(unittest.TestCase):
    def test_migrate(self):
        with TemporaryDirectory() as tmp, patch("cellxgene_schema.migrate.DEPRECATED_FEATURE_IDS", ["DUMMY"]):
            result_h5ad = tmp + "result.h5ad"
            test_h5ad = tmp + "test.h5ad"

            adata_with_lables_unmigrated.copy().write_h5ad(test_h5ad, compression="gzip")
            migrate(
                input_file=test_h5ad,
                output_file=result_h5ad,
                collection_id="",
                dataset_id="",
            )

            adata = anndata.read_h5ad(result_h5ad)
            assert not any(adata.var.index.isin(["DUMMY"]))
            raw_adata = anndata.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs)
            assert not any(raw_adata.var.index.isin(["DUMMY"]))
