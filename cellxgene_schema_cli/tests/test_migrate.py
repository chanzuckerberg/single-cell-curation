from tempfile import TemporaryDirectory
from unittest.mock import patch

import anndata
from cellxgene_schema.migrate import migrate
from fixtures.examples_validate import adata_with_lables_unmigrated


class TestMigrate:
    def test_migrate(self):
        test_ONTOLOGY_TERM_MAPS = {
            "assay": {
                "assay:1": "assay:2",
            },
            "cell_type": {
                "cell_type:1": "cell_type:2",
            },
            "development_stage": {
                "development_stage:1": "development_stage:2",
            },
            "disease": {
                "disease:1": "disease:2",
            },
            "organism": {
                "organism:1": "organism:2",
            },
            "self_reported_ethnicity": {
                "sre:1": "sre:2",
            },
            "sex": {"sex:1": "sex:2"},
            "tissue": {"tissue:1": "tissue:2"},
        }
        with TemporaryDirectory() as tmp, patch("cellxgene_schema.migrate.DEPRECATED_FEATURE_IDS", ["DUMMY"]), patch(
            "cellxgene_schema.migrate.ONTOLOGY_TERM_MAPS", test_ONTOLOGY_TERM_MAPS
        ):
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
