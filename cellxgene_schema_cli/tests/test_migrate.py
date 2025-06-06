from tempfile import TemporaryDirectory
from unittest.mock import patch

import anndata
from cellxgene_schema.migrate import migrate
from fixtures.examples_validate import adata_with_labels_unmigrated


class TestMigrate:
    def test_migrate(self):
        test_ONTOLOGY_OBS_TERM_MAPS = {
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
            "self_reported_ethnicity": {
                "sre:1": "sre:2",
            },
            "sex": {"sex:1": "sex:2"},
            "tissue": {"tissue:1": "tissue:2"},
        }
        test_ONTOLOGY_UNS_TERM_MAPS = {
            "organism": {
                "NCBITaxon:9606": "NCBITaxon:9606 migrated",
            }
        }
        with TemporaryDirectory() as tmp, patch("cellxgene_schema.migrate.DEPRECATED_FEATURE_IDS", ["DUMMY"]), patch(
            "cellxgene_schema.migrate.ONTOLOGY_TERM_OBS_MAPS", test_ONTOLOGY_OBS_TERM_MAPS
        ), patch("cellxgene_schema.migrate.GENCODE_MAPPER", {"ENSSASG00005000004": "ENSSASG00005000004_NEW"}), patch(
            "cellxgene_schema.migrate.ONTOLOGY_TERM_UNS_MAPS", test_ONTOLOGY_UNS_TERM_MAPS
        ):
            result_h5ad = tmp + "result.h5ad"
            test_h5ad = tmp + "test.h5ad"
            adata_with_labels_unmigrated.copy().write_h5ad(test_h5ad, compression="gzip")

            # Verify regular adata is what we expect before migration
            assert any(adata_with_labels_unmigrated.var.index.isin(["DUMMY"]))
            assert any(adata_with_labels_unmigrated.var.index.isin(["ENSSASG00005000004"]))
            assert not any(adata_with_labels_unmigrated.var.index.isin(["ENSSASG00005000004_NEW"]))
            assert adata_with_labels_unmigrated.X.shape == (2, 2)

            # Verify raw adata is what we expect before migration
            adata_raw_with_labels_unmigrated = anndata.AnnData(
                adata_with_labels_unmigrated.raw.X,
                var=adata_with_labels_unmigrated.raw.var,
                obs=adata_with_labels_unmigrated.obs,
            )
            assert any(adata_raw_with_labels_unmigrated.var.index.isin(["DUMMY"]))
            assert any(adata_raw_with_labels_unmigrated.var.index.isin(["ENSSASG00005000004"]))
            assert not any(adata_raw_with_labels_unmigrated.var.index.isin(["ENSSASG00005000004_NEW"]))
            assert adata_raw_with_labels_unmigrated.X.shape == (2, 2)
            assert "organism_ontology_term_id" not in adata_raw_with_labels_unmigrated.uns

            migrate(
                input_file=test_h5ad,
                output_file=result_h5ad,
                collection_id="",
                dataset_id="",
            )

            # Verify regular adata is what we expect after migration
            adata = anndata.read_h5ad(result_h5ad)
            assert not any(adata.var.index.isin(["DUMMY"]))
            assert not any(adata.var.index.isin(["ENSSASG00005000004"]))
            assert any(adata.var.index.isin(["ENSSASG00005000004_NEW"]))
            assert adata.X.shape == (2, 1)

            # Verify raw adata is what we expect after migration
            raw_adata = anndata.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs)
            assert not any(raw_adata.var.index.isin(["DUMMY"]))
            assert not any(raw_adata.var.index.isin(["ENSSASG00005000004"]))
            assert any(raw_adata.var.index.isin(["ENSSASG00005000004_NEW"]))
            assert raw_adata.X.shape == (2, 1)

            # Verify organism ontology term was mapped
            assert adata.uns["organism_ontology_term_id"] == "NCBITaxon:9606 migrated"
