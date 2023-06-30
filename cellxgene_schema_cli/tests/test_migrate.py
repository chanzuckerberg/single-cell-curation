import unittest
from tempfile import NamedTemporaryFile

from cellxgene_schema.migrate import migrate
from cellxgene_schema.validate import validate
from fixtures.examples_validate import h5ad_valid


class TestMigrate(unittest.TestCase):
    def test_migrate(self):
        with NamedTemporaryFile() as tmp:
            migrate(
                input_file=h5ad_valid,
                output_file=tmp.name,
                collection_id="",
                dataset_id="",
            )

            success, errors, is_seurat_convertible = validate(tmp.name)

            self.assertTrue(success)
            self.assertListEqual(errors, [])
            self.assertTrue(is_seurat_convertible)
