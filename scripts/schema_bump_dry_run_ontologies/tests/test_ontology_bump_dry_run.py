import json
from tempfile import NamedTemporaryFile
from unittest import TestCase
from unittest.mock import Mock, patch

from scripts.schema_bump_dry_run_ontologies import ontology_bump_dry_run


class TestOntologyBumpDryRun(TestCase):
    def setUp(self):
        mock_onto_map = {
            "EFO": {
                "EFO:0000001": {
                    "label": "non-deprecated term",
                    "deprecated": False,
                },
                "EFO:0000002": {
                    "label": "obsolete term with replacement",
                    "deprecated": True,
                    "replaced_by": "EFO:0000001",
                },
                "EFO:0000003": {
                    "label": "obsolete term without replacement, with comment",
                    "deprecated": True,
                    "comments": ["note: replace with EFO:0000001"],
                },
                "EFO:0000004": {
                    "label": "obsolete term without replacement, with consider",
                    "deprecated": True,
                    "consider": ["EFO:0000001", "EFO:0000010"],
                },
                "EFO:0000005": {
                    "label": "obsolete term without replacement, with comment and consider",
                    "deprecated": True,
                    "comments": ["Consider removing this term entirely", "Or using one from a different ontology"],
                    "consider": ["EFO:0000001", "EFO:0000010"],
                },
                "EFO:0000006": {
                    "label": "obsolete term with replacement from a different ontology",
                    "deprecated": True,
                    "replaced_by": "CL:0000006",
                    "comments": ["Comment adding context to replacement"],
                },
                "EFO:0000007": {
                    "label": "obsolete term with replacement and comment",
                    "deprecated": True,
                    "replaced_by": "EFO:0000070",
                    "comments": ["Comment adding context to replacement"],
                },
            }
        }
        self.public_datasets = [
            {
                "collection_id": "public_coll_0",
                "dataset_id": "public_ds_0",
                "assay": [{"ontology_term_id": "EFO:0000001"}, {"ontology_term_id": "unknown"}],
            }
        ]
        self.private_collections = [
            {
                "collection_id": "private_coll_0",
                "revision_of": None,
                "datasets": [
                    {
                        "dataset_id": "private_ds_0",
                        "assay": [{"ontology_term_id": "EFO:0000001"}],
                        "processing_status": "SUCCESS",
                    },
                    {"dataset_id": "failed_ds_0", "processing_status": "FAILURE"},
                ],
            },
            {
                "collection_id": "public_coll_0_revision",
                "revision_of": "public_coll_0",
                "datasets": [
                    {
                        "dataset_id": "public_ds_0",
                        "assay": [{"ontology_term_id": "EFO:0000001"}],
                        "processing_status": "SUCCESS",
                    }
                ],
            },
        ]
        self.patcher = patch.multiple(
            "scripts.schema_bump_dry_run_ontologies.ontology_bump_dry_run",
            load_ontology_map=Mock(return_value=mock_onto_map),
            get_headers=Mock(),
            fetch_public_datasets=Mock(return_value=self.public_datasets),
            fetch_private_collections=Mock(return_value=self.private_collections),
            fetch_private_dataset=Mock(side_effect=self.mock_fetch_private_dataset),
        )
        self.mock_load_ontology_map = self.patcher.start()
        self.expected_replaced_by_map = {
            "assay": {},
            "cell_type": {},
            "development_stage": {},
            "disease": {},
            "organism": {},
            "self_reported_ethnicity": {},
            "sex": {},
            "tissue": {},
        }

    def mock_fetch_private_dataset(self, *args):
        collection_id = args[2]
        dataset_id = args[3]
        for c in self.private_collections:
            if c["collection_id"] == collection_id:
                for ds in c["datasets"]:
                    if ds["dataset_id"] == dataset_id:
                        return ds

    def test_no_deprecated_terms_no_revisions(self):
        self.private_collections.pop()
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            "fixtures/no_deprecated_terms_no_revisions_expected", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            self.assertListEqual(list(expected), list(tmp))
            self.assertDictEqual(self.expected_replaced_by_map, json.load(tmp_json))

    def test_no_deprecated_terms_with_open_revisions(self):
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            "fixtures/no_deprecated_terms_with_open_revisions", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            self.assertListEqual(list(expected), list(tmp))
            self.assertDictEqual(self.expected_replaced_by_map, json.load(tmp_json))

    def test_with_replaced_by(self):
        self.public_datasets[0]["assay"].extend(
            [{"ontology_term_id": "EFO:0000002"}, {"ontology_term_id": "EFO:0000007"}]
        )
        self.private_collections[0]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000002"})
        self.private_collections[1]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000002"})
        self.expected_replaced_by_map["assay"]["EFO:0000002"] = "EFO:0000001"
        self.expected_replaced_by_map["assay"]["EFO:0000007"] = "EFO:0000070"
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            "fixtures/with_replaced_by_expected", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            self.assertListEqual(list(expected), list(tmp))
            self.assertDictEqual(self.expected_replaced_by_map, json.load(tmp_json))

    def test_with_comments_and_considers(self):
        self.public_datasets[0]["assay"].append({"ontology_term_id": "EFO:0000003"})
        self.private_collections[0]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000004"})
        self.private_collections[1]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000005"})
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            "fixtures/with_comments_and_considers_expected", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            self.assertListEqual(list(expected), list(tmp))
            self.assertDictEqual(self.expected_replaced_by_map, json.load(tmp_json))

    def test_with_replaced_by_diff_ontology(self):
        self.public_datasets[0]["assay"].append({"ontology_term_id": "EFO:0000006"})
        self.private_collections[0]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000006"})
        self.private_collections[1]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000006"})
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            "fixtures/with_replaced_by_diff_ontology_expected", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            self.assertListEqual(list(expected), list(tmp))
            self.assertDictEqual(self.expected_replaced_by_map, json.load(tmp_json))

    def test_group_datasets_by_collection(self):
        self.public_datasets[0]["assay"].append({"ontology_term_id": "EFO:0000002"})
        self.public_datasets.append(
            {
                "collection_id": "public_coll_0",
                "dataset_id": "public_ds_1",
                "assay": [{"ontology_term_id": "EFO:0000002"}],
            }
        )
        self.private_collections[0]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000005"})
        self.private_collections[0]["datasets"].append(
            {
                "dataset_id": "private_ds_1",
                "assay": [{"ontology_term_id": "EFO:0000005"}],
                "processing_status": "SUCCESS",
            }
        )
        self.private_collections[1]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000006"})
        self.private_collections[1]["datasets"].append(
            {
                "dataset_id": "revision_ds_0",
                "assay": [{"ontology_term_id": "EFO:0000006"}],
                "processing_status": "SUCCESS",
            }
        )
        self.expected_replaced_by_map["assay"]["EFO:0000002"] = "EFO:0000001"
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            "fixtures/group_datasets_by_collection_expected", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            self.assertListEqual(list(expected), list(tmp))
            self.assertDictEqual(self.expected_replaced_by_map, json.load(tmp_json))
