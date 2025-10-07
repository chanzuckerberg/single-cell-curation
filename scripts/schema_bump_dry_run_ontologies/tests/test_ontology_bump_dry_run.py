import json
import os
from tempfile import NamedTemporaryFile
from unittest.mock import Mock, patch

import pytest

from scripts.schema_bump_dry_run_ontologies import ontology_bump_dry_run

FIXTURES_ROOT = os.path.join(os.path.dirname(__file__), "fixtures")


@pytest.fixture
def mock_ontology_metadata():  # type: ignore
    return {
        "EFO:0000001": {
            "label": "non-deprecated term",
            "deprecated": False,
            "consider": None,
            "term_tracker": None,
            "comments": [],
        },
        "EFO:0000002": {
            "label": "obsolete term with replacement",
            "deprecated": True,
            "replaced_by": "EFO:0000001",
            "consider": None,
            "term_tracker": None,
            "comments": [],
        },
        "EFO:0000003": {
            "label": "obsolete term without replacement, with comment",
            "deprecated": True,
            "consider": None,
            "term_tracker": None,
            "comments": ["note: replace with EFO:0000001"],
        },
        "EFO:0000004": {
            "label": "obsolete term without replacement, with consider",
            "deprecated": True,
            "consider": ["EFO:0000001", "EFO:0000010"],
            "term_tracker": None,
            "comments": [],
        },
        "EFO:0000005": {
            "label": "obsolete term without replacement, with comment and consider",
            "deprecated": True,
            "consider": ["EFO:0000001", "EFO:0000010"],
            "term_tracker": None,
            "comments": ["Consider removing this term entirely", "Or using one from a different ontology"],
        },
        "EFO:0000006": {
            "label": "obsolete term with replacement from a different ontology",
            "deprecated": True,
            "replaced_by": "CL:0000006",
            "consider": None,
            "term_tracker": None,
            "comments": ["Comment adding context to replacement"],
        },
        "EFO:0000007": {
            "label": "obsolete term with replacement and comment",
            "deprecated": True,
            "replaced_by": "EFO:0000070",
            "consider": None,
            "term_tracker": None,
            "comments": ["Comment adding context to replacement"],
        },
        "EFO:0000008": {
            "label": "obsolete term with replacement and term tracker",
            "deprecated": True,
            "replaced_by": "EFO:0000080",
            "consider": None,
            "term_tracker": "www.fake-github-link.com/repo/example-1",
            "comments": [],
        },
        "EFO:0000009": {
            "label": "obsolete term with considers and term tracker",
            "deprecated": True,
            "consider": ["EFO:0000090"],
            "term_tracker": "www.fake-github-link.com/repo/example-2",
            "comments": [],
        },
        "HANCESTRO:0000001": {
            "label": "non-deprecated term",
            "deprecated": False,
            "consider": None,
            "term_tracker": None,
            "comments": [],
        },
        "HANCESTRO:0000002": {
            "label": "obsolete term with replacement",
            "deprecated": True,
            "replaced_by": "HANCESTRO:0000003",
            "consider": None,
            "term_tracker": None,
            "comments": [],
        },
    }


@pytest.fixture
def mock_ontology_parser(mock_ontology_metadata):  # type: ignore
    mock_ontology_parser = Mock()
    mock_ontology_parser.is_term_deprecated.side_effect = lambda term_id: mock_ontology_metadata[term_id]["deprecated"]
    mock_ontology_parser.get_term_metadata.side_effect = lambda term_id: mock_ontology_metadata[term_id]
    mock_ontology_parser.get_term_replacement.side_effect = lambda term_id: mock_ontology_metadata[term_id].get(
        "replaced_by", None
    )
    return Mock(return_value=mock_ontology_parser)


@pytest.fixture
def public_datasets():  # type: ignore
    return [
        {
            "collection_id": "public_coll_0",
            "dataset_id": "public_ds_0",
            "assay": [{"ontology_term_id": "EFO:0000001"}, {"ontology_term_id": "unknown"}],
        }
    ]


@pytest.fixture
def private_collections():  # type: ignore
    return [
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


@pytest.fixture
def mock_fetch_private_dataset(private_collections):  # type: ignore
    def mock_fetch_private_dataset(*args):  # type: ignore
        collection_id = args[2]
        dataset_id = args[3]
        for c in private_collections:
            if c["collection_id"] == collection_id:
                for ds in c["datasets"]:
                    if ds["dataset_id"] == dataset_id:
                        return ds

    return Mock(side_effect=mock_fetch_private_dataset)


@pytest.fixture
def expected_replaced_by_map():  # type: ignore
    return {
        "assay": {},
        "cell_type": {},
        "development_stage": {},
        "disease": {},
        "organism": {},
        "self_reported_ethnicity": {},
        "sex": {},
        "tissue": {},
    }


@pytest.fixture
def setup(mock_ontology_parser, public_datasets, private_collections, mock_fetch_private_dataset, expected_replaced_by_map):  # type: ignore
    patcher = patch.multiple(
        "scripts.schema_bump_dry_run_ontologies.ontology_bump_dry_run",
        OntologyParser=mock_ontology_parser,
        get_headers=Mock(),
        fetch_public_datasets=Mock(return_value=public_datasets),
        fetch_private_collections=Mock(return_value=private_collections),
        fetch_private_dataset=mock_fetch_private_dataset,
    )
    patcher.start()
    yield public_datasets, private_collections, expected_replaced_by_map
    patcher.stop()


class TestSplitTerm:
    def test_split_term_no_delimiter(self):
        assert ontology_bump_dry_run.split_term("HANCESTRO:0000001") == ["HANCESTRO:0000001"]

    @pytest.mark.parametrize("delimiter", [",", "||"])
    def test_split_term_delimiter(self, delimiter):
        assert ontology_bump_dry_run.split_term(f"HANCESTRO:0000001{delimiter}HANCESTRO:0000002") == [
            "HANCESTRO:0000001",
            "HANCESTRO:0000002",
        ]

    def test_split_term_mixed_delimiters_raises(self):
        with pytest.raises(ValueError):
            ontology_bump_dry_run.split_term("HANCESTRO:0000001,HANCESTRO:0000002||HANCESTRO:0000003")


class TestOntologyBumpDryRun:
    def test_no_deprecated_terms_no_revisions(self, setup):  # type: ignore
        public_datasets, private_collections, expected_replaced_by_map = setup
        private_collections.pop()
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            f"{FIXTURES_ROOT}/no_deprecated_terms_no_revisions_expected", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            assert list(expected) == list(tmp)
            assert expected_replaced_by_map == json.load(tmp_json)

    def test_no_deprecated_terms_with_open_revisions(self, setup):  # type: ignore
        _, _, expected_replaced_by_map = setup
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            f"{FIXTURES_ROOT}/no_deprecated_terms_with_open_revisions", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            assert list(expected) == list(tmp)
            assert expected_replaced_by_map == json.load(tmp_json)

    def test_with_multi_delimited_onto_id_list(self, setup):  # type: ignore
        public_datasets, private_collections, expected_replaced_by_map = setup
        private_collections.pop()
        public_datasets[0]["self_reported_ethnicity"] = [{"ontology_term_id": "HANCESTRO:0000001,HANCESTRO:0000002"}]
        expected_replaced_by_map["self_reported_ethnicity"]["HANCESTRO:0000002"] = "HANCESTRO:0000003"
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            f"{FIXTURES_ROOT}/with_comma_delimited_onto_id_list_expected", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            assert list(expected) == list(tmp)
            assert expected_replaced_by_map == json.load(tmp_json)

    def test_with_replaced_by(self, setup):  # type: ignore
        public_datasets, private_collections, expected_replaced_by_map = setup
        public_datasets[0]["assay"].extend([{"ontology_term_id": "EFO:0000002"}, {"ontology_term_id": "EFO:0000007"}])
        private_collections[0]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000002"})
        private_collections[1]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000002"})
        expected_replaced_by_map["assay"]["EFO:0000002"] = "EFO:0000001"
        expected_replaced_by_map["assay"]["EFO:0000007"] = "EFO:0000070"
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            f"{FIXTURES_ROOT}/with_replaced_by_expected", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            assert list(expected) == list(tmp)
            assert expected_replaced_by_map == json.load(tmp_json)

    def test_with_comments_and_considers(self, setup):  # type: ignore
        public_datasets, private_collections, expected_replaced_by_map = setup
        public_datasets[0]["assay"].append({"ontology_term_id": "EFO:0000003"})
        private_collections[0]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000004"})
        private_collections[1]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000005"})
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            f"{FIXTURES_ROOT}/with_comments_and_considers_expected", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            assert list(expected) == list(tmp)
            assert expected_replaced_by_map == json.load(tmp_json)

    def test_with_replaced_by_diff_ontology(self, setup):  # type: ignore
        public_datasets, private_collections, expected_replaced_by_map = setup
        public_datasets[0]["assay"].append({"ontology_term_id": "EFO:0000006"})
        private_collections[0]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000006"})
        private_collections[1]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000006"})
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            f"{FIXTURES_ROOT}/with_replaced_by_diff_ontology_expected", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            assert list(expected) == list(tmp)
            assert expected_replaced_by_map == json.load(tmp_json)

    def test_group_datasets_by_collection(self, setup):  # type: ignore
        public_datasets, private_collections, expected_replaced_by_map = setup
        public_datasets[0]["assay"].append({"ontology_term_id": "EFO:0000002"})
        public_datasets.append(
            {
                "collection_id": "public_coll_0",
                "dataset_id": "public_ds_1",
                "assay": [{"ontology_term_id": "EFO:0000002"}],
            }
        )
        private_collections[0]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000005"})
        private_collections[0]["datasets"].append(
            {
                "dataset_id": "private_ds_1",
                "assay": [{"ontology_term_id": "EFO:0000005"}],
                "processing_status": "SUCCESS",
            }
        )
        private_collections[1]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000006"})
        private_collections[1]["datasets"].append(
            {
                "dataset_id": "revision_ds_0",
                "assay": [{"ontology_term_id": "EFO:0000006"}],
                "processing_status": "SUCCESS",
            }
        )
        expected_replaced_by_map["assay"]["EFO:0000002"] = "EFO:0000001"
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            f"{FIXTURES_ROOT}/group_datasets_by_collection_expected", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            assert list(expected) == list(tmp)
            assert expected_replaced_by_map == json.load(tmp_json)

    def test_term_tracker(self, setup):  # type: ignore
        public_datasets, private_collections, expected_replaced_by_map = setup
        public_datasets[0]["assay"].append({"ontology_term_id": "EFO:0000008"})
        private_collections[0]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000009"})
        private_collections[1]["datasets"][0]["assay"].append({"ontology_term_id": "EFO:0000008"})
        expected_replaced_by_map["assay"]["EFO:0000008"] = "EFO:0000080"
        with NamedTemporaryFile() as tmp, NamedTemporaryFile() as tmp_json, open(
            f"{FIXTURES_ROOT}/with_term_tracker_expected", "rb"
        ) as expected:
            ontology_bump_dry_run.dry_run(tmp.name, tmp_json.name)
            assert list(expected) == list(tmp)
            assert expected_replaced_by_map == json.load(tmp_json)
