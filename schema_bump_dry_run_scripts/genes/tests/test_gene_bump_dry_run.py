from unittest.mock import Mock, patch

import pytest
from cellxgene_schema.ontology import SupportedOrganisms

from schema_bump_dry_run_scripts.genes.gene_bump_dry_run import (
    compare_genes,
    generate_deprecated_private,
    generate_deprecated_public,
    generate_report,
    get_diff_map,
)


def test_get_diff_map(tmp_path):
    for key in SupportedOrganisms:
        with open(f"{tmp_path}/{key.name}_diff.txt", "w") as fp:
            fp.write("test")
    with patch("schema_bump_dry_run_scripts.genes.gene_bump_dry_run.ONTOLOGY_DIR", tmp_path):
        diff_map = get_diff_map()
    assert len(diff_map) == 4
    for key in diff_map:
        assert key in ["NCBITaxon:9606", "NCBITaxon:10090", "NCBITaxon:2697049", "NCBITaxon:32630"]


@pytest.fixture
def sample_report_data():
    return {
        "deprecated_public": {
            "existing_collection": {
                "num_datasets": 3,
                "deprecated_genes": [
                    "term1",
                    "term2",
                    "term3",
                    "term1",
                    "term2",
                    "term3",
                    "term1",
                    "term2",
                    "term3",
                    "term1",
                    "term2",
                    "term3",
                    "term1",
                    "term2",
                    "term3",
                    "term1",
                    "term2",
                    "term3",
                ],
            },
            "new_collection": {"num_datasets": 1, "deprecated_genes": ["term4"]},
            "collection4": {"num_datasets": 2, "deprecated_genes": ["gene1", "gene2"]},
        },
        "open_revisions": {
            "collection3": {"revision_of": "collection4", "num_datasets": 2, "deprecated_genes": ["gene1", "gene2"]},
            "collection5": {"num_datasets": 1, "deprecated_genes": ["gene3"]},
        },
        "non_auto_migrated": ["collection6", "collection7"],
    }


def test_generate_report(sample_report_data):
    expected_report = """## Deprecated Terms in Public Datasets:

Collection ID: existing_collection
Number of Affected Datasets: 3
Deprecated Terms: 
    term1, term2, term3, term1, term2, term3, term1, term2, term3, term1, term2,
    term3, term1, term2, term3, term1, term2, term3

Collection ID: new_collection
Number of Affected Datasets: 1
Deprecated Terms: 
    term4

Collection ID: collection4
Number of Affected Datasets: 2
Deprecated Terms: 
    gene1, gene2

## Deprecated Genes in Private Datasets:

Collection ID: collection3
Note--In A Revision of: collection4
Number of Affected Datasets: 2
Deprecated Terms: 
    gene1, gene2

Collection ID: collection5
Number of Affected Datasets: 1
Deprecated Terms: 
    gene3

## The Following Public Collections Will Not Be Auto-Migrated Due To Having an Open Revision:
collection6
collection7
"""
    assert generate_report(sample_report_data) == expected_report


def test_generate_report__with_empty_data():
    empty_data = {"deprecated_public": {}, "open_revisions": {}, "non_auto_migrated": []}
    expected_report = """## Deprecated Terms in Public Datasets:

## Deprecated Genes in Private Datasets:

## The Following Public Collections Will Not Be Auto-Migrated Due To Having an Open Revision:
"""
    assert generate_report(empty_data) == expected_report


def get_genes(dataset):
    # Mock implementation for testing purposes
    return dataset.get("get_genes_mock_response", [])


@pytest.fixture
def sample_diff_map():
    return {
        "organism-1": ["deprecated:1-1", "deprecated:1-2", "deprecated:1-3"],
        "organism-2": ["deprecated:2-1", "deprecated:2-2"],
    }


@pytest.fixture
def sample_deprecated_datasets():
    return {
        "existing_collection": {
            "genes": {"current:1-1", "current:1-2", "deprecated:1-1", "deprecated:1-2"},
            "num_datasets": 2,
            "deprecated_genes": {"deprecated:1-1", "deprecated:1-2"},
        }
    }


@pytest.fixture
def sample_dataset():
    return {
        "dataset_id": "dataset1",
        "collection_id": "new_collection",
        "organism": [{"ontology_term_id": "organism-1"}],
        "get_genes_mock_response": ["current:1-A", "current:1-B"],
    }


def test_compare_genes__with_no_deprecated_genes(sample_diff_map, sample_deprecated_datasets, sample_dataset):
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    expected_is_deprecated_genes_found = False

    with patch("schema_bump_dry_run_scripts.genes.gene_bump_dry_run.get_genes", get_genes), patch(
        "schema_bump_dry_run_scripts.genes.gene_bump_dry_run.logger"
    ) as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            sample_dataset, sample_diff_map, sample_deprecated_datasets
        )

    mock_logger.info.assert_called_once_with("Dataset dataset1 has no deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_empty_diff_map(sample_deprecated_datasets, sample_dataset):
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    expected_is_deprecated_genes_found = False

    with patch("schema_bump_dry_run_scripts.genes.gene_bump_dry_run.get_genes", get_genes), patch(
        "schema_bump_dry_run_scripts.genes.gene_bump_dry_run.logger"
    ) as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            sample_dataset, {}, sample_deprecated_datasets
        )

    mock_logger.info.assert_called_once_with(f"Dataset {sample_dataset['dataset_id']} has no deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_existing_collection_and_deprecated_genes(
    sample_diff_map, sample_deprecated_datasets, sample_dataset
):
    sample_dataset.update(
        **{
            "collection_id": "existing_collection",
            "get_genes_mock_response": ["deprecated:1-2", "deprecated:1-3", "current:1-2", "current:1-3"],
        }
    )
    expected_deprecated_datasets = {
        "existing_collection": {
            "num_datasets": 3,
            "deprecated_genes": {"deprecated:1-1", "deprecated:1-2", "deprecated:1-3"},
            "genes": {
                "deprecated:1-1",
                "deprecated:1-2",
                "deprecated:1-3",
                "current:1-1",
                "current:1-2",
                "current:1-3",
            },
        }
    }
    expected_is_deprecated_genes_found = True
    with patch("schema_bump_dry_run_scripts.genes.gene_bump_dry_run.get_genes", get_genes), patch(
        "schema_bump_dry_run_scripts.genes.gene_bump_dry_run.logger"
    ) as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            sample_dataset, sample_diff_map, sample_deprecated_datasets
        )

    mock_logger.info.assert_called_once_with(f"Dataset {sample_dataset['dataset_id']} has 2 deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_new_collection_and_deprecated_genes(
    sample_diff_map, sample_deprecated_datasets, sample_dataset
):
    sample_dataset.update(
        **{
            "collection_id": "new_collection",
            "get_genes_mock_response": ["deprecated:1-2", "deprecated:1-3", "current:1-2", "current:1-3"],
        }
    )
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    expected_deprecated_datasets["new_collection"] = {
        "num_datasets": 1,
        "deprecated_genes": {"deprecated:1-2", "deprecated:1-3"},
        "genes": {"deprecated:1-2", "deprecated:1-3", "current:1-2", "current:1-3"},
    }

    expected_is_deprecated_genes_found = True
    with patch("schema_bump_dry_run_scripts.genes.gene_bump_dry_run.get_genes", get_genes), patch(
        "schema_bump_dry_run_scripts.genes.gene_bump_dry_run.logger"
    ) as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            sample_dataset, sample_diff_map, sample_deprecated_datasets
        )

    mock_logger.info.assert_called_once_with(f"Dataset {sample_dataset['dataset_id']} has 2 deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_existing_collection_multiple_organisms_and_deprecated_genes(
    sample_dataset, sample_deprecated_datasets, sample_diff_map
):
    sample_dataset.update(
        **{
            "collection_id": "existing_collection",
            "organism": [{"ontology_term_id": "organism-1"}, {"ontology_term_id": "organism-2"}],
            "get_genes_mock_response": ["deprecated:1-3", "deprecated:2-1", "current:1-3", "current:2-1"],
        }
    )
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    expected_deprecated_datasets["existing_collection"] = {
        "num_datasets": 3,
        "deprecated_genes": {"deprecated:1-1", "deprecated:1-2", "deprecated:1-3", "deprecated:2-1"},
        "genes": {
            "current:1-1",
            "current:1-2",
            "deprecated:1-1",
            "deprecated:1-2",
            "deprecated:1-3",
            "deprecated:2-1",
            "current:1-3",
            "current:2-1",
        },
    }
    expected_is_deprecated_genes_found = True
    with patch("schema_bump_dry_run_scripts.genes.gene_bump_dry_run.get_genes", get_genes), patch(
        "schema_bump_dry_run_scripts.genes.gene_bump_dry_run.logger"
    ) as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            sample_dataset, sample_diff_map, sample_deprecated_datasets
        )

    mock_logger.info.assert_called_once_with(f"Dataset {sample_dataset['dataset_id']} has 2 deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_new_collection_multiple_organisms_and_deprecated_genes(
    sample_diff_map, sample_deprecated_datasets, sample_dataset
):
    sample_dataset.update(
        **{
            "collection_id": "new_collection",
            "organism": [{"ontology_term_id": "organism-1"}, {"ontology_term_id": "organism-2"}],
            "get_genes_mock_response": ["deprecated:1-1", "current:1-1", "deprecated:2-1", "current:2-1"],
        }
    )
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    expected_deprecated_datasets["new_collection"] = {
        "num_datasets": 1,
        "deprecated_genes": {"deprecated:1-1", "deprecated:2-1"},
        "genes": {"deprecated:1-1", "deprecated:2-1", "current:1-1", "current:2-1"},
    }

    expected_is_deprecated_genes_found = True
    with patch("schema_bump_dry_run_scripts.genes.gene_bump_dry_run.get_genes", get_genes), patch(
        "schema_bump_dry_run_scripts.genes.gene_bump_dry_run.logger"
    ) as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            sample_dataset, sample_diff_map, sample_deprecated_datasets
        )

    mock_logger.info.assert_called_once_with(f"Dataset {sample_dataset['dataset_id']} has 2 deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_generate_deprecated_private(sample_diff_map):
    # Mock data
    base_url = "https://example.com"
    dataset3 = {"dataset_id": "3", "collection_id": "new_collection", "revision_of": "collection_public"}
    dataset4 = {"dataset_id": "4", "collection_id": "collection_with_revision", "revision_of": None}

    # Expected output
    expected_deprecated_datasets = {
        "new_collection": {"revision_of": "collection_public"},
        "collection_with_revision": {},
    }
    expected_non_auto_migrated = ["collection_public"]

    fetch_private_datasets_mock = Mock(return_value=[(dataset3, "collection_public"), (dataset4, None)])
    compare_genes_mock = Mock(
        side_effect=[
            ({"new_collection": {}}, True),
            ({"new_collection": {"revision_of": "collection_public"}, "collection_with_revision": {}}, True),
        ]
    )
    with patch(
        "schema_bump_dry_run_scripts.genes.gene_bump_dry_run.fetch_private_datasets", fetch_private_datasets_mock
    ), patch("schema_bump_dry_run_scripts.genes.gene_bump_dry_run.compare_genes", compare_genes_mock):
        deprecated_datasets, non_auto_migrated = generate_deprecated_private(base_url, sample_diff_map)

    # Assert the results
    assert deprecated_datasets == expected_deprecated_datasets
    assert non_auto_migrated == expected_non_auto_migrated
    fetch_private_datasets_mock.assert_called_once_with(base_url)


def test_generate_deprecated_private__with_datasets_in_same_collectionn(sample_diff_map):
    # Mock data
    base_url = "https://example.com"
    dataset1 = {"dataset_id": "1", "collection_id": "existing_collection", "revision_of": None}
    dataset2 = {"dataset_id": "2", "collection_id": "existing_collection"}

    # Expected output
    expected_deprecated_datasets = {"existing_collection": {}}
    expected_non_auto_migrated = []

    fetch_private_datasets_mock = Mock(return_value=[(dataset1, None), (dataset2, None)])
    compare_genes_mock = Mock(side_effect=[({"existing_collection": {}}, False), ({"existing_collection": {}}, True)])
    with patch(
        "schema_bump_dry_run_scripts.genes.gene_bump_dry_run.fetch_private_datasets", fetch_private_datasets_mock
    ), patch("schema_bump_dry_run_scripts.genes.gene_bump_dry_run.compare_genes", compare_genes_mock):
        deprecated_datasets, non_auto_migrated = generate_deprecated_private(base_url, sample_diff_map)

    # Assert the results
    assert deprecated_datasets == expected_deprecated_datasets
    assert non_auto_migrated == expected_non_auto_migrated
    fetch_private_datasets_mock.assert_called_once_with(base_url)


def test_generate_deprecated_public__with_datasets_in_same_collectionn(sample_diff_map):
    # Mock data
    base_url = "https://example.com"
    dataset1 = {"dataset_id": "1", "collection_id": "existing_collection"}
    dataset2 = {"dataset_id": "2", "collection_id": "existing_collection"}

    # Expected output
    expected_deprecated_datasets = {"existing_collection": {}}

    fetch_public_datasets_mock = Mock(return_value=[(dataset1, None), (dataset2, None)])
    compare_genes_mock = Mock(side_effect=[({"existing_collection": {}}, False), ({"existing_collection": {}}, True)])
    with patch(
        "schema_bump_dry_run_scripts.genes.gene_bump_dry_run.fetch_public_datasets", fetch_public_datasets_mock
    ), patch("schema_bump_dry_run_scripts.genes.gene_bump_dry_run.compare_genes", compare_genes_mock):
        deprecated_datasets = generate_deprecated_public(base_url, sample_diff_map)

    # Assert the results
    assert deprecated_datasets == expected_deprecated_datasets
    fetch_public_datasets_mock.assert_called_once_with(base_url)


def test_generate_deprecated_public(sample_diff_map):
    # Mock data
    base_url = "https://example.com"
    dataset3 = {
        "dataset_id": "dataset3",
        "collection_id": "new_collection",
    }
    dataset4 = {
        "dataset_id": "dataset4",
        "collection_id": "existing_collection",
    }

    # Expected output
    expected_deprecated_datasets = {"new_collection": {}, "existing_collection": {}}

    fetch_public_datasets_mock = Mock(return_value=[(dataset3, "collection_public"), (dataset4, None)])
    compare_genes_mock = Mock(
        side_effect=[
            ({"new_collection": {}}, True),
            ({"new_collection": {}, "existing_collection": {}}, True),
        ]
    )
    with patch(
        "schema_bump_dry_run_scripts.genes.gene_bump_dry_run.fetch_public_datasets", fetch_public_datasets_mock
    ), patch("schema_bump_dry_run_scripts.genes.gene_bump_dry_run.compare_genes", compare_genes_mock):
        deprecated_datasets = generate_deprecated_public(base_url, sample_diff_map)

    # Assert the results
    assert deprecated_datasets == expected_deprecated_datasets
    fetch_public_datasets_mock.assert_called_once_with(base_url)
