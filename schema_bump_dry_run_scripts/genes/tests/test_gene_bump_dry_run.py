from unittest.mock import Mock, patch

import pytest

from schema_bump_dry_run_scripts.genes.gene_bump_dry_run import (
    compare_genes,
    generate_deprecated_private,
    generate_deprecated_public,
    generate_report,
    get_diff_map,
)


def test_get_diff_map():
    diff_map = get_diff_map()
    assert len(diff_map) == 4
    for key in diff_map:
        assert key in ["NCBITaxon:9606", "NCBITaxon:10090", "NCBITaxon:2697049", "NCBITaxon:32630"]


@pytest.fixture
def sample_report_data():
    return {
        "deprecated_public": {
            "collection1": {
                "num_datasets": 3,
                "deprecated_terms": [
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
            "collection2": {"num_datasets": 1, "deprecated_terms": ["term4"]},
            "collection4": {"num_datasets": 2, "deprecated_terms": ["gene1", "gene2"]},
        },
        "open_revisions": {
            "collection3": {"revision_of": "collection4", "num_datasets": 2, "deprecated_terms": ["gene1", "gene2"]},
            "collection5": {"num_datasets": 1, "deprecated_terms": ["gene3"]},
        },
        "non_auto_migrated": ["collection6", "collection7"],
    }


def test_generate_report(sample_report_data):
    expected_report = """## Deprecated Terms in Public Datasets:

Collection ID: collection1
Number of Affected Datasets: 3
Deprecated Terms: 
    term1, term2, term3, term1, term2, term3, term1, term2, term3, term1, term2,
    term3, term1, term2, term3, term1, term2, term3

Collection ID: collection2
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
        "organism1": ["gene1-1", "gene1-2", "gene1-3"],
        "organism2": ["gene2-1", "gene2-2"],
    }


@pytest.fixture
def sample_deprecated_datasets():
    return {"collection1": {"num_datasets": 2, "deprecated_terms": {"gene1-1", "gene1-2"}}}


def test_compare_genes__with_no_deprecated_genes(sample_diff_map, sample_deprecated_datasets):
    dataset = {
        "dataset_id": "dataset1",
        "collection_id": "collection2",
        "organism": [{"ontology_term_id": "organism2"}],
        "get_genes_mock_response": ["geneA", "geneB", "geneC"],
    }
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    expected_is_deprecated_genes_found = False

    with patch("gene_bump_dry_run.get_genes", get_genes), patch("gene_bump_dry_run.logger") as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            dataset, sample_diff_map, sample_deprecated_datasets
        )

    mock_logger.info.assert_called_once_with("Dataset dataset1 has no deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_empty_diff_map(sample_deprecated_datasets):
    dataset = {
        "dataset_id": "dataset3",
        "collection_id": "collection3",
        "organism": [{"ontology_term_id": "organism3"}],
        "get_genes_mock_response": ["gene3-1", "gene3-2", "gene3-3"],
    }
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    expected_is_deprecated_genes_found = False

    with patch("gene_bump_dry_run.get_genes", get_genes), patch("gene_bump_dry_run.logger") as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(dataset, {}, sample_deprecated_datasets)

    mock_logger.info.assert_called_once_with("Dataset dataset3 has no deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_existing_collection_and_deprecated_genes(sample_diff_map, sample_deprecated_datasets):
    dataset = {
        "dataset_id": "dataset2",
        "collection_id": "collection1",
        "organism": [{"ontology_term_id": "organism1"}],
        "get_genes_mock_response": ["gene1-1", "gene1-2", "gene1-3", "gene1-10"],
    }
    expected_deprecated_datasets = {
        "collection1": {"num_datasets": 3, "deprecated_terms": {"gene1-1", "gene1-2", "gene1-3"}}
    }
    expected_is_deprecated_genes_found = True
    with patch("gene_bump_dry_run.get_genes", get_genes), patch("gene_bump_dry_run.logger") as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            dataset, sample_diff_map, sample_deprecated_datasets
        )

    mock_logger.info.assert_called_once_with("Dataset dataset2 has 3 deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_new_collection_and_deprecated_genes(sample_diff_map, sample_deprecated_datasets):
    dataset = {
        "dataset_id": "dataset4",
        "collection_id": "collection2",
        "organism": [{"ontology_term_id": "organism1"}],
        "get_genes_mock_response": ["gene1-1", "gene1-2", "gene1-3", "gene1-10"],
    }
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    expected_deprecated_datasets.update(
        collection2={"num_datasets": 1, "deprecated_terms": {"gene1-1", "gene1-2", "gene1-3"}}
    )
    expected_is_deprecated_genes_found = True
    with patch("gene_bump_dry_run.get_genes", get_genes), patch("gene_bump_dry_run.logger") as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            dataset, sample_diff_map, sample_deprecated_datasets
        )

    mock_logger.info.assert_called_once_with("Dataset dataset4 has 3 deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_multiple_organisms_and_deprecated_genes(sample_diff_map, sample_deprecated_datasets):
    dataset = {
        "dataset_id": "dataset2",
        "collection_id": "collection1",
        "organism": [{"ontology_term_id": "organism1"}, {"ontology_term_id": "organism2"}],
        "get_genes_mock_response": ["gene1-1", "gene2-2", "geneA"],
    }
    expected_deprecated_datasets = {
        "collection1": {"num_datasets": 3, "deprecated_terms": {"gene1-1", "gene1-2", "gene2-2"}}
    }
    expected_is_deprecated_genes_found = True
    with patch("gene_bump_dry_run.get_genes", get_genes), patch("gene_bump_dry_run.logger") as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            dataset, sample_diff_map, sample_deprecated_datasets
        )

    mock_logger.info.assert_called_once_with("Dataset dataset2 has 2 deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_new_collection_multiple_organisms_and_deprecated_genes(
    sample_diff_map, sample_deprecated_datasets
):
    dataset = {
        "dataset_id": "dataset2",
        "collection_id": "collection2",
        "organism": [{"ontology_term_id": "organism1"}, {"ontology_term_id": "organism2"}],
        "get_genes_mock_response": ["gene1-1", "gene2-2", "geneA"],
    }
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    expected_deprecated_datasets.update(collection2={"num_datasets": 1, "deprecated_terms": {"gene1-1", "gene2-2"}})
    expected_is_deprecated_genes_found = True
    with patch("gene_bump_dry_run.get_genes", get_genes), patch("gene_bump_dry_run.logger") as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            dataset, sample_diff_map, sample_deprecated_datasets
        )

    mock_logger.info.assert_called_once_with("Dataset dataset2 has 2 deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_generate_deprecated_private(sample_diff_map):
    # Mock data
    base_url = "https://example.com"
    dataset3 = {"dataset_id": "3", "collection_id": "collection2", "revision_of": "collection_public"}
    dataset4 = {"dataset_id": "4", "collection_id": "collection3", "revision_of": None}

    # Expected output
    expected_deprecated_datasets = {"collection2": {"revision_of": "collection_public"}, "collection3": {}}
    expected_non_auto_migrated = ["collection_public"]

    fetch_private_datasets_mock = Mock(return_value=[(dataset3, "collection_public"), (dataset4, None)])
    compare_genes_mock = Mock(
        side_effect=[
            ({"collection2": {}}, True),
            ({"collection2": {"revision_of": "collection_public"}, "collection3": {}}, True),
        ]
    )
    with patch("gene_bump_dry_run.fetch_private_datasets", fetch_private_datasets_mock), patch(
        "gene_bump_dry_run.compare_genes", compare_genes_mock
    ):
        deprecated_datasets, non_auto_migrated = generate_deprecated_private(base_url, sample_diff_map)

    # Assert the results
    assert deprecated_datasets == expected_deprecated_datasets
    assert non_auto_migrated == expected_non_auto_migrated
    fetch_private_datasets_mock.assert_called_once_with(base_url)


def test_generate_deprecated_private__with_datasets_in_same_collectionn(sample_diff_map):
    # Mock data
    base_url = "https://example.com"
    dataset1 = {"dataset_id": "1", "collection_id": "collection1", "revision_of": None}
    dataset2 = {"dataset_id": "2", "collection_id": "collection1", "revision_of": None}

    # Expected output
    expected_deprecated_datasets = {"collection1": {}}
    expected_non_auto_migrated = []

    fetch_private_datasets_mock = Mock(return_value=[(dataset1, None), (dataset2, None)])
    compare_genes_mock = Mock(side_effect=[({"collection1": {}}, False), ({"collection1": {}}, True)])
    with patch("gene_bump_dry_run.fetch_private_datasets", fetch_private_datasets_mock), patch(
        "gene_bump_dry_run.compare_genes", compare_genes_mock
    ):
        deprecated_datasets, non_auto_migrated = generate_deprecated_private(base_url, sample_diff_map)

    # Assert the results
    assert deprecated_datasets == expected_deprecated_datasets
    assert non_auto_migrated == expected_non_auto_migrated
    fetch_private_datasets_mock.assert_called_once_with(base_url)


def test_generate_deprecated_public__with_datasets_in_same_collectionn(sample_diff_map):
    # Mock data
    base_url = "https://example.com"
    dataset1 = {"dataset_id": "1", "collection_id": "collection1"}
    dataset2 = {"dataset_id": "2", "collection_id": "collection1"}

    # Expected output
    expected_deprecated_datasets = {"collection1": {}}

    fetch_public_datasets_mock = Mock(return_value=[(dataset1, None), (dataset2, None)])
    compare_genes_mock = Mock(side_effect=[({"collection1": {}}, False), ({"collection1": {}}, True)])
    with patch("gene_bump_dry_run.fetch_public_datasets", fetch_public_datasets_mock), patch(
        "gene_bump_dry_run.compare_genes", compare_genes_mock
    ):
        deprecated_datasets = generate_deprecated_public(base_url, sample_diff_map)

    # Assert the results
    assert deprecated_datasets == expected_deprecated_datasets
    fetch_public_datasets_mock.assert_called_once_with(base_url)


def test_generate_deprecated_public(sample_diff_map):
    # Mock data
    base_url = "https://example.com"
    dataset3 = {
        "dataset_id": "3",
        "collection_id": "collection2",
    }
    dataset4 = {
        "dataset_id": "4",
        "collection_id": "collection3",
    }

    # Expected output
    expected_deprecated_datasets = {"collection2": {}, "collection3": {}}

    fetch_public_datasets_mock = Mock(return_value=[(dataset3, "collection_public"), (dataset4, None)])
    compare_genes_mock = Mock(
        side_effect=[
            ({"collection2": {}}, True),
            ({"collection2": {}, "collection3": {}}, True),
        ]
    )
    with patch("gene_bump_dry_run.fetch_public_datasets", fetch_public_datasets_mock), patch(
        "gene_bump_dry_run.compare_genes", compare_genes_mock
    ):
        deprecated_datasets = generate_deprecated_public(base_url, sample_diff_map)

    # Assert the results
    assert deprecated_datasets == expected_deprecated_datasets
    fetch_public_datasets_mock.assert_called_once_with(base_url)
