from unittest.mock import Mock, patch

import pytest

from cellxgene_schema_cli.cellxgene_schema.gencode import SupportedOrganisms
from scripts.schema_bump_dry_run_genes.gene_bump_dry_run import (
    compare_genes,
    generate_deprecated_private,
    generate_deprecated_public,
    generate_report,
    get_diff_map,
)


def test_get_diff_map(tmp_path):  # type: ignore
    for key in SupportedOrganisms:
        with open(f"{tmp_path}/{key.name}_diff.txt", "w") as fp:
            fp.write("test")
    with patch("scripts.schema_bump_dry_run_genes.gene_bump_dry_run.GENCODE_DIR", tmp_path):
        diff_map = get_diff_map()
    
    # one diff-map for each of our supported species
    assert len(diff_map) == len(SupportedOrganisms)

    # each species should have a diff map
    assert len(set([x.value for x in SupportedOrganisms]).difference([k for k in diff_map.keys()])) == 0


@pytest.fixture
def sample_report_data():  # type: ignore
    return {
        "deprecated_public": {
            "with multiple dataset groups": {
                "dataset_groups": [
                    {
                        "datasets": ["dataset_1", "dataset_2"],
                        "num_datasets": 2,
                        "num_genes": 5,
                        "deprecated_genes": ["term_1", "term_2", "term_3"],
                    },
                    {
                        "datasets": ["dataset_3"],
                        "num_datasets": 1,
                        "num_genes": 4,
                        "deprecated_genes": ["term_4", "term_5"],
                    },
                ]
            },
            "with single dataset group": {
                "dataset_groups": [
                    {
                        "datasets": [
                            "dataset_4",
                            "dataset_5",
                            "dataset_4",
                            "dataset_5",
                            "dataset_4",
                            "dataset_5",
                            "dataset_4",
                            "dataset_5",
                            "dataset_4",
                            "dataset_5",
                            "dataset_4",
                            "dataset_5",
                            "dataset_4",
                            "dataset_5",
                            "dataset_4",
                            "dataset_5",
                        ],
                        "num_datasets": 16,
                        "num_genes": 50,
                        "deprecated_genes": [
                            "term_6",
                            "term_6",
                            "term_6",
                            "term_6",
                            "term_6",
                            "term_6",
                            "term_6",
                            "term_6",
                            "term_6",
                            "term_6",
                            "term_6",
                            "term_6",
                            "term_6",
                            "term_6",
                        ],
                    }
                ]
            },
        },
        "open_revisions": {
            "with revision": {
                "revision_of": "collection_id_1",
                "dataset_groups": [
                    {
                        "datasets": ["dataset_6"],
                        "num_datasets": 1,
                        "num_genes": 6,
                        "deprecated_genes": ["term_7", "term_8", "term_9", "term_10"],
                    }
                ],
            },
            "no revision": {
                "dataset_groups": [
                    {
                        "datasets": ["dataset_7", "dataset_8", "dataset_9"],
                        "num_datasets": 3,
                        "num_genes": 2,
                        "deprecated_genes": ["term1"],
                    }
                ]
            },
        },
        "non_auto_migrated": ["collection_id_2", "collection_id_4"],
    }


def test_generate_report(sample_report_data):  # type: ignore
    expected_report = """## Deprecated Features in Public Collections:

Collection ID: with multiple dataset groups
    Affected Datasets:
        dataset_1, dataset_2
    Number of Affected Datasets: 2
    Number of Deprecated Features in Group: 3
    Number of Features in Group: 5
    Deprecated Features:
        term_1, term_2, term_3

    Affected Datasets:
        dataset_3
    Number of Affected Datasets: 1
    Number of Deprecated Features in Group: 2
    Number of Features in Group: 4
    Deprecated Features:
        term_4, term_5

Collection ID: with single dataset group
    Affected Datasets:
        dataset_4, dataset_5, dataset_4, dataset_5, dataset_4, dataset_5, dataset_4,
        dataset_5, dataset_4, dataset_5, dataset_4, dataset_5, dataset_4, dataset_5,
        dataset_4, dataset_5
    Number of Affected Datasets: 16
    Number of Deprecated Features in Group: 14
    Number of Features in Group: 50
    Deprecated Features:
        term_6, term_6, term_6, term_6, term_6, term_6, term_6, term_6, term_6,
        term_6, term_6, term_6, term_6, term_6

## Deprecated Features in Private Collections:

Collection ID: with revision
Note--In A Revision of: collection_id_1
    Affected Datasets:
        dataset_6
    Number of Affected Datasets: 1
    Number of Deprecated Features in Group: 4
    Number of Features in Group: 6
    Deprecated Features:
        term_7, term_8, term_9, term_10

Collection ID: no revision
    Affected Datasets:
        dataset_7, dataset_8, dataset_9
    Number of Affected Datasets: 3
    Number of Deprecated Features in Group: 1
    Number of Features in Group: 2
    Deprecated Features:
        term1

## The Following Public Collections Will Not Be Auto-Migrated Due To Having an Open Revision:
collection_id_2
collection_id_4
"""
    # Assert results
    assert generate_report(sample_report_data) == expected_report


def test_generate_report__with_empty_data():  # type: ignore
    empty_data = {"deprecated_public": {}, "open_revisions": {}, "non_auto_migrated": []}  # type: ignore
    expected_report = """## Deprecated Features in Public Collections:

## Deprecated Features in Private Collections:

## The Following Public Collections Will Not Be Auto-Migrated Due To Having an Open Revision:
"""
    assert generate_report(empty_data) == expected_report


def get_genes(dataset):  # type: ignore
    # Mock implementation for testing purposes
    return dataset.get("get_genes_mock_response", [])


@pytest.fixture
def sample_diff_map():  # type: ignore
    return {
        "organism-1": ["deprecated:1-1", "deprecated:1-2", "deprecated:1-3"],
        "organism-2": ["deprecated:2-1", "deprecated:2-2"],
    }


@pytest.fixture
def sample_deprecated_datasets():  # type: ignore
    group = {
        "datasets": ["dataset1"],
        "num_datasets": 1,
        "num_deprecated_genes": 2,
        "num_genes": 4,
        "deprecated_genes": {"deprecated:1-1", "deprecated:1-2"},
    }
    return {
        "existing_collection": {"dataset_groups": {(*sorted(group["deprecated_genes"]), group["num_genes"]): group}}  # type: ignore
    }


@pytest.fixture
def sample_dataset():  # type: ignore
    return {
        "dataset_id": "dataset1",
        "collection_id": "new_collection",
        "organism": [{"ontology_term_id": "organism-1"}],
    }


def test_compare_genes__with_no_deprecated_genes(sample_diff_map, sample_deprecated_datasets, sample_dataset):  # type: ignore
    # Expected results
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    expected_is_deprecated_genes_found = False

    with patch("scripts.schema_bump_dry_run_genes.gene_bump_dry_run.get_genes", get_genes), patch(
        "scripts.schema_bump_dry_run_genes.gene_bump_dry_run.logger"
    ) as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            sample_dataset, sample_diff_map, sample_deprecated_datasets
        )

    # Assert the results
    mock_logger.debug.assert_called_once_with("Dataset dataset1 has no deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_empty_diff_map(sample_deprecated_datasets, sample_dataset):  # type: ignore
    # Expected results
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    expected_is_deprecated_genes_found = False

    with patch("scripts.schema_bump_dry_run_genes.gene_bump_dry_run.get_genes", get_genes), patch(
        "scripts.schema_bump_dry_run_genes.gene_bump_dry_run.logger"
    ) as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            sample_dataset, {}, sample_deprecated_datasets
        )

    # Assert the results
    mock_logger.debug.assert_called_once_with(f"Dataset {sample_dataset['dataset_id']} has no deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_existing_collection_deprecated_genes_and_new_dataset_group(  # type: ignore
    sample_diff_map, sample_deprecated_datasets, sample_dataset
):
    # Mock data
    sample_dataset.update(
        **{
            "dataset_id": "dataset2",
            "collection_id": "existing_collection",
            "get_genes_mock_response": ["deprecated:1-2", "deprecated:1-3", "current:1-2", "current:1-3"],
        }
    )

    # Expected results
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    dataset_group_key = ("deprecated:1-2", "deprecated:1-3", 4)
    expected_deprecated_datasets["existing_collection"]["dataset_groups"][dataset_group_key] = {
        "datasets": ["dataset2"],
        "deprecated_genes": {"deprecated:1-2", "deprecated:1-3"},
        "num_datasets": 1,
        "num_genes": 4,
    }
    expected_is_deprecated_genes_found = True

    with patch("scripts.schema_bump_dry_run_genes.gene_bump_dry_run.get_genes", get_genes), patch(
        "scripts.schema_bump_dry_run_genes.gene_bump_dry_run.logger"
    ) as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            sample_dataset, sample_diff_map, sample_deprecated_datasets
        )

    # Assert results
    mock_logger.debug.assert_called_once_with(f"Dataset {sample_dataset['dataset_id']} has 2 deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_existing_collection_deprecated_genes_and_existing_dataset_group(  # type: ignore
    sample_diff_map, sample_deprecated_datasets, sample_dataset
):
    # Mock data
    sample_dataset.update(
        **{
            "dataset_id": "dataset2",
            "collection_id": "existing_collection",
            "get_genes_mock_response": ["deprecated:1-2", "deprecated:1-1", "current:1-2", "current:1-3"],
        }
    )

    # Expected output
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    dataset_group_key = ("deprecated:1-1", "deprecated:1-2", 4)
    expected_deprecated_datasets["existing_collection"]["dataset_groups"][dataset_group_key]["num_datasets"] += 1
    expected_deprecated_datasets["existing_collection"]["dataset_groups"][dataset_group_key]["datasets"].append(
        "dataset2"
    )
    expected_is_deprecated_genes_found = True

    with patch("scripts.schema_bump_dry_run_genes.gene_bump_dry_run.get_genes", get_genes), patch(
        "scripts.schema_bump_dry_run_genes.gene_bump_dry_run.logger"
    ) as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            sample_dataset, sample_diff_map, sample_deprecated_datasets
        )

    # Assert the results
    mock_logger.debug.assert_called_once_with(f"Dataset {sample_dataset['dataset_id']} has 2 deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_new_collection_and_deprecated_genes(sample_diff_map, sample_dataset):  # type: ignore
    # Mock data
    sample_dataset.update(
        **{
            "collection_id": "new_collection",
            "get_genes_mock_response": ["deprecated:1-2", "deprecated:1-3", "current:1-2", "current:1-3"],
        }
    )

    # Expected output
    dataset_group_key = ("deprecated:1-2", "deprecated:1-3", 4)
    expected_deprecated_datasets = {
        "new_collection": {
            "dataset_groups": {
                dataset_group_key: {
                    "datasets": ["dataset1"],
                    "num_datasets": 1,
                    "deprecated_genes": {"deprecated:1-2", "deprecated:1-3"},
                    "num_genes": 4,
                }
            }
        }
    }
    expected_is_deprecated_genes_found = True
    with patch("scripts.schema_bump_dry_run_genes.gene_bump_dry_run.get_genes", get_genes), patch(
        "scripts.schema_bump_dry_run_genes.gene_bump_dry_run.logger"
    ) as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(sample_dataset, sample_diff_map, {})

    # Assert the results
    mock_logger.debug.assert_called_once_with(f"Dataset {sample_dataset['dataset_id']} has 2 deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_existing_collection_multiple_organisms_and_deprecated_genes(  # type: ignore
    sample_dataset, sample_deprecated_datasets, sample_diff_map
):
    # Mock data
    sample_dataset.update(
        **{
            "collection_id": "existing_collection",
            "dataset_id": "dataset2",
            "organism": [{"ontology_term_id": "organism-1"}, {"ontology_term_id": "organism-2"}],
            "get_genes_mock_response": ["deprecated:1-3", "deprecated:2-1", "current:1-3", "current:2-1"],
        }
    )

    # Excpected output
    expected_deprecated_datasets = sample_deprecated_datasets.copy()
    dataset_group_key = ("deprecated:1-3", "deprecated:2-1", 4)
    expected_deprecated_datasets["existing_collection"][dataset_group_key] = {
        "datasets": ["dataset2"],
        "num_datasets": 1,
        "deprecated_genes": {"deprecated:1-3", "deprecated:2-1"},
        "num_genes": 4,
    }
    expected_is_deprecated_genes_found = True
    with patch("scripts.schema_bump_dry_run_genes.gene_bump_dry_run.get_genes", get_genes), patch(
        "scripts.schema_bump_dry_run_genes.gene_bump_dry_run.logger"
    ) as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(
            sample_dataset, sample_diff_map, sample_deprecated_datasets
        )

    # Assert the results
    mock_logger.debug.assert_called_once_with(f"Dataset {sample_dataset['dataset_id']} has 2 deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_compare_genes__with_new_collection_multiple_organisms_and_deprecated_genes(sample_diff_map, sample_dataset):  # type: ignore
    # Mock data
    sample_dataset.update(
        **{
            "collection_id": "new_collection",
            "dataset_id": "dataset2",
            "organism": [{"ontology_term_id": "organism-1"}, {"ontology_term_id": "organism-2"}],
            "get_genes_mock_response": ["deprecated:1-1", "current:1-1", "deprecated:2-1", "current:2-1"],
        }
    )
    dataset_group_key = ("deprecated:1-1", "deprecated:2-1", 4)

    # Expected output
    expected_deprecated_datasets = {
        "new_collection": {
            "dataset_groups": {
                dataset_group_key: {
                    "datasets": ["dataset2"],
                    "num_datasets": 1,
                    "deprecated_genes": {"deprecated:1-1", "deprecated:2-1"},
                    "num_genes": 4,
                }
            }
        }
    }

    expected_is_deprecated_genes_found = True
    with patch("scripts.schema_bump_dry_run_genes.gene_bump_dry_run.get_genes", get_genes), patch(
        "scripts.schema_bump_dry_run_genes.gene_bump_dry_run.logger"
    ) as mock_logger:
        actual_deprecated_datasets, is_deprecated_genes_found = compare_genes(sample_dataset, sample_diff_map, {})

    # Assert the results
    mock_logger.debug.assert_called_once_with(f"Dataset {sample_dataset['dataset_id']} has 2 deprecated genes")
    assert actual_deprecated_datasets == expected_deprecated_datasets
    assert is_deprecated_genes_found == expected_is_deprecated_genes_found


def test_generate_deprecated_private__with_revision_collection(sample_diff_map):  # type: ignore
    # Mock data
    base_url = "https://example.com"
    dataset3 = {"dataset_id": "3", "collection_id": "collection_with_revision", "revision_of": "collection_public"}

    fetch_private_datasets_mock = Mock(return_value=[(dataset3, "collection_public")])
    compare_genes_mock = Mock(
        side_effect=[
            ({"collection_with_revision": {"dataset_groups": {"a": "test"}, "revision_of": "collection_public"}}, True),
        ]
    )

    # Expected output
    expected_deprecated_datasets = {
        "collection_with_revision": {"dataset_groups": ["test"], "revision_of": "collection_public"},
    }
    expected_non_auto_migrated = {"collection_public"}

    with patch(
        "scripts.schema_bump_dry_run_genes.gene_bump_dry_run.fetch_private_datasets", fetch_private_datasets_mock
    ), patch("scripts.schema_bump_dry_run_genes.gene_bump_dry_run.compare_genes", compare_genes_mock):
        deprecated_datasets, non_auto_migrated = generate_deprecated_private(base_url, sample_diff_map)

    # Assert the results
    assert deprecated_datasets == expected_deprecated_datasets
    assert non_auto_migrated == expected_non_auto_migrated  # type: ignore
    fetch_private_datasets_mock.assert_called_once_with(base_url)


def test_generate_deprecated_private__with_datasets_in_same_collectionn(sample_diff_map):  # type: ignore
    # Mock data
    base_url = "https://example.com"
    dataset1 = {"dataset_id": "1", "collection_id": "existing_collection"}
    dataset2 = {"dataset_id": "2", "collection_id": "existing_collection"}

    fetch_private_datasets_mock = Mock(return_value=[(dataset1, None), (dataset2, None)])
    compare_genes_mock = Mock(
        side_effect=[({}, False), ({"existing_collection": {"dataset_groups": {"a": "test", "b": "test2"}}}, True)]
    )

    # Expected output
    expected_deprecated_datasets = {"existing_collection": {"dataset_groups": ["test", "test2"]}}
    expected_non_auto_migrated = set()  # type: ignore

    with patch(
        "scripts.schema_bump_dry_run_genes.gene_bump_dry_run.fetch_private_datasets", fetch_private_datasets_mock
    ), patch("scripts.schema_bump_dry_run_genes.gene_bump_dry_run.compare_genes", compare_genes_mock):
        deprecated_datasets, non_auto_migrated = generate_deprecated_private(base_url, sample_diff_map)

    # Assert the results
    assert deprecated_datasets == expected_deprecated_datasets
    assert non_auto_migrated == expected_non_auto_migrated
    fetch_private_datasets_mock.assert_called_once_with(base_url)


def test_generate_deprecated_public__with_datasets_in_same_collectionn(sample_diff_map):  # type: ignore
    # Mock data
    base_url = "https://example.com"
    dataset1 = {"dataset_id": "1", "collection_id": "existing_collection"}
    dataset2 = {"dataset_id": "2", "collection_id": "existing_collection"}

    fetch_public_datasets_mock = Mock(return_value=[(dataset1, None), (dataset2, None)])
    compare_genes_mock = Mock(
        side_effect=[({}, False), ({"existing_collection": {"dataset_groups": {"a": "test", "b": "test2"}}}, True)]
    )

    # Expected output
    expected_deprecated_datasets = {"existing_collection": {"dataset_groups": ["test", "test2"]}}

    with patch(
        "scripts.schema_bump_dry_run_genes.gene_bump_dry_run.fetch_public_datasets", fetch_public_datasets_mock
    ), patch("scripts.schema_bump_dry_run_genes.gene_bump_dry_run.compare_genes", compare_genes_mock):
        deprecated_datasets = generate_deprecated_public(base_url, sample_diff_map)

    # Assert the results
    assert deprecated_datasets == expected_deprecated_datasets
    fetch_public_datasets_mock.assert_called_once_with(base_url)


def test_generate_deprecated_public(sample_diff_map):  # type: ignore
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
    expected_deprecated_datasets = {
        "new_collection": {"dataset_groups": ["test1"]},
        "existing_collection": {"dataset_groups": ["test2"]},
    }

    fetch_public_datasets_mock = Mock(return_value=[(dataset3, "collection_public"), (dataset4, None)])
    compare_genes_mock = Mock(
        side_effect=[
            ({}, False),
            (
                {
                    "new_collection": {"dataset_groups": {"a": "test1"}},
                    "existing_collection": {"dataset_groups": {"b": "test2"}},
                },
                False,
            ),
        ]
    )
    with patch(
        "scripts.schema_bump_dry_run_genes.gene_bump_dry_run.fetch_public_datasets", fetch_public_datasets_mock
    ), patch("scripts.schema_bump_dry_run_genes.gene_bump_dry_run.compare_genes", compare_genes_mock):
        deprecated_datasets = generate_deprecated_public(base_url, sample_diff_map)

    # Assert the results
    assert deprecated_datasets == expected_deprecated_datasets
    fetch_public_datasets_mock.assert_called_once_with(base_url)
