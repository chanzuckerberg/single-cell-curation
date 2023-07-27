import os
from unittest import mock

import pytest
from jinja2 import Template

from scripts.migration_assistant.generate_script import generate_script, get_deprecated_feature_ids, get_template


# Define a fixture for the template
@pytest.fixture
def template():
    return Template("{{ schema_version }}")


# Define the test cases
def test_generate_script__without_gencode_changes(template, tmpdir):
    # Prepare the input data
    ontology_term_map = {
        "assay": {},
        "cell_type": {},
        "development_stage": {},
        "disease": {},
        "organism": {},
        "self_reported_ethnicity": {},
        "sex": {},
        "tissue": {},
    }
    gencode_term_map = []
    mock_target_file = tmpdir + "/migrate.py"
    with mock.patch("scripts.migration_assistant.generate_script.target_file", mock_target_file):
        # Execute the function
        generate_script(get_template(), ontology_term_map, gencode_term_map)

    # Verify the output
    expected_output = """
import anndata as ad

from . import utils

# fmt: off
# ONTOLOGY TERMS TO UPDATE ACROSS ALL DATASETS IN CORPUS
# Initialization is AUTOMATED for newly deprecated terms that have 'Replaced By' terms in their ontology files

# Curators should review the monthly 'Curator Report' and add deprecated term replacements to corresponding map if
# 'Replaced By' is not available for a deprecated term.

# If Curators have non-deprecated term changes to apply to all datasets in the corpus where applicable,
# add them here.
ONTOLOGY_TERM_MAPS = {
    "assay": {
    },
    "cell_type": {
    },
    "development_stage": {
    },
    "disease": {
    },
    "organism": {
    },
    "self_reported_ethnicity": {
    },
    "sex": {
    },
    "tissue": {
    },
}

DEPRECATED_FEATURE_IDS = [
]
# fmt: on


def migrate(input_file, output_file, collection_id, dataset_id):
    print(f"Converting {input_file} into {output_file}")

    dataset = ad.read_h5ad(input_file, backed="r")
    if dataset.raw is not None and DEPRECATED_FEATURE_IDS:
        dataset = dataset.to_memory()

    # AUTOMATED, DO NOT CHANGE
    for ontology_name, deprecated_term_map in ONTOLOGY_TERM_MAPS.items():
        utils.replace_ontology_term(dataset.obs, ontology_name, deprecated_term_map)

    # CURATOR-DEFINED, DATASET-SPECIFIC UPDATES
    # Use the template below to define dataset and collection specific ontology changes. Will only apply to dataset
    # if it matches a condition.
    # If no such changes are needed, leave blank
    # Examples:
    # if dataset_id == "<dataset_1_id>":
    #   <no further logic necessary>
    #   utils.replace_ontology_term(df, <ontology_name>, {"term_to_replace": "replacement_term", ...})
    # elif dataset_id == "<dataset_2_id>":
    #   <custom transformation logic beyond scope of util functions>
    # elif collection_id == "<collection_1_id>":
    #   <no further logic necessary>
    #   utils.replace_ontology_term(df, <ontology_name>, {"term_to_replace": "replacement_term", ...})
    # elif collection_id == "<collection_2_id>":
    #   <custom transformation logic beyond scope of replace_ontology_term>
    # ...

    # AUTOMATED, DO NOT CHANGE -- IF GENCODE UPDATED, DEPRECATED FEATURE FILTERING ALGORITHM WILL GO HERE.
    if DEPRECATED_FEATURE_IDS:
        dataset = utils.remove_deprecated_features(dataset, DEPRECATED_FEATURE_IDS)

    dataset.write(output_file, compression="gzip")"""

    with open(os.path.join(mock_target_file), "r") as fp:
        actual_output = fp.read().strip()

    assert actual_output == expected_output[1:]


def test_generate_script__with_automated_replaced_by_map(template, tmpdir):
    ontology_term_map = {
        "assay": {"EFO:0000002": "EFO:0000001"},
        "cell_type": {"CL:0000002": "CL:0000001", "CL:0000004": "CL:0000003"},
        "development_stage": {},
        "disease": {},
        "organism": {},
        "self_reported_ethnicity": {},
        "sex": {},
        "tissue": {},
    }
    gencode_term_map = []
    mock_target_file = tmpdir + "/migrate.py"
    with mock.patch("scripts.migration_assistant.generate_script.target_file", mock_target_file):
        # Execute the function
        generate_script(get_template(), ontology_term_map, gencode_term_map)

        # Verify the output
        expected_output = """
import anndata as ad

from . import utils

# fmt: off
# ONTOLOGY TERMS TO UPDATE ACROSS ALL DATASETS IN CORPUS
# Initialization is AUTOMATED for newly deprecated terms that have 'Replaced By' terms in their ontology files

# Curators should review the monthly 'Curator Report' and add deprecated term replacements to corresponding map if
# 'Replaced By' is not available for a deprecated term.

# If Curators have non-deprecated term changes to apply to all datasets in the corpus where applicable,
# add them here.
ONTOLOGY_TERM_MAPS = {
    "assay": {
        "EFO:0000002": "EFO:0000001", # AUTOMATED
    },
    "cell_type": {
        "CL:0000002": "CL:0000001", # AUTOMATED
        "CL:0000004": "CL:0000003", # AUTOMATED
    },
    "development_stage": {
    },
    "disease": {
    },
    "organism": {
    },
    "self_reported_ethnicity": {
    },
    "sex": {
    },
    "tissue": {
    },
}

DEPRECATED_FEATURE_IDS = [
]
# fmt: on


def migrate(input_file, output_file, collection_id, dataset_id):
    print(f"Converting {input_file} into {output_file}")

    dataset = ad.read_h5ad(input_file, backed="r")
    if dataset.raw is not None and DEPRECATED_FEATURE_IDS:
        dataset = dataset.to_memory()

    # AUTOMATED, DO NOT CHANGE
    for ontology_name, deprecated_term_map in ONTOLOGY_TERM_MAPS.items():
        utils.replace_ontology_term(dataset.obs, ontology_name, deprecated_term_map)

    # CURATOR-DEFINED, DATASET-SPECIFIC UPDATES
    # Use the template below to define dataset and collection specific ontology changes. Will only apply to dataset
    # if it matches a condition.
    # If no such changes are needed, leave blank
    # Examples:
    # if dataset_id == "<dataset_1_id>":
    #   <no further logic necessary>
    #   utils.replace_ontology_term(df, <ontology_name>, {"term_to_replace": "replacement_term", ...})
    # elif dataset_id == "<dataset_2_id>":
    #   <custom transformation logic beyond scope of util functions>
    # elif collection_id == "<collection_1_id>":
    #   <no further logic necessary>
    #   utils.replace_ontology_term(df, <ontology_name>, {"term_to_replace": "replacement_term", ...})
    # elif collection_id == "<collection_2_id>":
    #   <custom transformation logic beyond scope of replace_ontology_term>
    # ...

    # AUTOMATED, DO NOT CHANGE -- IF GENCODE UPDATED, DEPRECATED FEATURE FILTERING ALGORITHM WILL GO HERE.
    if DEPRECATED_FEATURE_IDS:
        dataset = utils.remove_deprecated_features(dataset, DEPRECATED_FEATURE_IDS)

    dataset.write(output_file, compression="gzip")"""

        with open(os.path.join(mock_target_file), "r") as fp:
            actual_output = fp.read().strip()

        assert actual_output == expected_output[1:]


def test_generate_script__with_gencode_changes(template, tmpdir):
    # Prepare the input data
    ontology_term_map = {
        "assay": {},
        "cell_type": {},
        "development_stage": {},
        "disease": {},
        "organism": {},
        "self_reported_ethnicity": {},
        "sex": {},
        "tissue": {},
    }
    mock_target_file = tmpdir + "/migrate.py"
    with mock.patch("scripts.migration_assistant.generate_script.target_file", mock_target_file):
        # Execute the function
        generate_script(
            get_template(),
            ontology_term_map,
            [
                "ENSG00000223972",
                "ENSG00000227232",
            ],
        )

    # Verify the output
    expected_output = """
import anndata as ad

from . import utils

# fmt: off
# ONTOLOGY TERMS TO UPDATE ACROSS ALL DATASETS IN CORPUS
# Initialization is AUTOMATED for newly deprecated terms that have 'Replaced By' terms in their ontology files

# Curators should review the monthly 'Curator Report' and add deprecated term replacements to corresponding map if
# 'Replaced By' is not available for a deprecated term.

# If Curators have non-deprecated term changes to apply to all datasets in the corpus where applicable,
# add them here.
ONTOLOGY_TERM_MAPS = {
    "assay": {
    },
    "cell_type": {
    },
    "development_stage": {
    },
    "disease": {
    },
    "organism": {
    },
    "self_reported_ethnicity": {
    },
    "sex": {
    },
    "tissue": {
    },
}

DEPRECATED_FEATURE_IDS = [
    "ENSG00000223972",
    "ENSG00000227232",
]
# fmt: on


def migrate(input_file, output_file, collection_id, dataset_id):
    print(f"Converting {input_file} into {output_file}")

    dataset = ad.read_h5ad(input_file, backed="r")
    if dataset.raw is not None and DEPRECATED_FEATURE_IDS:
        dataset = dataset.to_memory()

    # AUTOMATED, DO NOT CHANGE
    for ontology_name, deprecated_term_map in ONTOLOGY_TERM_MAPS.items():
        utils.replace_ontology_term(dataset.obs, ontology_name, deprecated_term_map)

    # CURATOR-DEFINED, DATASET-SPECIFIC UPDATES
    # Use the template below to define dataset and collection specific ontology changes. Will only apply to dataset
    # if it matches a condition.
    # If no such changes are needed, leave blank
    # Examples:
    # if dataset_id == "<dataset_1_id>":
    #   <no further logic necessary>
    #   utils.replace_ontology_term(df, <ontology_name>, {"term_to_replace": "replacement_term", ...})
    # elif dataset_id == "<dataset_2_id>":
    #   <custom transformation logic beyond scope of util functions>
    # elif collection_id == "<collection_1_id>":
    #   <no further logic necessary>
    #   utils.replace_ontology_term(df, <ontology_name>, {"term_to_replace": "replacement_term", ...})
    # elif collection_id == "<collection_2_id>":
    #   <custom transformation logic beyond scope of replace_ontology_term>
    # ...

    # AUTOMATED, DO NOT CHANGE -- IF GENCODE UPDATED, DEPRECATED FEATURE FILTERING ALGORITHM WILL GO HERE.
    if DEPRECATED_FEATURE_IDS:
        dataset = utils.remove_deprecated_features(dataset, DEPRECATED_FEATURE_IDS)

    dataset.write(output_file, compression="gzip")"""
    mock.patch("scripts.migration_assistant.generate_script.get_current_version", return_value=expected_output)

    with open(os.path.join(mock_target_file), "r") as fp:
        actual_output = fp.read().strip()

    assert actual_output == expected_output[1:]


@pytest.fixture
def organisms():
    return ["apple", "dog", "mouse"]


def test_get_deprecated_feature_ids(tmp_path, organisms):
    expected_deprecated_feature_ids = []
    for organism in organisms:
        with open(f"{tmp_path}/{organism}_diff.txt", "w") as fp:
            organism_feature_ids = [f"{organism}:{i}" for i in range(4)]
            for feature_id in organism_feature_ids:
                fp.write(feature_id + "\n")
            expected_deprecated_feature_ids.extend(organism_feature_ids)
    with mock.patch("scripts.migration_assistant.generate_script.env.ONTOLOGY_DIR", tmp_path):
        actual_deprecated_features = get_deprecated_feature_ids()
    expected_deprecated_feature_ids.sort()
    actual_deprecated_features.sort()
    assert expected_deprecated_feature_ids == actual_deprecated_features


def test_get_deprecated_feature_ids__no_files(tmp_path):
    with mock.patch("scripts.migration_assistant.generate_script.env.ONTOLOGY_DIR", tmp_path):
        actual_deprecated_features = get_deprecated_feature_ids()
    assert actual_deprecated_features == []


def test_get_deprecated_feature_ids__empty_feature_files(tmp_path, organisms):
    for organism in organisms:
        with open(f"{tmp_path}/{organism}_diff.txt", "w") as fp:
            fp.write("")
    with mock.patch("scripts.migration_assistant.generate_script.env.ONTOLOGY_DIR", tmp_path):
        actual_deprecated_features = get_deprecated_feature_ids()
    assert actual_deprecated_features == []
