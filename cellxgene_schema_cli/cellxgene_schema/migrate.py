import re

import anndata as ad
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

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

    if "schema_version" in dataset.uns:
        del dataset.uns["schema_version"]

    dataset.obs["tissue_type"] = "tissue"
    dataset.obs["tissue_type"] = np.where(
        dataset.obs.tissue_ontology_term_id.str.contains("(cell culture)"), "cell culture", dataset.obs["tissue_type"]
    )
    dataset.obs["tissue_type"] = np.where(
        dataset.obs.tissue_ontology_term_id.str.contains("(organoid)"), "organoid", dataset.obs["tissue_type"]
    )
    dataset.obs["tissue_ontology_term_id"] = np.where(
        dataset.obs.tissue_ontology_term_id.str.contains("(organoid)"),
        dataset.obs["tissue_ontology_term_id"].str.replace(r" \(organoid\)", ""),
        dataset.obs["tissue_ontology_term_id"],
    )
    dataset.obs["tissue_ontology_term_id"] = np.where(
        dataset.obs.tissue_ontology_term_id.str.contains("(cell culture)"),
        dataset.obs["tissue_ontology_term_id"].str.replace(r" \(cell culture\)", ""),
        dataset.obs["tissue_ontology_term_id"],
    )

    dataset.obs["tissue_ontology_term_id"] = pd.Categorical(dataset.obs.tissue_ontology_term_id)
    dataset.obs["tissue_type"] = pd.Categorical(dataset.obs.tissue_type)

    # Colors Validation, remove _colors key if invalid
    # Mapping from obs column name to number of unique categorical values
    category_mapping = {}

    # Check for categorical dtypes in the dataframe directly
    for column_name in dataset.obs.columns:
        column = dataset.obs[column_name]
        if column.dtype.name == "category":
            category_mapping[column_name] = column.nunique()

    def color_key_is_valid(key: str, value: str) -> bool:
        column_name = key.replace("_colors", "")
        obs_unique_values = category_mapping.get(column_name)
        if not obs_unique_values:
            return False
        if value is None or not isinstance(value, np.ndarray):
            return False
        if not all(isinstance(color, str) for color in value):
            return False
        if len(value) < obs_unique_values:
            return False
        all_hex_colors = all(re.match(r"^#([0-9a-fA-F]{6})$", color) for color in value)
        all_css4_colors = all(color in mcolors.CSS4_COLORS for color in value)
        if not (all_hex_colors or all_css4_colors):
            return False
        return True

    for key, value in dataset.uns.items():
        if key.endswith("_colors") and not color_key_is_valid(key, value):
            del dataset.uns[key]

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

    dataset.write(output_file, compression="gzip")
