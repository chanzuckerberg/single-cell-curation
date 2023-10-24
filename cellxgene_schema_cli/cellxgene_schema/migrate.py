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

    if collection_id == "ced320a1-29f3-47c1-a735-513c7084d508":
        utils.replace_ontology_term(dataset.obs, "self_reported_ethnicity",
            {"HANCESTRO:0487": "HANCESTRO:0598", "HANCESTRO:0496": "HANCESTRO:0597"})

    if dataset.uns["title"] == "adipose from neck - all nuclei":
        utils.map_ontology_term(dataset.obs, "self_reported_ethnicity", "donor_id",
            {"CREEK003": "HANCESTRO:0005,HANCESTRO:0013"})

    if collection_id == "1ca90a2d-2943-483d-b678-b809bf464c30":
        utils.map_ontology_term(dataset.obs, "self_reported_ethnicity", "donor_id",
            {"H20.33.018": "HANCESTRO:0013,HANCESTRO:0014",
            "H20.33.034": "HANCESTRO:0005",
            "H21.33.037": "HANCESTRO:0005"})

    if collection_id == "6f6d381a-7701-4781-935c-db10d30de293":
        new_eth = "HANCESTRO:0005,HANCESTRO:0008"
        study = "homosapiens_None_2023_None_sikkemalisa_"
        donor = "_d10_1101_2022_03_10_4837472020-3173-NC004"
        if dataset_id == "066943a2-fdac-4b29-b348-40cede398e4e":
            utils.map_ontology_term(dataset.obs, "self_reported_ethnicity", "donor_id",
                {study + "001" + donor: new_eth})
        elif dataset_id == "9f222629-9e39-47d0-b83f-e08d610c7479":
            utils.map_ontology_term(dataset.obs, "self_reported_ethnicity", "donor_id",
                {study + "002" + donor: new_eth})

    if collection_id == "7d7cabfd-1d1f-40af-96b7-26a0825a306d":
        donors_to_update = ["Rep_C_1012","Rep_C_1017","Rep_C_1033","Rep_C_1036","Rep_C_1037",
            "Rep_C_1039","Rep_C_1050","Rep_C_1053","Rep_C_1055","Rep_C_1059","Rep_C_1060",
            "Rep_C_1064","Rep_C_1076","Rep_C_1078","Rep_C_1094","Rep_C_1095","Rep_C_1107",
            "Rep_C_1143","Rep_C_1151","Rep_C_1154","Rep_C_1161"]
        update_map = {d:"HANCESTRO:0014" for d in donors_to_update}
        utils.map_ontology_term(dataset.obs, "self_reported_ethnicity", "donor_id", update_map)

    if collection_id == "4d74781b-8186-4c9a-b659-ff4dc4601d91":
        if dataset_id == "b07fb54c-d7ad-4995-8bb0-8f3d8611cabe":
            utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0000234": "CL:0000113"})

    if "multiethnic" in dataset.obs["self_reported_ethnicity"].unique():
        utils.replace_ontology_term(dataset.obs, "self_reported_ethnicity", {"multiethnic": "unknown"})

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
