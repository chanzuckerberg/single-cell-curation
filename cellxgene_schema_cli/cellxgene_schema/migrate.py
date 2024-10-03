import json
import os

import anndata as ad

from . import utils

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# fmt: off
# ONTOLOGY TERMS TO UPDATE ACROSS ALL DATASETS IN CORPUS
# Initialization is AUTOMATED for newly deprecated terms that have 'Replaced By' terms in their ontology files

# Curators should review the monthly 'Curator Report' and add deprecated term replacements to corresponding map if
# 'Replaced By' is not available for a deprecated term.

# If Curators have non-deprecated term changes to apply to all datasets in the corpus where applicable,
# add them here.

# We have used this mapping for automigrating terms, not sure if this is accounted for in other migration logic
with open(os.path.join(BASE_DIR, "migrate_files/automigrate_terms.json"), "r") as file:
    DEV_STAGE_AUTO_MIGRATE_MAP = json.load(file)

ONTOLOGY_TERM_MAPS = {
    "assay": {
    },
    "cell_type": {
        "CL:4023070": "CL:4023064", # AUTOMATED
    },
    "development_stage": DEV_STAGE_AUTO_MIGRATE_MAP,
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

# Dictionary for CURATOR-DEFINED remapping of deprecated feature IDs, if any, to new feature IDs.
GENCODE_MAPPER = {}

with open(os.path.join(BASE_DIR, "migrate_files/donor_updates.json"), "r") as file:
    DONOR_DEV_STAGE_MAP = json.load(file)

with open(os.path.join(BASE_DIR, "migrate_files/title_donor_updates.json"), "r") as file:
    TITLE_DONOR_DEV_STAGE_MAP = json.load(file)
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

    # logic for mapping donor updates
    if collection_id in DONOR_DEV_STAGE_MAP:
        collection_donor_dev_map = DONOR_DEV_STAGE_MAP[collection_id]
        utils.map_ontology_term(dataset.obs, "development_stage", "donor_id", collection_donor_dev_map)

    # private dataset titles in title_donor_updates.json should be unique in corpus
    if dataset.uns["title"] in TITLE_DONOR_DEV_STAGE_MAP:
        dataset_donor_dev_map = TITLE_DONOR_DEV_STAGE_MAP[dataset.uns["title"]]
        utils.map_ontology_term(dataset.obs, "development_stage", "donor_id", dataset_donor_dev_map)

    # https://github.com/chanzuckerberg/single-cell-curation/issues/958

    if collection_id == "c114c20f-1ef4-49a5-9c2e-d965787fb90c":
        utils.replace_ontology_term(dataset.obs, "assay", {"EFO:0010550": "EFO:0030028"})

    if collection_id == "45d5d2c3-bc28-4814-aed6-0bb6f0e11c82":
        utils.replace_ontology_term(dataset.obs, "assay", {"EFO:0010550": "EFO:0030028"})

    if collection_id == "962df42d-9675-4d05-bc75-597ec7bf4afb":
        utils.replace_ontology_term(dataset.obs, "development_stage", {"unknown": "HsapDv:0000264"})

    if collection_id == "48d354f5-a5ca-4f35-a3bb-fa3687502252":
        utils.map_ontology_term(dataset.obs, "development_stage", "var_time", {"P7": "MmusDv:0000117"})

    if collection_id == "613f5480-4957-4f80-b804-0e2b85ac454c":
        utils.map_ontology_term(
            dataset.obs,
            "development_stage",
            "age",
            {
                "P4": "MmusDv:0000114",
                "P7": "MmusDv:0000117",
            },
        )

    if collection_id == "d86517f0-fa7e-4266-b82e-a521350d6d36":
        utils.map_ontology_term(
            dataset.obs,
            "development_stage",
            "Dataset",
            {
                "RomanovDev10x": "MmusDv:0000144",
                "Mickelsen10x": "MmusDv:0000149",
                "Flynn10x": "MmusDv:0000149",
            },
        )
        utils.map_ontology_term(
            dataset.obs,
            "development_stage",
            "SRA_ID",
            {
                "SRR5164436": "MmusDv:0000178",
                "SRR5164437": "MmusDv:0000178",
                "SRR5164438": "MmusDv:0000149",
            },
        )

    # dev migration revealed column 'feature_type' in var and raw.var, titles should be unique in corpus
    titles_with_feature_type = [
        "Healthy reference",
        "Extended - Neural cells",
        "Extended - Endothelial cells",
        "Extended+ - 18485 genes",
        "Healthy - Small intestine (adult/pediatric)",
        "Healthy - Mesenchymal (adult/pediatric)",
        "Healthy - Stomach",
        "Healthy - Large Intestine (adult/pediatric)",
        "Healthy - Mesenchymal (first trimester)",
        "Extended - All genes",
        "Healthy - Small Intestine (first trimester)",
        "Healthy - Neural cells",
        "Healthy - T and NK cells",
        "Healthy - Myeloid cells",
        "Extended - Small intestine (adult/pediatric)",
        "Extended - T and NK cells",
        "Healthy - Salivary gland",
        "Extended - Large Intestine (adult/pediatric)",
        "Healthy - Large Intestine (first trimester)",
        "Healthy - Oesophagus",
        "Extended - Stomach",
        "Healthy - Mesenchymal (second trimester)",
        "Extended - Small intestine (adult/pediatric)- 18485genes",
        "Healthy - Oral mucosa",
        "Extended - B and B plasma cells",
        "Extended - Myeloid_with_neutrophils",
        "Healthy - B and B plasma cells",
        "Healthy - Small Intestine (second trimester)",
        "Healthy - Endothelial cells",
        "Extended - Myeloid cells",
        "Extended - Mesenchymal (adult/pediatric)",
        "Healthy - Large Intestine (second trimester)",
    ]

    if dataset.uns["title"] in titles_with_feature_type or collection_id == "e5f58829-1a66-40b5-a624-9046778e74f5":
        if dataset.raw:
            dataset.raw.var.drop(columns="feature_type", inplace=True)

        dataset.var.drop(columns="feature_type", inplace=True)

    if GENCODE_MAPPER:
        dataset = utils.remap_deprecated_features(adata=dataset, remapped_features=GENCODE_MAPPER)

    # AUTOMATED, DO NOT CHANGE -- IF GENCODE UPDATED, DEPRECATED FEATURE FILTERING ALGORITHM WILL GO HERE.
    if DEPRECATED_FEATURE_IDS:
        dataset = utils.remove_deprecated_features(adata=dataset, deprecated=DEPRECATED_FEATURE_IDS)

    dataset.write(output_file, compression="gzip")
