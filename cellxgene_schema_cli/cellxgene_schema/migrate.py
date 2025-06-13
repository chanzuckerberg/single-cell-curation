import dask

from . import utils

# fmt: off
# ONTOLOGY TERMS TO UPDATE ACROSS ALL DATASETS IN CORPUS
# Initialization is AUTOMATED for newly deprecated terms that have 'Replaced By' terms in their ontology files

# Curators should review the monthly 'Curator Report' and add deprecated term replacements to corresponding map if
# 'Replaced By' is not available for a deprecated term.

# If Curators have non-deprecated term changes to apply to all datasets in the corpus where applicable,
# add them here.
ONTOLOGY_TERM_OBS_MAPS = {
    "assay": {
    },
    "cell_type": {
    },
    "development_stage": {
    },
    "disease": {
    },
    "self_reported_ethnicity": {
    },
    "sex": {
    },
    "tissue": {
    },
}

ONTOLOGY_TERM_UNS_MAPS = {
    "organism": {
    },
}

DEPRECATED_FEATURE_IDS = [
]

# Dictionary for CURATOR-DEFINED remapping of deprecated feature IDs, if any, to new feature IDs.
GENCODE_MAPPER = {}
# fmt: on


def migrate(input_file, output_file, collection_id, dataset_id):
    print(f"Converting {input_file} into {output_file}")

    dataset = utils.read_h5ad(input_file)

    # Migrate organism_ontology_term_id for 6.0.0
    dataset = utils.move_column_from_obs_to_uns(adata=dataset, column_name="organism_ontology_term_id")

    # Migrate delimter from `,` to ` || ` for self_reported_ethnicity_ontology_term_id
    utils.replace_delimiter(
        dataframe=dataset.obs,
        old_delimiter=",",
        new_delimiter=" || ",
        column_name="self_reported_ethnicity_ontology_term_id",
    )

    # AUTOMATED, DO NOT CHANGE
    for ontology_name, deprecated_term_map in ONTOLOGY_TERM_OBS_MAPS.items():
        utils.replace_ontology_term(dataset.obs, ontology_name, deprecated_term_map)
    for ontology_name, deprecated_term_map in ONTOLOGY_TERM_UNS_MAPS.items():
        dataset = utils.replace_ontology_term_uns(dataset, ontology_name, deprecated_term_map)

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

    if (
        dataset.uns["title"]
        == "Time-resolved single-cell analysis of Brca1 associated mammary tumourigenesis reveals aberrant differentiation of luminal progenitors"
    ):
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0000451": "CL:4047054"})

    if collection_id == "126afc71-47fb-4e9d-8aaf-9d9f61e0ac77":
        utils.map_ontology_term(
            dataset.obs,
            "disease",
            "disease_group",
            {
                "Healthy": "PATO:0000461",
                "HIV": "MONDO:0005109",
                "VL_HIV": "MONDO:0005109 || MONDO:0005445",
                "aL_HIV": "MONDO:0005109 || MONDO:0011989",
                "cVL_HIV": "MONDO:0005109 || MONDO:0005445",
                "pVL_HIV": "MONDO:0005109 || MONDO:0005445",
            },
        )

    if collection_id == "eeba7193-4d32-46bd-a21b-089936d60601":
        utils.map_ontology_term(
            dataset.obs,
            "disease",
            "sample_id",
            {
                "mraf1": "MONDO:0004981 || MONDO:1030008",
                "mraf2": "MONDO:0004981 || MONDO:1030008",
                "mraf3": "MONDO:0004981 || MONDO:1030008",
                "mraf4": "MONDO:0004981 || MONDO:1030008",
                "mraf7": "MONDO:0004981 || MONDO:1030008",
                "mraf13": "MONDO:0004981 || MONDO:1030008",
                "mraf14": "MONDO:0004981 || MONDO:1030008",
            },
        )

    # errors found in first dev migration run
    # 'organism' already in uns
    if dataset_id == "de2c780c-1747-40bd-9ccf-9588ec186cee":
        del dataset.uns["organism"]

    # deprecated CL term found in 3 datasets in one collection
    if dataset_id in [
        "773b9b2e-70c8-40be-8cbb-e7b5abab360d",
        "e5233a94-9e43-418c-8209-6f1400c31530",
        "7bb64315-9e5a-41b9-9235-59acf9642a3e",
    ]:
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0000651": "CL:0002181"})

    # 2 other datasets with deprecated CL term
    if dataset_id in ["e40c6272-af77-4a10-9385-62a398884f27", "d6dfdef1-406d-4efb-808c-3c5eddbfe0cb"]:
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0000212": "CL:0000677"})

    # feature_is_filtered set incorrectly for dataset with no raw.X matrix, should be all False
    if dataset_id == "ee195b7d-184d-4dfa-9b1c-51a7e601ac11":
        dataset.var["feature_is_filtered"] = False

    # 11 datasets with "organism_ontology_term_id_colors" in uns, all can be removed
    # seems easy enough to quickly check for key, can also list out the 11 dataset_ids
    if "organism_ontology_term_id_colors" in dataset.uns:
        del dataset.uns["organism_ontology_term_id_colors"]

    if GENCODE_MAPPER:
        dataset = utils.remap_deprecated_features(adata=dataset, remapped_features=GENCODE_MAPPER)

    # AUTOMATED, DO NOT CHANGE -- IF GENCODE UPDATED, DEPRECATED FEATURE FILTERING ALGORITHM WILL GO HERE.
    if DEPRECATED_FEATURE_IDS:
        dataset = utils.remove_deprecated_features(adata=dataset, deprecated=DEPRECATED_FEATURE_IDS)

    with dask.config.set(scheduler="single-threaded"):
        dataset.write_h5ad(output_file, compression="gzip")
