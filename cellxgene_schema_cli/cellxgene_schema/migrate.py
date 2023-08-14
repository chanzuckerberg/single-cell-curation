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
        "CL:0000663": "CL:1000147", # AUTOMATED
        "CL:0000391": "CL:0000394", # AUTOMATED
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
    if collection_id == "a48f5033-3438-4550-8574-cdff3263fdfd":
        utils.replace_ontology_term(dataset.obs, "assay", {"EFO:0008913": "EFO:0700010"})
    elif collection_id == "64b24fda-6591-4ce1-89e7-33eb6c43ad7b": #REMOVE RACHADELE'S COMMENT - DUPLICATED
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0000677": "CL:4030026"})
    elif collection_id == "edb893ee-4066-4128-9aec-5eb2b03f8287":
        utils.replace_ontology_term(dataset.obs, "assay", {"EFO:0010183": "EFO:0700011"})
    elif collection_id == "48259aa8-f168-4bf5-b797-af8e88da6637" and dataset_id in [\
            "0ba636a1-4754-4786-a8be-7ab3cf760fd6", "975e13b6-bec1-4eed-b46a-9be1f1357373"]:
            utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0011026": "CL:0009116"})
    elif dataset.uns["title"] == "Single nucleus transcriptomic profiling of human healthy hamstring tendon":
        utils.replace_ontology_term(dataset.obs, "tissue", {"UBERON:0000043": "UBERON:8480009"}) #remove coll link in git
    elif collection_id == "1d1c7275-476a-49e2-9022-ad1b1c793594" and dataset_id == \
            "5cdbb2ea-c622-466d-9ead-7884ad8cb99f":
            utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0000561": "CL:4030027", "CL:1001509": "CL:4030028"})
    elif collection_id == "939769a8-d8d2-4d01-abfc-55699893fd49":
        if dataset_id == "f8c77961-67a7-4161-b8c2-61c3f917b54f":
            update_map = {"1.0": "CL:0004219", "5.0": "CL:4030028", "11.0": "CL:0004232"}
            utils.map_ontology_term(dataset.obs, "cell_type", "author_cell_type", update_map)
        elif dataset_id == "cec9f9a5-8832-437d-99af-fb8237cde54b":
            update_map = {"1.0": "CL:4023033", "4.0": "CL:4023033", "2.0": "CL:4023032", "5.0": "CL:4023032"}
            utils.map_ontology_term(dataset.obs, "cell_type", "author_cell_type", update_map)
        elif dataset_id == "d95ab381-2b7c-4885-b168-0097ed4e397f":
            update_map = {"2.0": "CL:0003050"}
            utils.map_ontology_term(dataset.obs, "cell_type", "author_cell_type", update_map)
    elif collection_id == "a18474f4-ff1e-4864-af69-270b956cee5b":
        utils.replace_ontology_term(dataset.obs, "disease", {"MONDO:0004298": "MONDO:0100190"})
    elif dataset.uns["title"] in [
        "All annotated cells from infection 1",
        "All annotated cells from infection 2",
        "All annotated cells from infection 3",
        "All annotated cells from infection 4"
        ]:
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0002144": "CL:4028001", "CL:0000669": "CL:0009089"})
    elif dataset.uns["title"] == "Transcriptional responses of the human dorsal striatum in opioid use disorder implicates cell type-specifc programs":
            update_map = {"D2-Matrix": "CL:4023029", "D2-Striosome": "CL:4023029", "D1-Matrix": "CL:4023026", "D1-Striosome": "CL:4023026"}
            utils.map_ontology_term(dataset.obs, "cell_type", "author_cell_type", update_map)
    elif dataset.uns["title"] in [
        'A single-cell multi-omic atlas spanning the adult rhesus macaque brain',
        'A single-cell multi-omic atlas spanning the adult rhesus macaque brain (1.5 million cell subset)',
        'A single-cell multi-omic atlas spanning the adult rhesus macaque brain (cell class "vascular cells" subset)',
        ]:
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0002319": "CL:4023072"})
    elif collection_id == "e3aa612b-0d7d-4d3f-bbea-b8972a74dd4b" and dataset_id == \
            "812fa7bd-db15-4357-b2c9-efc8e1eb0450":
            utils.replace_ontology_term(dataset.obs, "tissue", {"UBERON:0000956": "UBERON:8440075"})
            utils.replace_ontology_term(dataset.obs, "assay", {"EFO:0010184": "EFO:0700016"})
    elif collection_id == "9d63fcf1-5ca0-4006-8d8f-872f3327dbe9" and dataset_id in [\
            "c05e6940-729c-47bd-a2a6-6ce3730c4919","12194ced-8086-458e-84a8-e2ab935d8db1"]:
            update_map = {"COP_A": "CL:4023059", "COP_B": "CL:4023059", "COP_C": "CL:4023059"}
            utils.map_ontology_term(dataset.obs, "cell_type", "author_cell_type", update_map)

    # AUTOMATED, DO NOT CHANGE -- IF GENCODE UPDATED, DEPRECATED FEATURE FILTERING ALGORITHM WILL GO HERE.
    if DEPRECATED_FEATURE_IDS:
        dataset = utils.remove_deprecated_features(dataset, DEPRECATED_FEATURE_IDS)

    dataset.write(output_file, compression="gzip")
# fmt: on
