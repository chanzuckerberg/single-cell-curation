import json
import os

import dask

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
        "EFO:0010961" : "EFO:0022857", # AUTOMATED VISIUM TERM UPDATE
    },
    "cell_type": { 
        "CL:0000402": "CL:0000099", # AUTOMATED
        "CL:0010003": "CL:0000322", # AUTOMATED
        "CL:0000555": "CL:4023161", # AUTOMATED
    },
    "development_stage": DEV_STAGE_AUTO_MIGRATE_MAP,
    "disease": {
    },
    "self_reported_ethnicity": {
    },
    "sex": {
    },
    "tissue": {
        "CL:0010003": "CL:0000322", # AUTOMATED
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

with open(os.path.join(BASE_DIR,'migrate_files/non_csr_list.txt'), 'r') as f:
    text = f.read()
    non_csr_list = text.split("\n")
    f.close()

with open(os.path.join(BASE_DIR, 'migrate_files/private_non_csr.json'), 'r') as f:
   private_non_csr_list = [tuple(x) for x in json.load(f)]
   f.close()

# Dictionary for CURATOR-DEFINED remapping of deprecated feature IDs, if any, to new feature IDs.
GENCODE_MAPPER = {}


def migrate(input_file, output_file, collection_id, dataset_id):
    print(f"Converting {input_file} into {output_file}")

    dataset = utils.read_h5ad(input_file)

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

    if GENCODE_MAPPER:
        dataset = utils.remap_deprecated_features(adata=dataset, remapped_features=GENCODE_MAPPER)

    # AUTOMATED, DO NOT CHANGE -- IF GENCODE UPDATED, DEPRECATED FEATURE FILTERING ALGORITHM WILL GO HERE.
    if DEPRECATED_FEATURE_IDS:
        dataset = utils.remove_deprecated_features(adata=dataset, deprecated=DEPRECATED_FEATURE_IDS)

    # Schema 5.2 Migration
    # https://github.com/chanzuckerberg/single-cell-curation/issues/958
    '''

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

    '''

    # Schema 5.3 Migration
    # https://github.com/chanzuckerberg/single-cell-curation/issues/1023

    # Checking pre-identified datasets with non-csr matrices for sparsity
    # If sparsity is high enough, converting to sparse_csr
    if dataset_id in non_csr_list:
        dataset = utils.check_non_csr_matrices(dataset)

    # Same check, but for private datasets
    if dataset.uns['title'] in list(zip(*private_non_csr_list))[0]:
        dataset_donors = dataset.obs['donor_id'].unique()
        for d in dataset_donors:
            if (dataset.uns["title"],d) in private_non_csr_list:
                dataset = utils.check_non_csr_matrices(dataset)

    # Public dataset changes
    if collection_id == "33f48a52-31d8-4cc8-bd00-1e89c659a87f":
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0000065": "CL:4052001"})
        utils.replace_ontology_term(dataset.obs, "tissue", {"CL:0000065":"CL:4052001"})

    if collection_id == "ef9f0fae-c12c-4f1a-bb67-63346c2ae0f1":
        utils.map_ontology_term(dataset.obs,
         "cell_type",
         "author_cell_type", 
         {
            "Type IIX": "CL:4052026",
            "NMJ":"CL:4052025",
            "MTJ":"CL:4052024"
         })

    if collection_id == "964581c3-b5ac-494d-9bd4-0dcb6fb058da":
        utils.map_ontology_term(dataset.obs,
         "cell_type",
         "author_cell_type", 
         {
            "Intercalated": "CL:4052048",
            "Striated": "CL:4052049"
         })

    if collection_id == "350fe77e-cc64-482e-a42a-18bec9e4a629":
        utils.map_ontology_term(dataset.obs, "cell_type", "cluster.names", {"Early Theca": "CL:4052010"})

    if collection_id == "1e313b15-4aca-4e5a-94a3-1093c4c42abd":
        utils.replace_ontology_term(dataset.obs, "assay", {"EFO:0008953":"EFO:0022845"})

    if collection_id == "4c6eaf5c-6d57-4c76-b1e9-60df8c655f1e":
        utils.map_ontology_term(dataset.obs,
         "development_stage",
         "donor_id", 
         {
            "MMD_22_20137": "HsapDv:0000274",
            "A23_5566": "HsapDv:0000274",
            "19-D011": "HsapDv:0000274",
            "D018-13": "HsapDv:0000274",
            "19-D006": "HsapDv:0000274"
         })

    if collection_id == "8f126edf-5405-4731-8374-b5ce11f53e82":
        utils.map_ontology_term(dataset.obs,
         "development_stage",
         "donor_id", 
         {
            "S00006": "HsapDv:0000274",
            "S00063": "HsapDv:0000274"
         })

    if collection_id == '606dc65c-ebd8-4f93-ac86-d0789c04863f':
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0002503":"CL:4052030"})

    # Private dataset changes
    # Changes for datasets in private collection "HTAN CHOP - A longitudinal single-cell atlas of pediatric high-grade glioma."

    if dataset.uns["title"] == "Longitudinal single-nucleus RNA-seq atlas of pediatric high-grade glioma":
        utils.replace_ontology_term(dataset.obs, "disease", {"MONDO:0100342":"MONDO:1010030"})

    if dataset.uns["title"] == "Longitudinal single-nucleus RNA-seq atlas of neoplastic cells in pediatric high-grade glioma":
        utils.replace_ontology_term(dataset.obs, "disease", {"MONDO:0100342":"MONDO:1010030"})

    if dataset.uns["title"] == "Longitudinal single-nucleus RNA-seq atlas of myeloid cells in pediatric high-grade glioma":
        utils.replace_ontology_term(dataset.obs, "disease", {"MONDO:0100342":"MONDO:1010030"})

    # Changes for datasets in private collection "Human Breast Cancer Single Cell Atlas"

    if dataset.uns["title"] == "Global Atlas" and 'BC17086-12_Tumor' in dataset.obs['donor_id'].unique():
        utils.map_ontology_term(dataset.obs, "development_stage", "donor_id", {"TBB330": "HsapDv:0000274"})

    if dataset.uns["title"] == "Immune Compartment" and 'BC17086-12_Tumor' in dataset.obs['donor_id'].unique():
        utils.map_ontology_term(dataset.obs, "development_stage", "donor_id", {"TBB330": "HsapDv:0000274"})

    if dataset.uns["title"] == "Epithelial Compartment" and 'BC17086-12_Tumor' in dataset.obs['donor_id'].unique():
        utils.map_ontology_term(dataset.obs, "development_stage", "donor_id", {"TBB330": "HsapDv:0000274"})

    if dataset.uns["title"] == "Stromal Compartment" and 'wu_natgen_CID44041' in dataset.obs['donor_id'].unique():
        utils.map_ontology_term(dataset.obs, "development_stage", "donor_id", {"TBB330": "HsapDv:0000274"})

    #  Changes for dataset in private collection "Single cell RNA sequencing of the human embryonic meninges at 5-13 weeks post conception"

    if dataset.uns["title"] == "All_cells" and 'XHU364' in dataset.obs['donor_id'].unique():
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0011026":"CL:4042021"})

    if dataset.uns["title"] == "Fibroblasts" and 'XHU364' in dataset.obs['donor_id'].unique():
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0011026":"CL:4042021"})

    if dataset.uns["title"] == "Erythropoietic" and 'XHU364' in dataset.obs['donor_id'].unique():
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0011026":"CL:4042021"})

    if dataset.uns["title"] == "Epithelial" and 'XDD400' in dataset.obs['donor_id'].unique():
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0011026":"CL:4042021"})

    if dataset.uns["title"] == "Neural" and 'XHU364' in dataset.obs['donor_id'].unique():
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0011026":"CL:4042021"})

    if dataset.uns["title"] == "Immune" and 'XHU364' in dataset.obs['donor_id'].unique():
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0011026":"CL:4042021"})

    if dataset.uns["title"] == "Vascular_endothelial" and 'XHU364' in dataset.obs['donor_id'].unique():
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0011026":"CL:4042021"})

    if dataset.uns["title"] == "Vascular_perivascular" and 'XHU364' in dataset.obs['donor_id'].unique():
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0011026":"CL:4042021"})

    if dataset.uns["title"] == "Neural_crest" and 'XHU364' in dataset.obs['donor_id'].unique():
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0011026":"CL:4042021"})

    with dask.config.set(scheduler="single-threaded"):
        dataset.write_h5ad(output_file, compression="gzip")
