import anndata as ad
import pandas as pd
from . import utils
import dask

# fmt: off
# ONTOLOGY TERMS TO UPDATE ACROSS ALL DATASETS IN CORPUS
ONTOLOGY_TERM_MAPS = {
    "assay": {
        "EFO:0010961" : "EFO:0022857", # AUTOMATED VISIUM TERM UPDATE
    },
    "cell_type": { 
        "CL:0000402": "CL:0000099", # AUTOMATED
        "CL:0010003": "CL:0000322", # AUTOMATED
        "CL:0000555": "CL:4023161", # AUTOMATED
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
        "CL:0010003": "CL:0000322", # AUTOMATED
    },
}

df = pd.read_csv('migrate_files/non_csr_list.csv')
NON_CSR_DATASETS = df['dataset_id'].to_list()

# fmt: on

# Dictionary for CURATOR-DEFINED remapping of deprecated feature IDs, if any, to new feature IDs.
GENCODE_MAPPER = {}


def migrate(input_file, output_file, collection_id, dataset_id):
    print(f"Converting {input_file} into {output_file}")

    dataset = utils.read_h5ad(input_file)

    # AUTOMATED, DO NOT CHANGE
    for ontology_name, deprecated_term_map in ONTOLOGY_TERM_MAPS.items():
        utils.replace_ontology_term(dataset.obs, ontology_name, deprecated_term_map)


    # https://github.com/chanzuckerberg/single-cell-curation/issues/1023

    # Checking pre-identified datasets with non-csr matrices for sparsity
    # If sparsity is high enough, converting to sparse_csr
    if dataset_id in NON_CSR_DATASETS:
        dataset = dataset.to_memory()
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
            "MRJ":"CL:4052024"
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

    # Private dataset changes

    # Changes for dataset in private collection "Mouse Post-Flu Time Series"
    if dataset.uns["title"] == "Mouse Post-Flu Time Series":
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0002503":"CL:4052030"})

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
