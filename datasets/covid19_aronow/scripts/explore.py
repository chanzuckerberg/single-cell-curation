import scanpy as sc

dataset = sc.read('./data/original/PBMC_merged_normalized_addRaw_1222.h5ad', 'r')

# Tissue
dataset.obs[['tissue_ontology_term_id', 'tissue']].drop_duplicates()

# assay
dataset.obs[['disease_ontology_term_id', 'disease']].drop_duplicates()

# cell type
dataset.obs[['cell_type_ontology_term_id', 'Cell.class']].drop_duplicates()

# ethnicity
dataset.obs[['ethnicity_ontology_term_id']].drop_duplicates()

# ethnicity
dataset.obs[['development_stage_ontology_term_id', 'age(y)']].drop_duplicates()
