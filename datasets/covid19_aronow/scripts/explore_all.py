import scanpy as sc
import os

working_dir = './data/original/'

dataset_files = [os.path.join(working_dir, i) for i in os.listdir(working_dir) if 'h5ad' in i]

for i in dataset_files:
    
    print(i)
    
    dataset = sc.read(i, 'r')
    
    print(dataset.obs[['sex']].drop_duplicates())
    
    print()
    print()
    print()
    print()
    
    # Raw layer
    print(dataset.raw)
    print('raw' in dataset.layers)
    
    # cell type
    print(dataset.obs[['cell_type_ontology_term_id']].drop_duplicates())

    # Tissue
    print(dataset.obs[['tissue_ontology_term_id']].drop_duplicates())

    # disease
    print(dataset.obs[['disease_ontology_term_id']].drop_duplicates())
    
    # asaay
    print(dataset.obs[['assay_ontology_term_id']].drop_duplicates())

    # ethnicity
    print(dataset.obs[['ethnicity_term_ontology_id']].drop_duplicates())

    # dev
    print(dataset.obs[['development_stage_ontology_term_id']].drop_duplicates())
    
    # sex
    print(dataset.obs[['sex']].drop_duplicates())
    
    # embeddings
    print(dataset.obsm)
    
    print()
    print()
