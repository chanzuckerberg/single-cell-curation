import scanpy as sc
import pandas as pd
import ontology, os, sys


in_dir = sys.argv[1]
out_dir = sys.argv[2]
files = ['COVID-19_PBMC_Arunachalam_et_al.h5ad', 'COVID-19_PBMC_Guo_et_al.h5ad', 'COVID-19_PBMC_Lee_et_al.h5ad', 'COVID-19_PBMC_Schulte-Schrepping_et_al.h5ad', 'COVID-19_PBMC_Wilk_et_al.h5ad', 'B_cell_subcluster.h5ad',  'T_NK_Subcluster_0121.h5ad', 'monocyte_subcluster_0126.h5ad', 'platelets_subcluster.h5ad', 'PBMC_merged_normalized_addRaw_coordinatesFixed_0126.h5ad']
files = ['T_NK_Subcluster_0121.h5ad']

for c_file in files:
    
    print("Working with:", c_file)
    
    dataset = sc.read(os.path.join(in_dir, c_file))

    columns = ['tissue', 'disease', 'cell_type', 'self_reported_ethnicity', 'development_stage']
    suffix = '_ontology_term_id'

    dataset.uns['version'] = {'corpora_schema_version': '1.1.0', 'corpora_encoding_version': '0.1.0'}
    dataset.uns['organism_ontology_term_id'] = "NCBITaxon:9606"
    dataset.uns['organism'] = ontology.get_ontology_label(dataset.uns['organism_ontology_term_id'])
    
    # correct typo in key
    if 'layer_description' in dataset.uns:
        dataset.uns['layer_descriptions'] = dataset.uns['layer_description']
        del dataset.uns['layer_description']
        
    # Eliminate empty clusters in T_NK_Subcluster_0121.h5ad
    if c_file == 'T_NK_Subcluster_0121.h5ad':
        #to_keep = dataset.obs.loc[dataset.obs['sub_cluster'] != '', ].index
        #dataset = dataset[to_keep,:]
        dataset.obs.drop(columns='sub_cluster', inplace=True)

    # Add labels for ontologies
    for column in columns:
        dataset.obs[column + suffix].cat.add_categories("", inplace=True)
        uniq = dataset.obs[[column + suffix]].drop_duplicates()
        uniq[column] = ''
        for i in range(len(uniq.index)):
            try:
                x = ontology.get_ontology_label(uniq.iloc[i,0])
            except:
                x = uniq.iloc[i,0]
            finally:
                uniq.iloc[i,1] = x
                
                
        uniq.set_index(column + suffix, inplace=True)
        
        if column in dataset.obs:
            dataset.obs.rename({column: column + '_original'}, axis=1, inplace=True)

        dataset.obs = dataset.obs.join(uniq, on=column + suffix)
        dataset.obs.loc[[not ":" in i for i in dataset.obs[column+suffix]], column+suffix] = ""
        dataset.obs[column + suffix].cat.remove_unused_categories(inplace=True)
        

    dataset.write(os.path.join(out_dir, c_file), compression='gzip')

