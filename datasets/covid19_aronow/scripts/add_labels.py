import scanpy as sc
import pandas as pd
import ontology 

dataset = sc.read('./data/original/PBMC_merged_normalized_addRaw_coordinatesFixed_0126.h5ad')
pca = pd.read_csv('./data/original/pca_coordinates.txt', sep='\t', index_col=0)
tsne = pd.read_csv('./data/original/tsne_coordinates.txt', sep='\t', index_col=0)
umap = pd.read_csv('./data/original/umap_coordinates.txt', sep='\t', index_col=0)

# Reorder
dataset=dataset[pca.index]
tsne=tsne.loc[pca.index,]
umap=umap.loc[pca.index,]

dataset.obsm['X_pca'] = pca.to_numpy()
dataset.obsm['X_tsne'] = tsne.to_numpy()
dataset.obsm['X_umap'] = umap.to_numpy()

columns = ['tissue', 'disease', 'cell_type', 'ethnicity', 'development_stage']
suffix = '_ontology_term_id'

# Correct issues
#dataset.obs['sex'].cat.add_categories(["male", "unknown"], inplace=True)
#dataset.obs.loc[dataset.obs['sex'] == 'na', 'sex'] = 'unknown'
#dataset.obs.loc[dataset.obs['sex'] == 'M', 'sex'] = 'male'
#dataset.obs['sex'].cat.remove_categories(["M", "na"], inplace=True)

#dataset.uns['layer_descriptions']['raw.X'] = dataset.uns['layer_descriptions']['raw']
#del dataset.uns['layer_descriptions']['raw'] 
dataset.uns['version'] = {'corpora_schema_version': '1.1.0', 'corpora_encoding_version': '0.1.0'}
dataset.uns['organism_ontology_term_id'] = "NCBITaxon:9606"

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
            #uniq.iloc[i,0] = ""
            ##Rename ontology id to ""
            #dataset.obs.loc[dataset.obs[column + suffix] == x, column + suffix] = ""
            #dataset.obs[column + suffix].cat.remove_categories(x, inplace=True)
        finally:
            uniq.iloc[i,1] = x
            
            
    uniq.set_index(column + suffix, inplace=True)
    
    if column in dataset.obs:
        dataset.obs.rename({column: column + '_original'}, axis=1, inplace=True)

    dataset.obs = dataset.obs.join(uniq, on=column + suffix)
    dataset.obs.loc[[not ":" in i for i in dataset.obs[column+suffix]], column+suffix] = ""
    dataset.obs[column + suffix].cat.remove_unused_categories(inplace=True)
    

dataset.write('./data/remixed/PBMC_merged_normalized_addRaw_coordinatesFixed_0126.h5ad', compression='gzip')

#for column in columns:
#    uniq = dataset.obs[[column + suffix, column]].drop_duplicates()
#    print(uniq)
