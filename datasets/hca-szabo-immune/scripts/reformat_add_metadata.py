# Reformats data and appends metadata. Outputs h5ad

import pandas as pd
import scanpy as sc
import sys

def main():
    
    # Files
    mtx_file_counts = "./data/original/counts.tsv.gz"
    sample_metadata_file = "./data/original/41467_2019_12464_MOESM9_ESM.txt.gz"
    out_file = "./data/transformed/1_reformatted/hca_immune_szabo.h5ad"
    
    # Read metadata
    print('Reading metadata')
    metadata = pd.read_csv(sample_metadata_file, sep='\t')
    
    # Fix tissues
    metadata.loc[metadata.loc[:,'tissue'] == 'LG','tissue'] = 'lung'
    metadata.loc[metadata.loc[:,'tissue'] == 'BM','tissue'] = 'bone marrow'
    metadata.loc[metadata.loc[:,'tissue'] == 'LN','tissue'] = 'lymph node'
    metadata.loc[metadata.loc[:,'tissue'] == 'BL','tissue'] = 'blood'
    
    # Make index
    metadata['index'] = metadata.loc[:,'barcode'] + '_' + metadata.loc[:,'donor'] + '_' + metadata.loc[:,'tissue'] + '_' + metadata.loc[:,'stimulation_status']
    metadata.set_index('index', inplace=True)
    metadata['status'] = metadata['stimulation_status'].astype('string') + '_' + metadata['cd4cd8_status']
    
    umap = metadata[['umap_x', 'umap_y']].to_numpy()
    metadata.drop(['umap_x', 'umap_y'], axis=1, inplace=True)
    
    # Read mtx file as anndata
    print('Reading count matrix')
    dataset_counts = sc.read_text(mtx_file_counts, delimiter='\t', first_column_names=True)
    dataset_counts = dataset_counts.T
    
    # Reshaping
    print('Appending metadata')
    dataset_counts = dataset_counts[metadata.index,]
    dataset_counts.obs = metadata
    dataset_counts.obsm['X_umap'] = umap
    
    dataset_counts.raw = dataset_counts
    dataset_counts.var_names_make_unique()

    # Write
    dataset_counts.write(out_file, compression='gzip')
    
        
if __name__ == "__main__":
    main()
