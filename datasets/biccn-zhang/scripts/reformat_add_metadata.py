# Reformats data and appends metadata. Outputs h5ad

import pandas as pd
import scanpy as sc
import sys

def main():
    
    # Files
    mtx_file_counts = sys.argv[1]
    sample_metadata_file = sys.argv[2]
    umap_file = sys.argv[3]
    out_file = sys.argv[4]
    
    #mtx_file_counts = "./data/original/counts.h5ad"
    #sample_metadata_file = "./data/original/cell_metadata.csv"
    #umap_file = "./data/original/umap_embedding.csv"
    #out_file = "./data/transformed/1_reformatted/tassic.h5ad"
    
    # Read metadata
    print('Reading metadata')
    metadata = pd.read_csv(sample_metadata_file, index_col=0)
    umap = pd.read_csv(umap_file, index_col=0)
    
    # Read counts
    dataset = sc.read(mtx_file_counts)

    # Append metadata
    print('Appending metadata')
    dataset=dataset[umap.index,:]
    dataset.obsm['X_umap'] = umap.to_numpy()
    dataset.obs = metadata.loc[dataset.obs_names]
    dataset.raw=dataset
    
    # Write
    dataset.write(out_file, compression='gzip')
    
        
if __name__ == "__main__":
    main()
