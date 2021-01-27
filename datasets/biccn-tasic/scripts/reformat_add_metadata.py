# Reformats data and appends metadata. Outputs h5ad

import pandas as pd
import scanpy as sc
import sys

def main():
    
    # Files
    mtx_file_counts = sys.argv[1]
    mtx_file_rpkm = sys.argv[2]
    sample_metadata_file = sys.argv[3]
    out_file = sys.argv[4]
    
    #mtx_file_counts = "./data/original/genes_counts.csv"
    #mtx_file_rpkm = "./data/original/genes_rpkm.csv"
    #sample_metadata_file = "./data/original/cell_metadata.csv"
    #out_file = "./data/transformed/1_reformatted/tassic.h5ad"
    
    # Read metadata
    print('Reading metadata')
    metadata = pd.read_csv(sample_metadata_file, index_col=0)
    
    # Read mtx file as anndata
    print('Reading count matrix')
    dataset_counts = sc.read_csv(mtx_file_counts)
    print('Reading rpkm matrix')
    dataset_rpkm = sc.read_csv(mtx_file_rpkm)
    
    print('Transposing matrices')
    dataset_counts = dataset_counts.T
    dataset_rpkm = dataset_rpkm.T

    # Append metadata
    print('Appending metadata')
    dataset_rpkm.obs = metadata.loc[dataset_counts.obs_names]
    
    # Merging
    dataset_counts = dataset_counts[dataset_counts.obs_names, dataset_counts.var_names]
    dataset_rpkm.raw = dataset_counts
    
    # Write
    dataset_rpkm.write(out_file, compression='gzip')
    
        
if __name__ == "__main__":
    main()
