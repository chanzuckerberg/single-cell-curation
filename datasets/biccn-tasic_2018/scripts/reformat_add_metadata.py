# Reformats data and appends metadata. Outputs h5ad

import pandas as pd
import scanpy as sc
import sys

def main():
    
    # Files
    mtx_file_counts = sys.argv[1]
    sample_metadata_file = sys.argv[2]
    tsne_file = sys.argv[3]
    out_file = sys.argv[4]
    
    #mtx_file_counts = "./data/original/GSE115746_cells_exon_counts.csv.gz"
    #sample_metadata_file = "./data/original/Supplementary_Table_10_Full_Metadata2.txt"
    #tsne_file = "./data/original/tsne_tassic_2018_no_NA.csv"
    #out_file = "./data/transformed/1_reformatted/tassic_2018.h5ad"
    
    # Read metadata
    print('Reading metadata')
    metadata = pd.read_csv(sample_metadata_file, index_col=0, sep='\t', encoding="ISO-8859-1")
    tsne = pd.read_csv(tsne_file, index_col=1)
    tsne = tsne.iloc[:,1:3]
    
    # Read mtx file as anndata
    print('Reading count matrix')
    dataset_counts = sc.read_csv(mtx_file_counts)
    dataset_counts = dataset_counts.T
    
    # Reshaping
    cell_ids = set(metadata.index)
    cell_ids = list(cell_ids.intersection(set(dataset_counts.obs_names)))
    
    metadata = metadata.reindex(index=cell_ids)
    dataset_counts = dataset_counts[cell_ids,]
    tsne = tsne.reindex(index=list(metadata['seq_name']))

    # Appending 
    print('Appending metadata')
    dataset_counts.obs = metadata
    dataset_counts.obsm['X_tsne'] = tsne.to_numpy()
    
    dataset_counts = dataset_counts[list(dataset_counts.obs['class'] != 'Low Quality'),]
    dataset_counts.raw=dataset_counts
    
    # Write
    dataset_counts.write(out_file, compression='gzip')
    
        
if __name__ == "__main__":
    main()
