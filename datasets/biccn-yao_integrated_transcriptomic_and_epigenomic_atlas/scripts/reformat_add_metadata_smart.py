# Reformats data and appends metadata. Outputs h5ad
# Steps
# 1. Transposes mtx file to make it compatible with scanpy
# 2. Reads in mtx count file as anndata, and metadata files.
# 3. Joins metadata matrices
# 4. Appends metadata to anndata
# 5. Writes to file

# example
# python3 reformat_add_metadata_10X.py ./data/original/10X_nuclei_v3_AIBS/matrix.mtx.gz ./data/original/10X_nuclei_v3_AIBS/sample_metadata.csv ./data/original/10X_nuclei_v3_AIBS/barcode.tsv /data/original/10X_nuclei_v3_AIBS/features.tsv.gz ./data/original/paper_tables/Table_S4_UnimodalClusters/Table_S4f_snRNA_10X_v3_A_metadata.tsv.gz ./data/transformed/10X_nuclei_v2_AIBS

import pandas as pd
import anndata
import scanpy as sc
import sys

def main():
    
    # Files
    mtx_file = sys.argv[1]
    sample_metadata_file = sys.argv[2]
    cell_metadata_file = sys.argv[3]
    out_prefix = sys.argv[4]
    
    #mtx_file = "../data/original/SMARTer_cells_MOp/exon.counts_small.csv.gz"
    ##sample_metadata_file = "../data/original/SMARTer_cells_MOp/sample_metadata.csv.gz"
    #cell_metadata_file = "../data/original/paper_tables/Table_S4_UnimodalClusters/Table_S4a_scRNA_SMART_metadata.tsv.gz"
    #out_prefix = "../data/transformed/SMARTer_cells_MOp"
    
    # Read metadata
    print('Reading metadata')
    #metadata = pd.read_csv(sample_metadata_file, index_col=0)
    cell_metadata = pd.read_csv(cell_metadata_file, sep="\t", index_col=0)
    metadata = pd.read_csv(cell_metadata_file, sep="\t", index_col=0)
    
    # Eliminate shared columns
    cell_metadata.drop(list(set(metadata.columns) & set(cell_metadata.columns)), inplace=True, axis=1)
    metadata = metadata.join(cell_metadata)
    
    # Read mtx file as anndata
    print('Reading count matrix')
    dataset = sc.read_csv(mtx_file)
    print('Transposing matrix')
    dataset = dataset.T
    
    # Finding shared cells
    dataset = dataset[list(set(metadata.index) & set(dataset.obs_names)),:]

    # Append metadata
    print('Appending metadata')
    dataset.obs = metadata.loc[dataset.obs_names,:]
    
    # Write
    print('Writing to Disk')
    dataset.write(out_prefix + ".h5ad")
    
    # Write without cells that have no cell type label
    dataset =  dataset[dataset.obs[ pd.notnull(dataset.obs['subclass_label'])].index,]
    dataset.write(out_prefix + "_filtered.h5ad")
        
if __name__ == "__main__":
    main()
