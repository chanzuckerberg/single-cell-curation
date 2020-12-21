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
    barcodes_file = sys.argv[3]
    features_file = sys.argv[4]
    cell_metadata_file = sys.argv[5]
    out_prefix = sys.argv[6]
    
    #mtx_file = "../data/original/10X_nuclei_v3_Broad/matrix.mtx.gz"
    #sample_metadata_file = "../data/original/10X_nuclei_v3_Broad/sample_metadata.csv"
    #barcodes_file = "../data/original/10X_nuclei_v3_Broad/barcode.tsv"
    #features_file = "../data/original/10X_nuclei_v3_Broad/features.tsv.gz"
    #cell_metadata_file = "../data/original/paper_tables/Table_S4_UnimodalClusters/Table_S4e_snRNA_10X_v3_B_metadata.tsv.gz"
    #out_prefix = "../data/transformed/10X_nuclei_v2_Broad"
    #accessions_file = "../data/original/paper_tables/Table_S3_accession.xlsx"
    
    # Read metadata
    print('Reading metadata')
    metadata = pd.read_csv(sample_metadata_file, index_col=1)
    barcodes = pd.read_csv(barcodes_file, index_col=1, names=['row'], header=0)
    cell_metadata = pd.read_csv(cell_metadata_file, sep="\t", index_col=0)
    cell_metadata.drop(['nUMI', 'nGene', 'QC'], inplace=True, axis=1)
    # accessions = pd.read_excel(accessions_file, sheet_name=None)
    
    metadata = metadata.join(barcodes)
    metadata = metadata.join(cell_metadata)
    metadata = metadata.sort_values(by='row')
    
    # Read gene info
    features = pd.read_csv(features_file, sep="\t", header=None)
    var_names = anndata.utils.make_index_unique(pd.Index(features[1].values))
    
    # Read mtx file as anndata
    print('Reading count matrix')
    dataset = sc.read_mtx(mtx_file)
    print('Transposing matrix')
    dataset = dataset.T

    # Append metadata
    print('Appending metadata')
    dataset.obs = metadata
    dataset.obs_names = list(metadata.index)
    
    dataset.var_names = var_names
    dataset.var['gene_ids'] = features[0].values
    dataset.var['feature_types'] = features[2].values
    
    # Write
    print('Writing to Disk')
    dataset.write(out_prefix + ".h5ad")
    
    # Write without cells that have no cell type label
    dataset =  dataset[dataset.obs[ pd.notnull(dataset.obs['subclass_label'])].index,]
    dataset.write(out_prefix + "_filtered.h5ad")
        
if __name__ == "__main__":
    main()
