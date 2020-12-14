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
    #mtx_file = sys.argv[1]
    #sample_metadata_file = sys.argv[2]
    #features_file = sys.argv[3]
    #clusters_file = sys.argv[4]
    #out_file = sys.argv[5]
    
    working_dir = "../data/original/"
    
    mtx_file = working_dir + "gene_sorted-raw_matrix.mtx.gz"
    sample_metadata_file = working_dir + "cerebellum_metadata_updated_custom.tsv"
    features_file = working_dir + "genes_cerebellum.tsv"
    clusters_file = working_dir + "cerebellum_clusters.tsv"
    out_file = working_dir + "kozareva_transcriptomic_atlas_of_the_mouse_cerebellum.h5ad"
    
    # Read metadata
    print('Reading metadata')
    metadata = pd.read_csv(sample_metadata_file, index_col=0, sep="\t")
    metadata = metadata.iloc[1:, :]
    
    clusters = pd.read_csv(clusters_file, index_col=0, sep="\t")
    clusters = clusters.astype(float)
    clusters = clusters.loc[metadata.index, :]
    
    # Read gene info
    features = pd.read_csv(features_file, sep="\t", header=None)
    var_names = anndata.utils.make_index_unique(pd.Index(features[0].values))
    
    # Read mtx file as anndata
    print('Reading count matrix')
    dataset = sc.read_mtx(mtx_file)
    print('Transposing matrix')
    dataset = dataset.T

    # Append metadata
    print('Appending metadata')
    dataset.obs = metadata
    dataset.obs_names = list(metadata.index)
    
    dataset.obsm['X_tsne'] = clusters.to_numpy()
    
    dataset.var_names = var_names
    
    # Write
    print('Writing to Disk')
    dataset.write(out_file)
        
if __name__ == "__main__":
    main()
