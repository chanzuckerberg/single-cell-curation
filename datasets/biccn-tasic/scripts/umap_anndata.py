import umap.umap_
import scanpy as sc
import sys

in_file =  sys.argv[1]
out_file = sys.argv[2]

dataset = sc.read(in_file)

dataset_umap = dataset.copy()
dataset_umap.X = dataset.raw.X

# Normalize
sc.pp.normalize_total(dataset_umap)
sc.pp.log1p(dataset_umap)
sc.pp.highly_variable_genes(dataset_umap)
sc.pp.pca(dataset_umap)
sc.pp.neighbors(dataset_umap)
sc.tl.tsne(dataset_umap)
sc.tl.umap(dataset_umap)


dataset.obsm = dataset_umap.obsm
dataset.obsp = dataset_umap.obsp

dataset.write(out_file, compression='gzip')
