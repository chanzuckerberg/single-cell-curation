import scanpy as sc
import numpy as np
import sys

in_file =  sys.argv[1]
out_file = sys.argv[2]

dataset = sc.read(in_file)
dataset.X = np.nan_to_num(dataset.X, copy=True, nan=0.0, posinf=0.0, neginf=0.0)

dataset.obsm["X_tsne"] = dataset.obs[['l1-tsne_0', 'l1-tsne_1']].to_numpy()
dataset.obsm["X_umap"] = dataset.obs[['l1-umap_0', 'l1-umap_1']].to_numpy()
#dataset.obsm["X_tsne_2"]  = dataset.obs[['l2-tsne_0', 'l2-tsne_1']].to_numpy()
#dataset.obsm["X_umap_2"] = dataset.obs[['l2-umap_0', 'l2-umap_1']].to_numpy()
#dataset.obsm["X_tsne_3"]  = dataset.obs[['l3-tsne_0', 'l3-tsne_1']].to_numpy()
#dataset.obsm["X_umap_3"] = dataset.obs[['l3-umap_0', 'l3-umap_1']].to_numpy()
#dataset.obsm["X_tsne_4"] = dataset.obs[['class_tsne_0', 'class_tsne_1']].to_numpy()
#dataset.obsm["X_umap_4"] = dataset.obs[['class_umap_0', 'class_umap_0']].to_numpy()

dataset.obs.drop(inplace=True, columns=['l1-tsne_0', 'l1-tsne_1', 'l1-umap_0', 'l1-umap_1', 'l2-tsne_0', 'l2-tsne_1', 'l2-umap_0', 'l2-umap_1', 'l3-tsne_0', 'l3-tsne_1', 'l3-umap_0', 'l3-umap_1' ,'class_tsne_0', 'class_tsne_1', 'class_umap_0', 'class_umap_0'])

dataset.write(out_file, compression='gzip')
