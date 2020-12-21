import scanpy as sc
import pandas as pd
import sys

dataset_file = sys.argv[1]
umap_file = sys.argv[2]
out_file = sys.argv[3]

dataset = sc.read(dataset_file)
umap = pd.read_csv(umap_file, header=None, index_col=1)

# Create raw layer
dataset.raw = dataset

# Normalize
sc.pp.normalize_total(dataset)
sc.pp.log1p(dataset)

# Append umap
dataset.obsm['X_umap'] = umap.loc[dataset.obs_names, [2,3]].to_numpy()

dataset.write(out_file)
