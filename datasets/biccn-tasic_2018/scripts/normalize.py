import scanpy as sc
import sys

dataset_file = sys.argv[1]
out_file = sys.argv[2]

dataset = sc.read(dataset_file)

# Create raw layer
dataset.raw = dataset

# Normalize
sc.pp.normalize_total(dataset)
sc.pp.log1p(dataset)

dataset.write(out_file)
