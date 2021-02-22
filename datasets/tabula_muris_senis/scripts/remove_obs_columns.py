import sys
import scanpy as sc

cols = sys.argv[1].split(',')
dataset_file = sys.argv[2]

dataset = sc.read(dataset_file)

for i in cols:
    del dataset.obs[i]
    
dataset.write(dataset_file, compression='gzip')
