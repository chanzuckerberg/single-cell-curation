import scanpy as sc
import pandas as pd 
import sys

dataset_file = sys.argv[1]
outfile = sys.argv[2]

dataset = sc.read(dataset_file, 'r')

class_df = dataset.obs[['cell_type', 'cell_type', 'cell_type', 'cell_type']]
class_df.drop_duplicates(inplace=True)
class_df['final'] = True

class_df.to_csv(outfile, sep='\t', index=False)
