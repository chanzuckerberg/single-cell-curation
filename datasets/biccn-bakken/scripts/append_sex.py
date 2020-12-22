import scanpy as sc
import pandas as pd
import sys

in_file =  sys.argv[1]
sex_file = sys.argv[2]
out_file = sys.argv[3]

sex = pd.read_csv(sex_file, sep='\t', index_col=1, header=0)
original_columns=sex['original_columns'].values[0].split('@')
sex=sex.iloc[:,1]

dataset = sc.read(in_file)
# Correct the dataset that has a column 'donor' instead of 'donor_id'
if 'donor' in list(dataset.obs.columns):
    dataset.obs.rename({'donor': 'donor_id'}, axis=1, inplace=True)

dataset.obs = dataset.obs.join(sex, on=original_columns)
dataset.write(out_file)
