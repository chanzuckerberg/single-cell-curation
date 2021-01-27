import scanpy as sc
import pandas as pd
import sys

columns = sys.argv[1]
in_file =  sys.argv[2]
out_file = sys.argv[3]

columns = columns.split(',')
col_dict = {}
for i in columns:
    x = i.split(':')
    col_dict[x[0]] = x[1]

dataset = sc.read(in_file)
dataset.obs.rename(col_dict, axis=1, inplace=True)
dataset.write(out_file, compression='gzip')
