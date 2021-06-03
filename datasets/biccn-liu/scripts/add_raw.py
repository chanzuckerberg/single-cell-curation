import scanpy as sc
import pandas as pd
import sys

in_file =  sys.argv[1]
out_file = sys.argv[2]

dataset = sc.read(in_file)
dataset.raw = dataset

dataset.write(out_file, compression='gzip')
