import scanpy as sc
import sys


file_in = sys.argv[1]
file_out = sys.argv[2]

dataset = sc.read(file_in)

tissue = list(dataset.obs['tissue'].drop_duplicates())

if len(tissue) == 1:
    print('appending: ', tissue[0])
    dataset.uns['title'] = tissue[0].capitalize() + ' — ' + dataset.uns['title']
else:
    print('appending: all')
    dataset.uns['title'] = 'All — ' + dataset.uns['title']
    
dataset.write_h5ad(file_out, compression='gzip')
