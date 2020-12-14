# that need to be manually curated

import scanpy as sc 
import sys 

h5ad_in = sys.argv[1]
h5ad_out = sys.argv[2]

dataset = sc.read(h5ad_in)
dataset.obs['health_status'] = 'healthy'
dataset.obs.loc[ lambda d: d['Fetus_id'] == 'H27432', 'health_status'] = 'trisomy_18'

dataset.write(h5ad_out)
