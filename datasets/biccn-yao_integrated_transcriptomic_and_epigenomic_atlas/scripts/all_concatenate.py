import scanpy as sc
import anndata
import sys
import os
import re

working_dir = sys.argv[1]
out_file = sys.argv[2]

all_files = [os.path.join(working_dir, i) for i in os.listdir(working_dir)]

datasets = []
for i in range(0, len(all_files)):
    datasets.append(sc.read(all_files[i]))
    experiment = re.sub('_filtered.h5ad$', '', os.path.basename(all_files[i]))
    datasets[i].obs['BICCN_project'] = experiment
    

joined = anndata.concat(datasets, axis=0)
joined.uns = datasets[0].uns
joined.uns['title'] = 'An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types'

joined.write(out_file)
