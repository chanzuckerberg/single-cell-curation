import scanpy as sc
import pandas as pd
import sys

dataset_file = sys.argv[1]
ontology_table = sys.argv[2]
out_file = sys.argv[3]

dataset = sc.read(dataset_file)

# read ontology
ontology = pd.read_csv(ontology_table, index_col=1, sep="\t")
orginal_col = ontology.iloc[0,0]
ontology = ontology["ontology_label"].str.replace(" \(BICCN\)", "")
ontology = ontology.str.replace(" \(BICCN 2020\)", "")
ontology = ontology.rename('BICCN_subclass_label')

# Conver to numerical
to_convert = ["Length (bp)", "Soma depth (µm)", "Yield (pgµl)"]

for i in to_convert:
    dataset.obs[i] = dataset.obs[i].str.replace("?", "nan")
    dataset.obs[i] = dataset.obs[i].astype(float)

# merge to obs 
dataset.obs = dataset.obs.join(ontology, on=orginal_col)

dataset.write_h5ad(out_file, compression='gzip')
