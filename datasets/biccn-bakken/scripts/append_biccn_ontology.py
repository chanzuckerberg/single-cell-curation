import scanpy as sc
import pandas as pd
import sys

in_file =  sys.argv[1]
ontology_file = sys.argv[2]
out_file = sys.argv[3]

ontology = pd.read_csv(ontology_file, sep='\t', index_col=1, header=0)
original_columns=ontology['original_columns'].values[0].split('@')
ontology=ontology.loc[:,'BICCN_ontology_term_id']

dataset = sc.read(in_file)
dataset.obs = dataset.obs.join(ontology, on=original_columns)

dataset.write(out_file)
