import scanpy as sc
import pandas as pd
import sys

in_file =  sys.argv[1]
ontology_file = sys.argv[2]
out_file = sys.argv[3]

ontology = pd.read_csv(ontology_file, sep='\t', index_col=False, header=0)
ontology[['index1']] = ontology['original_values'].str.split('@', expand=True)
ontology.set_index(['index1'], inplace=True)
original_columns=ontology['original_columns'].values[0].split('@')
ontology=ontology.loc[:,'BICCN_ontology_term_id']

dataset = sc.read(in_file)
dataset.obs = dataset.obs.join(ontology, on=original_columns)

dataset.write(out_file, compression='gzip')
