import scanpy as sc
import pandas as pd
import sys

in_file =  sys.argv[1]
ontology_file = sys.argv[2]
out_file = sys.argv[3]

ontology = pd.read_csv(ontology_file, sep='\t', index_col=['BICCN_subclass_label', 'BICCN_class_label'], header=0)
ontology = ontology['BICCN_ontology_term_id']

dataset = sc.read(in_file)
dataset.obs = dataset.obs.join(ontology, on=['BICCN_subclass_label', 'BICCN_class_label'])

dataset.raw = dataset

dataset.write(out_file)
