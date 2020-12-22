import scanpy as sc
import pandas as pd
import sys

in_file =  sys.argv[1]
org_file = sys.argv[2]
out_file = sys.argv[3]

org = pd.read_csv(org_file, sep='\t', index_col=0, header=0)
org=org.loc[:,['organism_ontology_term_id', 'organism']]

dataset = sc.read(in_file)
dataset.obs = dataset.obs.join(org, on='orig.ident')
dataset.write(out_file)
