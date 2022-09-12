import anndata
import pandas as pd
import sys

dataset_patch_file = sys.argv[1]
dataset_10x_file = sys.argv[2]
ontology_10x_file = sys.argv[3]
out_file = sys.argv[4]

dataset_patch = anndata.read(dataset_patch_file)
dataset_10x = anndata.read(dataset_10x_file)

ontology = {}
with open(ontology_10x_file, 'r') as i:
    for line in i:
        line = line.strip().split('\t')
        ontology[line[0]] = line[1].replace(" (BICCN)", "")

# arraging 10x
to_keep = ['sex', 'cell_type', 'cell_type_ontology_term_id', 'tissue', 'tissue_ontology_term_id', 'assay', 'assay_ontology_term_id', 'ethnicity', 'ethnicity_term_ontology_id', 'disease', 'disease_ontology_term_id', 'BICCN_cluster_label', 'BICCN_subclass_label', 'BICCN_ontology_term_id', 'development_stage', 'development_stage_ontology_term_id']

dataset_10x.obs = dataset_10x.obs[to_keep]
dataset_10x.obs.loc[:,'BICCN_subclass_label'] = [ontology[i] for i in dataset_10x.obs['BICCN_ontology_term_id']]

# merging
dataset = dataset_patch.concatenate(dataset_10x)
dataset.obs.drop(columns='batch', inplace=True)
dataset.uns = dataset_patch.uns

dataset.write(out_file, compression='gzip')

