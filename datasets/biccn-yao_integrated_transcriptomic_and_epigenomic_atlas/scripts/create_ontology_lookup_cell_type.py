import scanpy as sc
import pandas as pd 
import os, sys

working_dir = sys.argv[1]
outfile = sys.argv[2]

dataset_files = [os.path.join( working_dir, i) for i in os.listdir(working_dir)]
datasets = [sc.read(i, 'r') for i in dataset_files]

class_df = pd.concat(axis=0, objs=[ i.obs[[ 'temp_class_label', 'BICCN_class_label', 'BICCN_subclass_label']].drop_duplicates() for i in datasets ])
class_df.drop_duplicates(inplace=True)
class_df.sort_values(by=['BICCN_class_label', 'BICCN_subclass_label'], inplace=True)
class_df[['ontology_term_id', 'ontology_term_name']] = 'NA'
class_df[['final']] = True


# Putting in ontology terms
class_df.loc[class_df['BICCN_subclass_label'] == 'Astro', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000127', 'astrocyte']
class_df.loc[class_df['BICCN_subclass_label'] == 'Endo', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000115', 'endothelial cell']
class_df.loc[class_df['BICCN_subclass_label'] == 'Macrophage', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000235', 'macrophage']
class_df.loc[class_df['BICCN_subclass_label'] == 'Micro', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000129', 'microgial cell']
class_df.loc[class_df['BICCN_subclass_label'] == 'OPC', ['ontology_term_id', 'ontology_term_name']] = ['CL:0002453', 'oligodendrocyte precursor cell']
class_df.loc[class_df['BICCN_subclass_label'] == 'Oligo', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000128', 'oligodendrocyte']
class_df.loc[class_df['BICCN_subclass_label'] == 'Peri', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000669', 'pericyte cell']
class_df.loc[class_df['BICCN_subclass_label'] == 'Prog/IP', ['ontology_term_id', 'ontology_term_name']] = ['prog', 'prog']
class_df.loc[class_df['BICCN_subclass_label'] == 'SMC', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000192', 'smooth muscle cell']
class_df.loc[class_df['BICCN_subclass_label'] == 'VLMC', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000708', 'leptomeningeal cel']

class_df.loc[class_df['BICCN_class_label'] == 'GABAergic', ['ontology_term_id', 'ontology_term_name']] = ["CL:0000617", 'GABAergic neuron']
class_df.loc[class_df['BICCN_class_label'] == 'Glutamatergic', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000679', 'glutamatergic neuron']

class_df.drop('BICCN_class_label', axis=1, inplace=True)

class_df.to_csv(outfile, sep='\t', index=False)
