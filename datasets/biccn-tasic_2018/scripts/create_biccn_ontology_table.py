import scanpy as sc
import pandas as pd 
import os, sys

dataset_file = sys.argv[1]
outfile = sys.argv[2]

dataset = sc.read(dataset_file, 'r')

class_df = dataset.obs[['subclass']].astype(str)
class_df.drop_duplicates(inplace=True)

class_df.sort_values(by=['subclass'], inplace=True)

class_df['original_columns'] = 'subclass'
class_df['original_values'] = class_df['subclass']
class_df['BICCN_ontology_term_id'] = "unknown"
class_df['ontology_label'] = ''
class_df['project_name'] = 'tassic_2018'
class_df['doi'] = 'https://doi.org/10.1038/s41586-018-0654-5'

# Putting in ontology terms
class_df.loc[class_df['subclass'] == 'L5 IT', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770157', 'L5 IT neuron (BICCN)']
class_df.loc[class_df['subclass'] == 'L5 PT', ['BICCN_ontology_term_id', 'ontology_label']] = ['L5 PT', 'L5 PT']
class_df.loc[class_df['subclass'] == 'CR', ['BICCN_ontology_term_id', 'ontology_label']] = ['CR', 'CR']
class_df.loc[class_df['subclass'] == 'L6 CT', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770162', 'L6 CT neuron (BICCN)']
class_df.loc[class_df['subclass'] == 'L6b', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770163', 'L6b neuron (BICCN)']
class_df.loc[class_df['subclass'] == 'L4', ['BICCN_ontology_term_id', 'ontology_label']] = ['L4', 'L4']
class_df.loc[class_df['subclass'] == 'L6 IT', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770158', 'L6 IT (BICCN)']
class_df.loc[class_df['subclass'] == 'NP', ['BICCN_ontology_term_id', 'ontology_label']] = ['NP', 'NP']
class_df.loc[class_df['subclass'] == 'L2/3 IT', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770156', 'L2/3 IT (BICCN)']
class_df.loc[class_df['subclass'] == 'Sncg', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770150', 'Sncg (BICCN)']
class_df.loc[class_df['subclass'] == 'Serpinf1', ['BICCN_ontology_term_id', 'ontology_label']] = ['Serpinf1', 'Serpinf1']
class_df.loc[class_df['subclass'] == 'Meis2', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770155', 'Meis2 neuron (BICCN)']
class_df.loc[class_df['subclass'] == 'Sst', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770152', 'Sst neuron (BICCN)']
class_df.loc[class_df['subclass'] == 'Pvalb', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770154', 'Pvalb neuron (BICCN)']
class_df.loc[class_df['subclass'] == 'Vip', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770151', 'Vip GABAergic neuron (BICCN)']
class_df.loc[class_df['subclass'] == 'Lamp5', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770149', 'Lamp5 neuron (BICCN)']
class_df.loc[class_df['subclass'] == 'SMC', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770144', 'Smooth muscle cell (BICCN)']
class_df.loc[class_df['subclass'] == 'Peri', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770145', 'Pericyte (BICCN)']
class_df.loc[class_df['subclass'] == 'Endo', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770142', 'Endothelial cell (BICCN)']
class_df.loc[class_df['subclass'] == 'VLMC', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770143', 'Vascular leptomeningeal cell (BICCN)']
class_df.loc[class_df['subclass'] == 'Macrophage', ['BICCN_ontology_term_id', 'ontology_label']] = ['Macrophage', 'Macrophage']
class_df.loc[class_df['subclass'] == 'Astro', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770141', 'Astrocyte (BICCN)']
class_df.loc[class_df['subclass'] == 'Oligo', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770140', 'Oligodendrocyte']
class_df.loc[class_df['subclass'] == 'No Class', ['BICCN_ontology_term_id', 'ontology_label']] = ['unknown', 'unknown']

class_df.drop(['subclass'], axis=1, inplace=True)


class_df.to_csv(outfile, sep='\t', index=False)
