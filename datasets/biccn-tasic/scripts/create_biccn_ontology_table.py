import scanpy as sc
import pandas as pd 
import os, sys

dataset_file = sys.argv[1]
outfile = sys.argv[2]

dataset = sc.read(dataset_file, 'r')

class_df = dataset.obs[['sub_class']].astype(str)
class_df.drop_duplicates(inplace=True)

class_df.sort_values(by=['sub_class'], inplace=True)

class_df['original_columns'] = 'sub_class'
class_df['original_values'] = class_df['sub_class']
class_df['BICCN_ontology_term_id'] = "unknown"
class_df['ontology_label'] = ''
class_df['project_name'] = 'tassic'
class_df['doi'] = 'https://doi.org/10.1038/nn.4216'

# Putting in ontology terms
class_df.loc[class_df['sub_class'] == 'Astrocyte', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770141', 'Astrocyte (BICCN)']
class_df.loc[class_df['sub_class'] == 'Chodl', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770152', 'Sst neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'Endothelial', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770142', 'Endothelial cell (BICCN)']
class_df.loc[class_df['sub_class'] == 'L2/3', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770156', 'L2/3 IT neuron']
class_df.loc[class_df['sub_class'] == 'L2/3_Syt10', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770156', 'L2/3 IT neuron']
class_df.loc[class_df['sub_class'] == 'L4', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770174', 'L4/5 IT (BICCN)']
class_df.loc[class_df['sub_class'] == 'L5', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770157', 'L5 IT neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'L5_Chrna6', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770157', 'L5 IT neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'L5a', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770157', 'L5 IT neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'L5a1', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770157', 'L5 IT neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'L5a2', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770157', 'L5 IT neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'L5b', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770157', 'L5 IT neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'L5b1', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770157', 'L5 IT neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'L5b2', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770157', 'L5 IT neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'L6a', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770158', 'L6 IT neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'L6a1', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770158', 'L6 IT neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'L6a2', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770158', 'L6 IT neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'L6b', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770163', 'L6b neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'Microglia', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770146', 'Microglia (BICCN)']
class_df.loc[class_df['sub_class'] == 'Ndnf', ['BICCN_ontology_term_id', 'ontology_label']] = ['Ndnf', '']
class_df.loc[class_df['sub_class'] == 'OPC', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770139', 'Oligodendrocyte progenitor cell (BICCN)']
class_df.loc[class_df['sub_class'] == 'Oligodendrocyte', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770140', 'Oligodendrocyte (BICCN)']
class_df.loc[class_df['sub_class'] == 'Pvalb', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770154', 'Pvalb neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'Sst', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770152', 'Sst neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'Sst_Cbln4', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770152', 'Sst neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'Sst_Chodl', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770152', 'Sst neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'Sst_Clstn2', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770152', 'Sst neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'Sst_Nmbr', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770152', 'Sst neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'Sst_Nr2f2', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770152', 'Sst neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'Th', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770152', 'Sst neuron (BICCN)']
class_df.loc[class_df['sub_class'] == 'Unknown', ['BICCN_ontology_term_id', 'ontology_label']] = ['unknown', '']
class_df.loc[class_df['sub_class'] == 'Vip', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770151', 'Vip GABAergic neuron (BICCN 2020)']
class_df.loc[class_df['sub_class'] == 'nan', ['BICCN_ontology_term_id', 'ontology_label']] = ['unknown', '']

class_df.drop(['sub_class'], axis=1, inplace=True)


class_df.to_csv(outfile, sep='\t', index=False)
