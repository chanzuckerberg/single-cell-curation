import scanpy as sc
import pandas as pd 
import os, sys

dataset_file = sys.argv[1]
outfile = sys.argv[2]

dataset = sc.read(dataset_file, 'r')

class_df = dataset.obs[['cell_type__custom', 'cell_type__ontology_label']].astype(str)
class_df.drop_duplicates(inplace=True)
class_df.sort_values(by=['cell_type__custom', 'cell_type__ontology_label'], inplace=True)

class_df['original_columns'] = 'cell_type__custom@cell_type__ontology_label'
class_df['original_values'] = class_df['cell_type__custom'] + "@" +  class_df['cell_type__ontology_label']
class_df['BICCN_ontology_term_id'] = class_df['cell_type__custom'].astype(str)
class_df['ontology_label'] = ''
class_df['project_name'] = 'kozareva'
class_df['doi'] = 'https://doi.org/10.1101/2020.03.04.976407'

# Putting in ontology terms
class_df.loc[class_df['cell_type__custom'] == 'astrocyte of the cerebellum', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770141', 'Astrocyte (BICCN)']
class_df.loc[class_df['cell_type__custom'] == 'oligodendrocyte', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770140', 'Oligodendrocyte (BICCN)']
class_df.loc[class_df['cell_type__custom'] == 'oligodendrocyte precursor cell', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770139', 'Oligodendrocyte progenitor cell (BICCN)']
class_df.loc[class_df['cell_type__custom'] == 'Endothelial mural', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770142', 'Endothelial cell (BICCN)']
class_df.loc[class_df['cell_type__custom'] == 'Endothelial mural', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770142', 'Endothelial cell (BICCN)']
class_df.loc[class_df['cell_type__custom'] == 'endothelial stalk cell', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770142', 'Endothelial cell (BICCN)']
class_df.loc[class_df['cell_type__custom'] == 'microglial cell', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770146', 'Microglia (BICCN)']
#class_df.loc[class_df['cell_type__custom'] == 'Purkinje layer interneuron (PLI)', ['BICCN_ontology_term_id', 'ontology_label']] = 
#class_df.loc[class_df['cell_type__custom'] == 'Molecular layer interneuron 1 (MLI1)', ['BICCN_ontology_term_id', 'ontology_label']] =
#class_df.loc[class_df['cell_type__custom'] == 'Molecular layer interneuron 2 (MLI2)', ['BICCN_ontology_term_id', 'ontology_label']] = 
#class_df.loc[class_df['cell_type__custom'] == 'Purkinje cell', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770152']
#class_df.loc[class_df['cell_type__custom'] == 'Choroid', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770151']
#class_df.loc[class_df['cell_type__custom'] == 'ependymal cell', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770156']
#class_df.loc[class_df['cell_type__custom'] == 'Bergmann glial cell', ['BICCN_ontology_term_id', 'ontology_label']] = ['ILX:0770160']
#class_df.loc[class_df['cell_type__custom'] == 'cerebellar Golgi cell', ['BICCN_ontology_term_id', 'ontology_label']] = 
#class_df.loc[class_df['cell_type__custom'] == 'Unipolar brush cell (UBC)', ['BICCN_ontology_term_id', 'ontology_label']] =
#class_df.loc[class_df['cell_type__custom'] == 'macrophage', ['BICCN_ontology_term_id', 'ontology_label']] =
#class_df.loc[class_df['cell_type__custom'] == 'fibroblast', ['BICCN_ontology_term_id', 'ontology_label']] = 
#class_df.loc[class_df['cell_type__custom'] == 'cerebellar granule cell', ['BICCN_ontology_term_id', 'ontology_label']] = 

class_df.drop(['cell_type__ontology_label', 'cell_type__custom'], axis=1, inplace=True)


class_df.to_csv(outfile, sep='\t', index=False)
