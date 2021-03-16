import scanpy as sc
import pandas as pd 
import sys

dataset_file = sys.argv[1]
outfile = sys.argv[2]

dataset = sc.read(dataset_file, 'r')


class_df = dataset.obs[['subclass']]
class_df.drop_duplicates(inplace=True)
class_df['lookup_name'] = class_df['subclass']
class_df['ontology_term_id'] = ''
class_df['ontology_term_name'] = ''
class_df['final'] = True

class_df.loc[class_df['subclass'] == 'L5 IT', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['subclass'] == 'L5 PT', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['subclass'] == 'CR', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['subclass'] == 'L6 CT', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['subclass'] == 'L6b', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['subclass'] == 'L4', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['subclass'] == 'L6 IT', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['subclass'] == 'NP', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['subclass'] == 'L2/3 IT', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['subclass'] == 'Sncg', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['subclass'] == 'Serpinf1', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['subclass'] == 'Meis2', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['subclass'] == 'Sst', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['subclass'] == 'Pvalb', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['subclass'] == 'Vip', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['subclass'] == 'Lamp5', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['subclass'] == 'SMC', 'ontology_term_id'] = 'CL:0000192'
class_df.loc[class_df['subclass'] == 'Peri', 'ontology_term_id'] = 'CL:0000669'
class_df.loc[class_df['subclass'] == 'Endo', 'ontology_term_id'] = 'CL:0000115'
class_df.loc[class_df['subclass'] == 'VLMC', 'ontology_term_id'] = 'CL:0000708'
class_df.loc[class_df['subclass'] == 'Macrophage', 'ontology_term_id'] = 'CL:0000235'
class_df.loc[class_df['subclass'] == 'Astro', 'ontology_term_id'] = 'CL:0000127'
class_df.loc[class_df['subclass'] == 'Oligo', 'ontology_term_id'] = 'CL:0000128'
class_df.loc[class_df['subclass'] == 'No Class', 'ontology_term_id'] = 'unknown'

class_df.to_csv(outfile, sep='\t', index=False)
