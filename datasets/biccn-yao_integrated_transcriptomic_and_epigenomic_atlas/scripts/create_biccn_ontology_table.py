import scanpy as sc
import pandas as pd 
import os, sys

working_dir = sys.argv[1]
outfile = sys.argv[2]

dataset_files = [os.path.join( working_dir, i) for i in os.listdir(working_dir)]
datasets = [sc.read(i, 'r') for i in dataset_files]

class_df = pd.concat(axis=0, objs=[ i.obs[[ 'BICCN_class_label', 'BICCN_subclass_label']].drop_duplicates() for i in datasets ])
class_df.drop_duplicates(inplace=True)
class_df.sort_values(by=['BICCN_class_label', 'BICCN_subclass_label'], inplace=True)
class_df[['BICCN_ontology_term_id']] = class_df[['BICCN_subclass_label']]


# Putting in ontology terms
class_df.loc[class_df['BICCN_subclass_label'] == 'Lamp5', ['BICCN_ontology_term_id']] = ['ILX:0770149']
class_df.loc[class_df['BICCN_subclass_label'] == 'Pvalb', ['BICCN_ontology_term_id']] = ['ILX:0770154']
class_df.loc[class_df['BICCN_subclass_label'] == 'Sncg', ['BICCN_ontology_term_id']] = ['ILX:0770150']
class_df.loc[class_df['BICCN_subclass_label'] == 'Sst', ['BICCN_ontology_term_id']] = ['ILX:0770152']
class_df.loc[class_df['BICCN_subclass_label'] == 'Vip', ['BICCN_ontology_term_id']] = ['ILX:0770151']
#class_df.loc[class_df['BICCN_subclass_label'] == 'CR', ['BICCN_ontology_term_id']] = ['']
class_df.loc[class_df['BICCN_subclass_label'] == 'L2/3 IT', ['BICCN_ontology_term_id']] = ['ILX:0770156']
class_df.loc[class_df['BICCN_subclass_label'] == 'L5 ET', ['BICCN_ontology_term_id']] = ['ILX:0770160']
class_df.loc[class_df['BICCN_subclass_label'] == 'L5 IT', ['BICCN_ontology_term_id']] = ['ILX:0770157']
class_df.loc[class_df['BICCN_subclass_label'] == 'L5/6 NP', ['BICCN_ontology_term_id']] = ['ILX:0770161']
class_df.loc[class_df['BICCN_subclass_label'] == 'L6 CT', ['BICCN_ontology_term_id']] = ['ILX:0770162']
class_df.loc[class_df['BICCN_subclass_label'] == 'L6 IT', ['BICCN_ontology_term_id']] = ['ILX:0770158']
class_df.loc[class_df['BICCN_subclass_label'] == 'L6 IT Car3', ['BICCN_ontology_term_id']] = ['ILX:0770159']
class_df.loc[class_df['BICCN_subclass_label'] == 'L6b', ['BICCN_ontology_term_id']] = ['ILX:0770163']
class_df.loc[class_df['BICCN_subclass_label'] == 'Astro', ['BICCN_ontology_term_id']] = ['ILX:0770141']
class_df.loc[class_df['BICCN_subclass_label'] == 'Endo', ['BICCN_ontology_term_id']] = ['ILX:0770142']
#class_df.loc[class_df['BICCN_subclass_label'] == 'Macrophage', ['BICCN_ontology_term_id']] = ['']
class_df.loc[class_df['BICCN_subclass_label'] == 'Micro', ['BICCN_ontology_term_id']] = ['ILX:0770146']
class_df.loc[class_df['BICCN_subclass_label'] == 'OPC', ['BICCN_ontology_term_id']] = ['ILX:0770139']
class_df.loc[class_df['BICCN_subclass_label'] == 'Oligo', ['BICCN_ontology_term_id']] = ['ILX:0770140']
class_df.loc[class_df['BICCN_subclass_label'] == 'Peri', ['BICCN_ontology_term_id']] = ['ILX:0770145']
#class_df.loc[class_df['BICCN_subclass_label'] == 'Prog/IP', ['BICCN_ontology_term_id']] = ['']
class_df.loc[class_df['BICCN_subclass_label'] == 'SMC', ['BICCN_ontology_term_id']] = ['ILX:0770144']
class_df.loc[class_df['BICCN_subclass_label'] == 'VLMC', ['BICCN_ontology_term_id']] = ['ILX:0770143']


class_df.to_csv(outfile, sep='\t', index=False)
