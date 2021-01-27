import scanpy as sc
import pandas as pd 
import sys

dataset_file = sys.argv[1]
outfile = sys.argv[2]

dataset = sc.read(dataset_file, 'r')

class_df = dataset.obs[['brain_region']]
class_df.drop_duplicates(inplace=True)
class_df['lookup_name'] = class_df['brain_region']
class_df['ontology_term_id'] = ''
class_df['ontology_term_name'] = ''
class_df['final'] = True

class_df.loc[class_df['brain_region'] == 'VISp', 'ontology_term_id'] = 'UBERON:0002436'
class_df.loc[class_df['brain_region'] == 'ALM', 'ontology_term_id'] = 'UBERON:0001870'

class_df.to_csv(outfile, sep='\t', index=False)
