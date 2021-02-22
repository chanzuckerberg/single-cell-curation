import scanpy as sc
import re
import ontology
import sys
import pandas as pd

field = sys.argv[1]
type_field = sys.argv[2]
dataset_file = sys.argv[3]
lookup_table_file = sys.argv[4]

# Read shendure big data

dataset = sc.read(dataset_file, 'r')
cell_types = dataset.obs.loc[:,[field]].drop_duplicates()
cell_types = list(cell_types.iloc[:,0])
cell_types = pd.DataFrame({field: cell_types, field + 'lookup_name' : cell_types, 'ontology_term_id' : '', 'ontology_term_name' : '',  'final' : 'True', 'type' : type_field})

# Creeate look up table

for i in range(cell_types.shape[0]):
    
    stripped_cell_type = cell_types.loc[i,field].lower()
    stripped_cell_type = re.sub('_', ' ', stripped_cell_type)
    
    try: 
        suggested_term=ontology.lookup_candidate_term(stripped_cell_type, ontology='cl', method='select')
    except:
        suggested_term=[['','']]
        
    
    if len(suggested_term) > 0:
        # store top hit
        cell_types.loc[i,'ontology_term_id'] = suggested_term[0][0]
        cell_types.loc[i,'ontology_term_name'] = suggested_term[0][1]
        
    cell_types.loc[i, field + 'lookup_name'] = stripped_cell_type
    
cell_types.to_csv(lookup_table_file, sep='\t', index=False)
