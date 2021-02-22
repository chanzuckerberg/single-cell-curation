# that need to be manually curated

import ontology
import re
import sys
import pandas as pd

# Read shendure big data

cell_types = pd.read_csv('./data/original/tabula_muris_senis_cell_types.csv')
cell_types[['id', 'term', 'other_ids', 'other_terms']] = ''

cache=dict()

# Creeate look up table

for i in range(cell_types.shape[0]):
    
    print(i)
    
    stripped_cell_type = cell_types.loc[i,'cell_ontology_class'].lower()
    
    
    if stripped_cell_type in cache:
        cell_types.loc[i, ['id', 'term', 'other_ids', 'other_terms']] = cache[stripped_cell_type]
        continue
    
    print('looking up term')
    suggested_term=ontology.lookup_candidate_term(stripped_cell_type, ontology='cl', method='select')
    
    
    # store top hit
    cell_types.loc[i,'id'] = suggested_term[0][0]
    cell_types.loc[i,'term'] = suggested_term[0][1]
    
    # store other hits
    if len(suggested_term) > 1:
        other_ids=[]
        other_terms=[]
        for j in range(1, len(suggested_term)):
            if "CL" in suggested_term[j][0]:
                other_ids.append(suggested_term[j][0])
                other_terms.append(suggested_term[j][1])
            
        cell_types.loc[i, 'other_ids'] = ";".join(other_ids)
        cell_types.loc[i, 'other_terms'] = ";".join(other_terms)
    
    cache[stripped_cell_type] = cell_types.loc[i,['id', 'term', 'other_ids', 'other_terms']]
    
cell_types.to_csv('./data/original/tabula_muris_senis_cell_types_with_ontology.csv')
