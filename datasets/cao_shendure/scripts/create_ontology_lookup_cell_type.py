# that need to be manually curated

import ontology
import re
import scanpy as sc 
import sys

# Read shendure big data
dataset = sc.read('./Survey_of_human_embryonic_development-processed.h5ad', 'r')

cell_types = list(dataset.obs['sub_cluster_name'].unique())

# Creeate look up table

print("original_cell_type", "stripped_cell_type", "ontology_term_id", "ontology_term_name", "final", sep="\t")
for cell_type in cell_types:
    
    stripped_cell_type=cell_type.split("-")[1]
    stripped_cell_type=re.sub("cells*", "", stripped_cell_type).lower()
    
    suggested_term=ontology.lookup_candidate_term(stripped_cell_type, ontology='cl', method='select')
    
    if len(suggested_term) == 1:
        final=True
    elif len(suggested_term) > 1:
        final=False
    else:
        final=False
        suggested_term=[("unkwown", "unkwown")]
        
    for i in suggested_term:
        
        if (i[1] == stripped_cell_type):
            final = True
        else:
            final = False
        
        print(cell_type, stripped_cell_type, *i, final, sep="\t")
        
