# that need to be manually curated

import ontology
import re
import scanpy as sc 

# Read shendure big data
dataset = sc.read('./Survey_of_human_embryonic_development-processed.h5ad', 'r')

tissues = list(dataset.obs['Organ'].unique())

# Creeate look up table

print("original_tisue", "stripped_tissue", "ontology_term_id", "ontology_term_name", "final", sep="\t")
for tissue in tissues:
    
    stripped_tissue = tissue.lower()
    
    suggested_term = ontology.lookup_candidate_term(stripped_tissue, ontology='uberon', method='select')
    
    if len(suggested_term) == 0:
        suggested_term = [("unkwown", "unkwown")]
        
    for i in suggested_term:
        
        if (i[1] == stripped_tissue): 
            final = True
        else:
            final = False
            
        print(tissue, stripped_tissue, *i, final, sep="\t")
