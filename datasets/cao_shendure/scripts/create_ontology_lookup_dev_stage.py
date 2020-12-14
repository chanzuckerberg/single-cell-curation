# that need to be manually curated

import ontology
import scanpy as sc 

# Read shendure big data
dataset = sc.read('./Survey_of_human_embryonic_development-processed.h5ad', 'r')

days = list(dataset.obs['Development_day'].unique())
days.sort()

print("original_day", "week", "ontology_term_id", "ontology_term_name", "final", sep="\t")
for day in days: 
    
    week = int(day/7) + 1
    week_string = str(week) + "th week post-fertilization human stage"
    
    suggested_term = ontology.lookup_candidate_term(week_string, ontology='HsapDv', method='select')
    
    if len(suggested_term) == 0:
        suggested_term = [("unkwown", "unkwown")]
        
    for i in suggested_term:
        
        if (i[1] == week_string): 
            final = True
        else:
            final = False
            
        print(day, week, *i, final, sep="\t")
