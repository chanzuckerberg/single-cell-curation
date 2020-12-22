import scanpy as sc
import pandas as pd 
import sys, os

dataset_file = sys.argv[1]
out_prefix = sys.argv[2]
cell_type = sys.argv[3]

dataset = sc.read(dataset_file, 'r')

# Do cell types
class_df = dataset.obs[['subclass_label']]
class_df['subclass_label_lookup'] = class_df['subclass_label']
class_df.drop_duplicates(inplace=True)

if cell_type != "non-neuronal":
    class_df['ontology_term_id'] = cell_type
    class_df['ontology_term_name'] = cell_type
else:
    class_df['ontology_term_id'] = ""
    class_df['ontology_term_name'] = ""
    class_df.loc[class_df['subclass_label'] == 'Astro', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000127', 'astrocyte']
    class_df.loc[class_df['subclass_label'] == 'Endo', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000115', 'endothelial cell']
    class_df.loc[class_df['subclass_label'] == 'Micro-PVM', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000129', 'microgial cell']
    class_df.loc[class_df['subclass_label'] == 'OPC', ['ontology_term_id', 'ontology_term_name']] = ['CL:0002453', 'oligodendrocyte precursor cell']
    class_df.loc[class_df['subclass_label'] == 'Oligo', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000128', 'oligodendrocyte']
    class_df.loc[class_df['subclass_label'] == 'Peri', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000669', 'pericyte cell']
    class_df.loc[class_df['subclass_label'] == 'SMC', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000192', 'smooth muscle cell']
    class_df.loc[class_df['subclass_label'] == 'VLMC', ['ontology_term_id', 'ontology_term_name']] = ['CL:0000708', 'leptomeningeal cel']

class_df['final'] = True

# Do organism
organism = dataset.obs[['orig.ident']]
organism.drop_duplicates(inplace=True)
organism['orig.ident_lookup'] = organism['orig.ident']
organism['organism_ontology_term_id'] = ""
organism['organism'] = ""
organism.loc[organism['orig.ident'] == 'human', ['organism_ontology_term_id', 'organism']] = ['NCBITaxon:10090', 'Mus musculus']
organism.loc[organism['orig.ident'] == 'mouse', ['organism_ontology_term_id', 'organism']] = ['NCBITaxon:9606', 'Homo sapiens']
organism.loc[organism['orig.ident'] == 'marmoset', ['organism_ontology_term_id', 'organism']] = ['NCBITaxon:9483', 'Callithrix jacchus']
organism.loc[organism['orig.ident'] == 'macaque', ['organism_ontology_term_id', 'organism']] = ['NCBITaxon:9539', 'Macaca']


class_df.to_csv(out_prefix + '_ontology_lookup_cell_type.tsv' , sep='\t', index=False)
organism.to_csv(out_prefix + '_ontology_lookup_organism.tsv' , sep='\t', index=False)
