import scanpy as sc
import pandas as pd 
import sys

dataset_file = sys.argv[1]
outfile = sys.argv[2]

dataset = sc.read(dataset_file, 'r')


class_df = dataset.obs[['sub_class']]
class_df.drop_duplicates(inplace=True)
class_df['sub_class_lookup_name'] = class_df['sub_class']
class_df['ontology_term_id'] = ''
class_df['ontology_term_name'] = ''
class_df['final'] = True

class_df.loc[class_df['sub_class'] == 'Astrocyte', 'ontology_term_id'] = 'CL:0000127'
class_df.loc[class_df['sub_class'] == 'Chodl', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['sub_class'] == 'Endothelial', 'ontology_term_id'] = 'CL:0000115'
class_df.loc[class_df['sub_class'] == 'L2/3', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L2/3_Syt10', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L4', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L5', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L5_Chrna6', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L5a', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L5a1', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L5a2', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L5b', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L5b1', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L5b2', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L6a', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L6a1', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L6a2', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'L6b', 'ontology_term_id'] = 'CL:0000679'
class_df.loc[class_df['sub_class'] == 'Microglia', 'ontology_term_id'] = 'CL:0000129'
class_df.loc[class_df['sub_class'] == 'Ndnf', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['sub_class'] == 'OPC', 'ontology_term_id'] = 'CL:0002453'
class_df.loc[class_df['sub_class'] == 'Oligodendrocyte', 'ontology_term_id'] = 'CL:0000128'
class_df.loc[class_df['sub_class'] == 'Pvalb', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['sub_class'] == 'Sst', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['sub_class'] == 'Sst_Cbln4', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['sub_class'] == 'Sst_Chodl', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['sub_class'] == 'Sst_Clstn2', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['sub_class'] == 'Sst_Nmbr', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['sub_class'] == 'Sst_Nr2f2', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['sub_class'] == 'Th', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['sub_class'] == 'Unknown', 'ontology_term_id'] = 'unknown'
class_df.loc[class_df['sub_class'] == 'Vip', 'ontology_term_id'] = 'CL:0000617'
class_df.loc[class_df['sub_class'] == 'nan', 'ontology_term_id'] = 'unknown'

class_df.to_csv(outfile, sep='\t', index=False)
