import scanpy as sc
import sys

in_file =  sys.argv[1]
column = sys.argv[2]
out_file = sys.argv[3]

# Remove nans
dataset = sc.read(in_file)
dataset = dataset[dataset.obs[column] != 'nan']

# Rename columns BICCN colums
dataset.obs.rename(inplace=True, columns={'cluster_id':'BICCN_cluster_id', 'cluster_label':'BICCN_cluster_label', 'subclass_label':'BICCN_subclass_label', 'class_label':'BICCN_class_label'})
for i in dataset.obs.columns:
    if i.lower() == 'sex':
        dataset.obs.rename(columns={i:'sex'}, inplace=True)
    if i.lower() == 'gender':
        dataset.obs.rename(columns={i:'sex'}, inplace=True)
        
# Create a temp column which will be used for ontology mapping
dataset.obs['temp_class_label'] = dataset.obs['BICCN_class_label'].astype(str) + dataset.obs['BICCN_subclass_label'].astype(str)


dataset.write(out_file)
