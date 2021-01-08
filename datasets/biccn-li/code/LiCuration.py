#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 13:14:40 2020

@author: mlombardo
"""

import scanpy as sc
import pandas as pd
import numpy as np

#creating anndata obejcts from the Li data


###############################################################################\

##########################################
#Load and process supplementary tables
#Region
supp1 = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/SupplementalTable1.tsv", sep = "\t")

#Cell Type
supp3 = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/SupplementalTable3.tsv", sep = "\t")
supp3['BICCN Ontology'] = supp3['BICCN Ontology'].str.split(' ').str[0]#Clean flagged annotations
supp3 = supp3[["Major Type", "BICCN Ontology"]]
supp3.columns = ["Major_Type", "BICCN_ontology_term_id"]
supp3 = supp3.drop_duplicates(subset="Major_Type")

#Read in extra table to add non-neuronal CL ontologies - no duplicates, don't need to drop
cl_ontologies = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/NonNCellTypes_CLOntologies_byMax.tsv", sep = "\t")

###############################################################################
#Non-Neuronal Cells
#Read in data and create the anndata object
nonN = sc.read_text("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/NonN/exprMatrix.tsv.gz")
#Transpose observations and variables since the original data matrix was tranposed somehow
nonN = nonN.transpose()

#Load observation metadata 
meta = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/NonN/meta.tsv", sep="\t")
#Annotate the observation metadata with info from supplementary tables
meta = meta.merge(supp3, how = 'left', left_on = 'L2cluster', right_on = "Major_Type")
meta = meta.merge(cl_ontologies, how = 'left', left_on = 'L2cluster', right_on = 'cell_type')

#Rename L2 cluster to BICCN cluster label
meta = meta.rename(columns = {'L2cluster':'BICCN_subclass_label'})

#Create BICCN type column
meta["BICCN_class_label"] = np.repeat("ILX:0770099", meta.shape[0])

#Drop extra columns
meta = meta.drop(labels = ["Major_Type"], axis = 1)

nonN.obs = meta



#Load umap coordinates and add to object (ensuring that they are numpy array, not pandads df)
umap = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/NonN/umap.coords.tsv.gz",
                   sep="\t",
                   index_col = 0, 
                   names = ["UMAP_1", "UMAP_2"])
nonN.obsm["X_umap"] = umap.to_numpy()

nonN.write("/Users/mlombardo/Documents/dev/czi/curation/Li2020/objects/non_neuronal_annotated.h5ad")

###############################################################################
#Gaba Neurons
#Read in data and create the anndata object
gaba = sc.read_text("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/GabaN/exprMatrix.tsv.gz")
#Transpose observations and variables since the original data matrix was tranposed somehow
gaba = gaba.transpose()

########################
#Load observation metadata 
meta = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/GabaN/meta.tsv", sep="\t")
#Annotate the observation metadata with info from supplementary tables
meta = meta.merge(supp3, how = 'left', left_on = 'L2cluster', right_on = "Major_Type")
#Rename L2 cluster to BICCN cluster label
meta = meta.rename(columns = {'L2cluster':'BICCN_subclass_label'})
meta["BICCN_class_label"] = np.repeat("ILX:0770098", meta.shape[0])
meta = meta.drop(labels = ["Major_Type"], axis = 1)


gaba.obs = meta
########################
#Load umap coordinates and add to object (ensuring that they are numpy array, not pandads df)
umap = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/GabaN/umap.coords.tsv.gz",
                   sep="\t",
                   index_col = 0, 
                   names = ["UMAP_1", "UMAP_2"])
gaba.obsm["X_umap"] = umap.to_numpy()

gaba.write("/Users/mlombardo/Documents/dev/czi/curation/Li2020/objects/gaba_neurons_annotated_latest.h5ad")


###############################################################################
#Gluta Neurons
#Read in data and create the anndata object
gluta = sc.read_text("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/GlutaN/exprMatrix.tsv.gz")
#Transpose observations and variables since the original data matrix was tranposed somehow
gluta = gluta.transpose()#May not need for this one?
#Load observation metadata 
meta = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/GlutaN/meta.tsv", sep="\t")

#Annotate the observation metadata with info from supplementary tables
meta = meta.merge(supp3, how = 'left', left_on = 'L2cluster', right_on = "Major_Type")

#Rename L2 cluster to BICCN cluster label
meta = meta.rename(columns = {'L2cluster':'BICCN_subclass_label'})
meta["BICCN_class_label"] = np.repeat("ILX:0770097", meta.shape[0])
meta = meta.drop(labels = ["Major_Type"], axis = 1)


gluta.obs = meta
#Load umap coordinates and add to object (ensuring that they are numpy array, not pandads df)
umap = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/GlutaN/umap.coords.tsv.gz",
                   sep="\t",
                   index_col = 0, 
                   names = ["UMAP_1", "UMAP_2"])
gluta.obsm["X_umap"] = umap.to_numpy()

gluta.write("/Users/mlombardo/Documents/dev/czi/curation/Li2020/objects/gluta_neurons_annotated_latest.h5ad")




##########################################
#Annotate the observation metadata with info from supplementary tables
#Region
#supp1 = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/SupplementalTable1.tsv", sep = "\t")
#Cell Type
#supp3 = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/SupplementalTable3.tsv", sep = "\t")
#supp3['BICCN Ontology'] = supp3['BICCN Ontology'].str.split(' ').str[0]
#supp3 = supp3[["Major Type", "BICCN Ontology"]]
#supp3.columns = ["Major_Type", "BICCN_CellType_Ontology"]
#supp3 = supp3.drop_duplicates(subset="Major_Type")

#Merge with metadata
#meta = meta.merge(supp3, how = 'left', left_on = 'L2cluster', right_on = "Major_Type")

