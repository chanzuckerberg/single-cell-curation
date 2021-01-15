#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 11:29:29 2021

@author: mlombardo
"""

import pandas as pd 

#Part One: Ensuring that Yao (https://doi.org/10.1101/2020.02.29.970558) ATAC data is subset of the Li (https://doi.org/10.1101/2020.05.10.087585) dataset

#Load study metadata from the Li and Yao papers
yao_meta = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Yao2020_epi/data/CEMBA_MOp.L2.cluster.meta.txt", sep = "\t")

LiGlutMeta = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/GlutaN/meta.tsv", sep = '\t')

LiGabaMeta = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/GabaN/meta.tsv", sep = '\t')

LiNonNMeta = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Li2020/data/NonN/meta.tsv", sep = '\t')


#Check if nuclei barcodes in Yao are present in any of the dataset in Li

yao_barcodes = yao_meta["barcode"]

#LiGlutCodes = LiGlutMeta['cellID']
#LiGlutCodes = [x.split('.')[1] for x in LiGlutCodes]

#LiGabaCodes = LiGabaMeta['cellID']
#LiGabaCodes = [x.split('.')[1] for x in LiGabaCodes]

#LiNonNCodes = LiNonNMeta['cellID']
#LiNonNCodes = [x.split('.')[1] for x in LiNonNCodes]


#yao_barcodes.isin(LiGlutCodes).value_counts()#True     59722, False    21474
#yao_barcodes.isin(LiGabaCodes).value_counts()#True     35906, False    45290
#yao_barcodes.isin(LiNonNCodes).value_counts()#True     43943, False    37253

#yao_barcodes = yao_barcodes.to_list()

###############3
#Identifying appropriate nuclei from the Li objects to carry into the Yao objects
####################

yao_merged_barcodes = yao_meta['sample'] + '.' + yao_meta['barcode']


LiGlutCodes = LiGlutMeta['cellID']
glut_subset = list(set(LiGlutCodes) & set(yao_merged_barcodes)) 

LiGabaCodes = LiGabaMeta['cellID']
gaba_subset = list(set(LiGabaCodes) & set(yao_merged_barcodes)) 

LiNonNCodes = LiNonNMeta['cellID']
nonN_subset = list(set(LiNonNCodes) & set(yao_merged_barcodes)) 

##########
#Get the Anndata objects that were created for Li2020

import scanpy as sc
import anndata as ad

#Read in remixed objects that were made for Li2020
glut_obj = sc.read_h5ad("/Users/mlombardo/Documents/dev/czi/curation/Li2020/objects/remixed/gluta_neurons_remixed.h5ad")
gaba_obj = sc.read_h5ad("/Users/mlombardo/Documents/dev/czi/curation/Li2020/objects/remixed/gaba_neurons_remixed.h5ad")
nonN_obj = sc.read_h5ad("/Users/mlombardo/Documents/dev/czi/curation/Li2020/objects/remixed/non_neuronal_remixed.h5ad")


glut_obj = glut_obj[glut_obj.obs['cellID'].isin(glut_subset)]
gaba_obj = gaba_obj[gaba_obj.obs['cellID'].isin(gaba_subset)]
nonN_obj = nonN_obj[nonN_obj.obs['cellID'].isin(nonN_subset)]

#Merge
merged_ad = ad.concat(adatas = [glut_obj, gaba_obj, nonN_obj], join = 'outer')

#Update the Yao metadata table with the merged sample and cell barcodes column - titled: 'cellID'
yao_meta_update = yao_meta.copy()
yao_meta_update['cellID_Yao'] = yao_meta_update['sample'] + '.' + yao_meta_update['barcode']

#Filter the the Yao embedddings for the cells that are in our merged anndata object
yao_embeddings = yao_meta_update[['cellID_Yao' ,'tsne1', 'tsne2', 'umap-1', 'umap-2']]
yao_embeddings_filtered = yao_embeddings.iloc[yao_embeddings['cellID_Yao'].isin(merged_ad.obs['cellID']).to_list()]
yao_embeddings_filtered.set_index('cellID_Yao', inplace = True)
yao_embeddings_filtered = yao_embeddings_filtered[merged_ad.obs['cellID'],]
#Ensure that the ordering of the indices/cellIDs matches that of the merged object
yao_embeddings_filtered = yao_embeddings_filtered.reindex(merged_ad.obs['cellID'])

#create the separate tsne and umap datafrmaes/numpy matrices
umap = yao_embeddings_filtered[['umap-1', 'umap-2']]
umap = umap.rename(columns = {"umap-1": "UMAP_1", "umap-2" : "UMAP_2"})
umap = umap.to_numpy()

tsne = yao_embeddings_filtered[['tsne1','tsne2']]
tsne = tsne.rename(columns = {"tsne1" : "TSNE_1", "tsne2" : "TSNE_2"})
tsne = tsne.to_numpy()

merged_ad.obsm["X_umap"] = umap#assign umap embeddings to merged object
merged_ad.obsm["X_tsne"] = tsne#assign tsne embeddings to merged object

#Test embeddings
sc.pl.umap(merged_ad)
sc.pl.tsne(merged_ad)

#Test known marker
sc.pl.tsne(merged_ad, color = 'TSHZ2')


#Save the object
merged_ad.write_h5ad("/Users/mlombardo/Documents/dev/czi/curation/Yao2020_epi/objects/Yao2020ATAC_allCells.h5ad")




