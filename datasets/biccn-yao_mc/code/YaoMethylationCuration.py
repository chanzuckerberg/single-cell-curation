#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 11:45:50 2021

@author: mlombardo
"""

import mygene
import scanpy as sc
import pandas as pd


mg = mygene.MyGeneInfo()

#Ch methylation - reading in an object that was obtained from Hanqing Liu (author)
mch = sc.read_h5ad("/Users/mlombardo/Documents/dev/czi/curation/Yao2020_epi/mC/objects/MOp.gene.CHN.h5ad")

#Rmove Ensembl version number from Gene IDs
mch_features = mch.var.index.str.split(pat = '.').to_list()
mch_features = [x[0] for x in mch_features]

#Gene symbols from Ensembl IDs - query using my gene
mch_gene_symbols = mg.querymany(mch_features ,
                                scopes='ensembl.gene',
                                fields='symbol',
                                species='mouse')

#create a dataframe to reference
mch_gene_symbols_df = pd.DataFrame(mch_gene_symbols)
mch_gene_symbols_df = mch_gene_symbols_df.set_index('query')
mch_gene_symbols_df = mch_gene_symbols_df.fillna(False)

#remove duplicate hits (ensembl id's that map to more than one gene, keep the first instance)
#could possibly update to consider the score that is returned by mygene
is_duplicate = mch_gene_symbols_df.index.duplicated(keep="first")#only keeps first observation, preserves order
not_duplicate = ~is_duplicate
mch_gene_symbols_df = mch_gene_symbols_df[not_duplicate]

####
#Update the var index with the only the ensembl ids we were able to find corresponding gene names for

found_index = [i for i, x in enumerate(mch_gene_symbols_df['notfound'].to_list()) if not x]
full_var = mch.var.index.to_list()
keep_var = [full_var[i] for i in found_index]

mch_mapped = mch[:,keep_var]
mch_mapped_features = mch_mapped.var.index.str.split(pat = '.').to_list()
mch_mapped_features = [x[0] for x in mch_mapped_features]

###

#pull the corresponding gene symbol in the correct order using the ordered list of ensembl ids
mch_gene_symbols = [mch_gene_symbols_df.loc[x, 'symbol'] for x in mch_mapped_features]

#add gene symbols to var df
mch_mapped.var['gene_symbol'] = mch_gene_symbols

#Change var names from ensembl to gene symbols
mch_mapped.var_names = mch_mapped.var['gene_symbol']

#drop the gene symbol column from the var df
mch_mapped.var = mch_mapped.var.drop(['gene_symbol'], axis = 1)

#make the var names unique
mch_mapped.var_names_make_unique()

#Plot tSNE
sc.pl.tsne(mch_mapped, color= 'Tshz2')



############################
############################
############################


#Cg methylation - reading in an object that was obtained from Hanqing Liu (author)
mcg = sc.read_h5ad("/Users/mlombardo/Documents/dev/czi/curation/Yao2020_epi/mC/objects/MOp.gene.CGN.h5ad")

#Rmove Ensembl version number from Gene IDs
mcg_features = mcg.var.index.str.split(pat = '.').to_list()
mcg_features = [x[0] for x in mcg_features]

#Gene symbols from Ensembl IDs - query using my gene
mcg_gene_symbols = mg.querymany(mcg_features ,
                                scopes='ensembl.gene',
                                fields='symbol',
                                species='mouse')

#create a dataframe to reference
mcg_gene_symbols_df = pd.DataFrame(mcg_gene_symbols)
mcg_gene_symbols_df = mcg_gene_symbols_df.set_index('query')
mcg_gene_symbols_df = mcg_gene_symbols_df.fillna(False)

#remove duplicate hits (ensembl id's that map to more than one gene, keep the first instance)
#could possibly update to consider the score
is_duplicate = mcg_gene_symbols_df.index.duplicated(keep="first")#only keeps first observation, preserves order
not_duplicate = ~is_duplicate
mcg_gene_symbols_df = mcg_gene_symbols_df[not_duplicate]

####
#Update the var index with the only the ensembl ids we were able to find corresponding gene names for

found_index = [i for i, x in enumerate(mcg_gene_symbols_df['notfound'].to_list()) if not x]
full_var = mcg.var.index.to_list()
keep_var = [full_var[i] for i in found_index]

mcg_mapped = mcg[:,keep_var]
mcg_mapped_features = mcg_mapped.var.index.str.split(pat = '.').to_list()
mcg_mapped_features = [x[0] for x in mcg_mapped_features]

###

#pull the corresponding gene symbol in the correct order using the ordered list of ensembl ids
mcg_gene_symbols = [mcg_gene_symbols_df.loc[x, 'symbol'] for x in mcg_mapped_features]

#add gene symbols to var df
mcg_mapped.var['gene_symbol'] = mcg_gene_symbols

#Change var names from ensembl to gene symbols
mcg_mapped.var_names = mcg_mapped.var['gene_symbol']

#drop the gene symbol column from the var df
mcg_mapped.var = mcg_mapped.var.drop(['gene_symbol'], axis = 1)

#make the var names unique
mcg_mapped.var_names_make_unique()

#Plot tSNE
sc.pl.tsne(mcg_mapped, color= 'Tshz2')


#######################
#######################
#######################

#Add the appropriate metadata to our var dataframes (for both cg and ch methylation)
#Load schema info
schema_info = pd.read_csv("/Users/mlombardo/Documents/dev/czi/curation/Yao2020_epi/mC/data/cell_ontology_lookup_mC.csv")

########
#chn
mch_obs_update = mch_mapped.obs.merge(right = schema_info, how = 'left',
                                      left_on = 'SubCluster', right_on = 'BICCN_cluster_label')
mch_mapped.obs = mch_obs_update
#Write updated object
mch_mapped.write_h5ad("/Users/mlombardo/Documents/dev/czi/curation/Yao2020_epi/mC/objects/updatedFeatures/MOp_gene_CHN.h5ad")

########
#cgn
mcg_obs_update = mcg_mapped.obs.merge(right = schema_info, how = 'left',
                                      left_on = 'SubCluster', right_on = 'BICCN_cluster_label')
mcg_mapped.obs = mcg_obs_update
#Write updated object
mcg_mapped.write_h5ad("/Users/mlombardo/Documents/dev/czi/curation/Yao2020_epi/mC/objects/updatedFeatures/MOp_gene_CGN.h5ad")
