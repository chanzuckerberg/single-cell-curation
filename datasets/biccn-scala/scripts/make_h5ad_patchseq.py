# adapted from https://github.com/berenslab/mini-atlas/blob/master/code/patch-seq-data-load.ipynb

import numpy as np
import scanpy as sc
import math
#import pylab as plt
#import seaborn as sns
import pandas as pd
import pickle
import scipy.sparse
import scipy
import time
import warnings
import os

OUTPUT_FOLDER = './data/transformed/1_refformated'

# Read metadata
meta = pd.read_csv('./data/original/m1_patchseq_meta_data.csv', sep='\t')
cells = meta['Cell'].values
meta.set_index("Cell", inplace=True)

# TRANSCRIPTOMIC DATA: read counts

data_exons = pd.read_csv('./data/original/m1_patchseq_exon_counts.csv.gz', na_filter=False, index_col=0)
exonCounts = data_exons.values.transpose()
exonCounts = scipy.sparse.csr_matrix(exonCounts)
assert(all(cells==data_exons.columns))
genes = np.array(data_exons.index)

data_introns = pd.read_csv('./data/original/m1_patchseq_intron_counts.csv.gz', na_filter=False, index_col=0)
intronCounts = data_introns.values.transpose()
intronCounts = scipy.sparse.csr_matrix(intronCounts)
assert(all(cells==data_introns.columns))
assert(all(data_introns.index==data_exons.index))

# EPHYS DATA

ephysData = pd.read_csv('./data/original/m1_patchseq_ephys_features.csv', index_col=0)

# MORPH DATA

morphometrics = pd.read_csv('./data/original/m1_patchseq_morph_features.csv', index_col=0)

        
# Read tsnes
tsne_general = pd.read_csv("./data/original/coverage-tsnes.csv", index_col=0).dropna()
tsne_features = pd.read_csv("./data/original/morpho-electric-tsnes.csv", index_col=0)

# Correct "/" in column names
for i in meta.columns:
    if "/" in i:
        meta.rename(columns={i:i.replace('/', '')}, inplace=True)
        

general = sc.AnnData(exonCounts + intronCounts, obs=meta, var=pd.DataFrame(genes, columns=['gene_symbol']).set_index("gene_symbol"))
#general['exons'] = sc.AnnData(exonCounts, obs=meta, var=pd.DataFrame(genes, columns=['gene_symbol']).set_index("gene_symbol"))
#general['inrons'] = sc.AnnData(intronCounts, obs=meta, var=pd.DataFrame(genes, columns = ['gene_symbol']).set_index("gene_symbol"))

# Save merged data with 10x AIBS
for neuron_type in tsne_general["Which t-SNE"].drop_duplicates():
    tsne = tsne_general.loc[tsne_general['Which t-SNE'] == neuron_type,["t-SNE position x", "t-SNE position y"]]
    
    x = general
    x = x[tsne.index,]
    x.obsm["X_tsne"] = tsne.to_numpy()
    neuron_type = neuron_type.replace("/", "")
    x.write_h5ad(os.path.join(OUTPUT_FOLDER, "patchseq_nonMorpho_" + neuron_type + ".h5ad"), compression='gzip')

# Save electro physiological h5ad
current_tsne = tsne_features[["Ephys t-SNE x", "Ephys t-SNE y"]].dropna()

x = general
x = x[current_tsne.index,]

x.obs = x.obs.join(ephysData)

x.obsm["X_tsne"] = current_tsne.to_numpy()
x.write_h5ad(os.path.join(OUTPUT_FOLDER, "patchseq_electro_phys.h5ad"), compression='gzip')

    
# Save morphological h5ad
current_tsne = tsne_features[["Morph t-SNE x", "Moprh t-SNE y"]].dropna()

x = general
x = x[current_tsne.index,]

x.obs = x.obs.join(morphometrics)

x.obsm["X_tsne"] = current_tsne.to_numpy()
x.write_h5ad(os.path.join(OUTPUT_FOLDER, "patchseq_morphological.h5ad"), compression='gzip')
    
# Save morphoelectro h5ad
current_tsne = tsne_features[["Morphoelectric t-SNE x", "Morphoelectric t-SNE y"]].dropna()

x = general
x = x[current_tsne.index,]

x.obs = x.obs.join(ephysData)
x.obs = x.obs.join(morphometrics)

x.obsm["X_tsne"] = current_tsne.to_numpy()
x.write_h5ad(os.path.join(OUTPUT_FOLDER, "patchseq_morpholphys.h5ad"), compression='gzip')
