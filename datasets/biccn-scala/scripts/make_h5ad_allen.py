import pickle 
import anndata

# Read files
h5ad_10x = anndata.read('./data/original/10X_cells_v2_AIBS.h5ad')
h5ad_10x.obs_names = [re.sub(r"(\w+\-\d+).+", r"\1", i) for i in h5ad_10x.obs_names]
h5ad_10x.X = h5ad_10x.raw.X
h5ad_10x.uns['layer_descriptions']: {"X": "raw"}
del h5ad_10x.raw

for i in list(h5ad_10x.obsm):
    del h5ad_10x.obsm[i]

with open('./data/original/10X_cells_v2_AIBS.pickle', 'rb') as file_in:
    data = pickle.load(file_in)
del data['neurons']

tsne = {}
with open('./data/original/10x-tsne-exc.pickle', 'rb') as file_in:
    tsne['exc'] = pickle.load(file_in)
    
with open('./data/original/10x-tsne-pvsst.pickle', 'rb') as file_in:
    tsne['pvsst'] = pickle.load(file_in)
    
with open('./data/original/10x-tsne-viplamp.pickle', 'rb') as file_in:
    tsne['viplamp'] = pickle.load(file_in)
    
    
# Make anndatas
for i in data:
    current = h5ad_10x[data[i]['cells'], ]
    current.obsm['X_tsne'] = tsne[i]
    current.write_h5ad('./data/transformed/1_refformated/10X_cells_v2_AIBS_' + i + '.h5ad', compression='gzip')
    
