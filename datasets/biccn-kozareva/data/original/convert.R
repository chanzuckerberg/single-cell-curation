library('Seurat')
library('loomR')

aaa <- readRDS('cb_annotated_object.RDS')
bbb <- UpdateSeuratObject(aaa)
ccc <- as.loom(bbb)
saveRDS(ccc, ‘cb_annotated_object.loom’)

