# Reformats from seurat to anndata
# If sct transform is present saves that to main layer and saves counts to X

library('sceasy')
library('reticulate')
library('Seurat')

main <- function(cmdArgs=commandArgs(T)) {
    
    conda_environment <- cmdArgs[1] # Empty srping if none
    input_file <- cmdArgs[2]
    output_file <- cmdArgs[3]
    
    if(nchar(conda_environment) > 0)
        use_condaenv(conda_environment)
    
    dataset <- readRDS(input_file)
    
    
    if('SCT' %in% names(dataset)) {
        ann <- convertFormat(dataset, from='seurat', to='anndata', outFile=NULL, assay='SCT', main_layer='data')
        mat <- Seurat::GetAssayData(object=dataset, assay='RNA', slot='counts')
        mat <- Matrix::t(mat)
        mat <- mat[py_to_r(ann$obs_names$tolist()), py_to_r(ann$var_names$tolist()) ]
        
        ann_raw <- ann$copy()
        ann_raw$X <- mat
        
        ann$raw <- ann_raw
    } else {
        ann <- convertFormat(dataset, from='seurat', to='anndata', outFile=NULL)
    }
    
    ann$write(output_file)
}

main()
