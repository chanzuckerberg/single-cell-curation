library('sceasy')
library('reticulate')
library('loomR')
library('Seurat')

main <- function(cmdArgs=commandArgs(T)) {
    
    from <- cmdArgs[1]
    to <- cmdArgs[2]
    conda_environment <- cmdArgs[3] # Empty srping if none
    input_file <- cmdArgs[4]
    output_file <- cmdArgs[5]
    
    if(nchar(conda_environment) > 0)
        use_condaenv(conda_environment)
    
    dataset <- readRDS(input_file)
    
    convertFormat(dataset, from=from, to=to, outFile=output_file)
    
    
}

main()
