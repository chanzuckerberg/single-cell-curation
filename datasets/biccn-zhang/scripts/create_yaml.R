# Create cellxgene schema yaml file

library("yaml")

main <- function(cmdArgs=commandArgs(T)) {
    
    cell_type_file <- cmdArgs[1]
    out_file <- cmdArgs[2]
    
    #cell_type_file <- '../data/misc/ontology_lookup_cell_type.tsv'
    
    # Create backbone of the yml with manual values
    yaml_final <- make_yaml_list()
    
    # Read suggested term tables
    cell_type <- read.table(cell_type_file, sep = '\t', stringsAsFactors = F, header = T)
    
    # Verify there's only on ontology term per item and convert to list
    cell_type_yaml = lookup_table_to_list(cell_type)
    
    # Appending to yaml
    yaml_final$obs <- c(yaml_final$obs, list(cell_type_ontology_term_id = list(BICCN_subclass_label = cell_type_yaml)))
    yaml_final <- as.yaml(yaml_final)
    
    # Eliminate single quotes in numbers
    #yaml_final <- gsub("\\'(\\d+)\\'", "\\1", yaml_final)
    
    writeLines(yaml_final, out_file)
    
}


#' Creates backbone of yaml
make_yaml_list <- function() {
    
    x <- list(
         #fixup_gene_symbols = list(X = "log1p", raw.X = "raw"),
         obs = list(assay_ontology_term_id = "EFO:0008992",
                    ethnicity_ontology_term_id = "na",
                    sex = "male",
                    disease_ontology_term_id =  'PATO:0000461',
                    tissue_ontology_term_id = 'UBERON:0001384',
                    development_stage_ontology_term_id = "EFO:0001272"
                    ),
         
         uns = list(version = list(corpora_schema_version="1.1.0", corpora_encoding_version="0.1.0"),
                    organism = "Mus musculus", 
                    organism_ontology_term_id = "NCBITaxon:10090",
                    layer_descriptions = list(X = "raw merfish (counts/cell volume)", raw.X = "raw"), 
                    publication_doi = "https://doi.org/10.1101/2020.06.04.105700", 
                    title = paste("Molecular, spatial and projection diversity of neurons in primary motor cortex revealed by in situ single-cell transcriptomics")
                    )
         )
    
    x
}
    


#' Takes the look up table and keeps the selected terms, throws an error if there's more than selected ontology terms per item
lookup_table_to_list <- function(x) {
    
    #if(! all(x[, 5, drop = T] %in% c('True', 'False')))
    #    stop('Column five from look up tables can only be "True" or "False"')
    
    x$final <- ifelse(x$final == 'True', T, F)
    
    validate_lookup_table(x)
    
    x <- x[x$final, ]
    
    results  <- setNames(x[, 3, drop = T], x[, 1, drop = T])
    results <- as.list(results)
    
    return(results)
    
}

validate_lookup_table <- function(x) {
    
    tallies <- tapply(x$final, x[,1,drop = T], sum)
    
    if(any(tallies!=1)){
        incorrect <- names(tallies[tallies!=1])
        stop("The following categories don't have exactly one ontology selected", paste(incorrect, collapse = ", "))
    }
}

main()
