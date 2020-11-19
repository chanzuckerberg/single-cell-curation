# Create yaml for shendure data from the tissue, cell_type, and develomental stage look up tables

library("yaml")

main <- function(cmdArgs=commandArgs(T)) {
    
    cell_type_file <- cmdArgs[1]
    tissue_file <- cmdArgs[2]
    dev_file <- cmdArgs[3]
    out_file <- cmdArgs[4]
    
    #cell_type_file <- "ontology_lookup_cell_type_curated.tsv"
    #tissue_file <- "ontology_lookup_tissue_curated.tsv"
    #dev_file <- "ontology_lookup_dev_stage_curated.tsv"
    #out_file <- "schema-shendure.yml"
    
    # Create backbone of the yml with manual values
    yaml_final <- make_yaml_list()
    
    # Read suggested term tables
    cell_type <- read.table(cell_type_file, sep = '\t', stringsAsFactors = F, header = T)
    tissue <- read.table(tissue_file, sep = '\t', stringsAsFactors = F, header = T)
    dev_stage <- read.table(dev_file, sep = '\t', stringsAsFactors = F, header = T)
    
    # Verify there's only on ontology term per item and convert to list
    cell_type_yaml = lookup_table_to_list(cell_type)
    tissue_yaml = lookup_table_to_list(tissue)
    dev_yaml = lookup_table_to_list(dev_stage)
    
    # Appending to yaml
    yaml_final$obs <- c(yaml_final$obs, list(tissue_ontology_term_id = list(Organ = tissue_yaml),
                                             cell_type_ontology_term_id = list(sub_cluster_name = cell_type_yaml),
                                             development_stage_ontology_term_id = list(Development_day = dev_yaml)
                                             )
                       )
    
    write_yaml(yaml_final, out_file)
    
}


#' Creates backbone of yaml
make_yaml_list <- function() {
    
    list(
         obs = list(assay_ontology_term_id = "EFO:0010550",
                    disease_ontology_term_id = "PATO:0000461",
                    ethnicity_ontology_term_id = "unknown",
                    sex = list(
                               Sex = list(
                                          M = 'male', 
                                          F = 'female')
                              ) 
                    ),
         
         uns = list(version = list(corpora_schema_version="1.1.0", corpora_encoding_version="0.1.0"),
                    organism = "Homo sapiens", 
                    organism_ontology_term_id = "NCBITaxon:9606",
                    layer_descriptions = list(X = "log1p"), 
                    publication_doi = "http://dx.doi.org/10.1126/science.aba7721", 
                    title = "Survey of human embryonic development"
                    )
         )
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
