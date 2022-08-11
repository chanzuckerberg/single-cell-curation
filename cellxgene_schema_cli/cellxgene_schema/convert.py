import anndata as ad
import pandas as pd

def convert(input_file, output_file):
    print(f"converting {input_file} into {output_file}")
    
    dataset = ad.read_h5ad(input_file)

    # Rename ethnicity_ontology_term_id field
    dataset.obs.rename(columns={"ethnicity_ontology_term_id":"self_reported_ethnicity_ontology_term_id"}, inplace=True)

    # Set schema version to 3.0.0
    dataset.uns["schema_version"] = "3.0.0"

    # Remove X_normalization
    del dataset.uns["X_normalization"]

    # Set suspension type
    # TODO

    dataset.write(output_file)

