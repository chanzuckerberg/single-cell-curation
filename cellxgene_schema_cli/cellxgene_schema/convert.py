import anndata as ad
import pandas as pd

def convert(input_file, output_file):
    print(f"converting {input_file} into {output_file}")
    
    dataset = ad.read_h5ad(input_file)

    # update obs entries to remove clarifying text from assay_ontology_term_id (for testing fixture, not required by 2to3)
    dataset.obs.loc[dataset.obs["assay_ontology_term_id"] == 'EFO:0009922 (sci-plex)', 'assay_ontology_term_id'] = 'EFO:0009922'
    dataset.obs["assay_ontology_term_id"].cat.remove_unused_categories(inplace=True)

    # Rename ethnicity_ontology_term_id field
    dataset.obs.rename(columns={"ethnicity_ontology_term_id":"self_reported_ethnicity_ontology_term_id"}, inplace=True)
    # add 'multiethnic' as a self_reported_ethnicity_ontology_term_id category (for testing fixture, not required by 2to3)
    dataset.obs["self_reported_ethnicity_ontology_term_id"].cat.add_categories("multiethnic", inplace=True)
    dataset.obs["self_reported_ethnicity_ontology_term_id"][3] = "multiethnic"

    # Set schema version to 3.0.0
    dataset.uns["schema_version"] = "3.0.0"

    dataset.write(output_file)

