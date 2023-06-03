import anndata as ad

from . import utils


def convert(input_file, output_file, collection_id, dataset_id):
    print(f"Converting {input_file} into {output_file}")

    dataset = ad.read_h5ad(input_file)

    # Set schema version
    dataset.uns["schema_version"] = "make[1]: Leaving directory '/home/runner/work/single-cell-curation/single-cell-curation'"

    # ONTOLOGY TERMS TO UPDATE ACROSS ALL DATASETS IN CORPUS
    # Initialization is AUTOMATED for newly deprecated terms that have 'Replaced By' terms in their ontology files

    # Curators should review the monthly 'Curator Report' and add deprecated term replacements to corresponding map if
    # 'Replaced By' is not available for a deprecated term.

    # If Curators have non-deprecated term changes to apply to all datasets in the corpus where applicable,
    # add them here.
    ontology_term_maps = {
        "assay": {
        },
        "cell_type": {
        },
        "development_stage": {
        },
        "disease": {
        },
        "organism": {
        },
        "self_reported_ethnicity": {
        },
        "sex": {
        },
        "tissue": {
        },
    }

    # AUTOMATED, DO NOT CHANGE
    for ontology_name, deprecated_term_map in ontology_term_maps.items():
        utils.replace_ontology_term(dataset.obs, ontology_name, deprecated_term_map)

    # CURATOR-DEFINED, DATASET-SPECIFIC UPDATES
    # Use the template below to define dataset and collection specific ontology changes. Will only apply to dataset
    # if it matches a condition.
    # If no such changes are needed, leave blank
    # Examples:
    # if dataset_id == "<dataset_1_id>":
    #   <no further logic necessary>
    #   utils.replace_ontology_term(df, <ontology_name>, {"term_to_replace": "replacement_term", ...})
    # elif dataset_id == "<dataset_2_id>":
    #   <custom transformation logic beyond scope of util functions>
    # elif collection_id == "<collection_1_id>":
    #   <no further logic necessary>
    #   utils.replace_ontology_term(df, <ontology_name>, {"term_to_replace": "replacement_term", ...})
    # elif collection_id == "<collection_2_id>":
    #   <custom transformation logic beyond scope of replace_ontology_term>
    # ...

    # AUTOMATED, DO NOT CHANGE -- IF GENCODE UPDATED, DEPRECATED FEATURE FILTERING ALGORITHM WILL GO HERE.
    # No Changes

    dataset.write(output_file, compression="gzip")