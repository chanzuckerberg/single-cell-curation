import anndata as ad
import numpy as np

from . import ontology
from . import utils


def convert(input_file, output_file, collection_id, dataset_id):
    print(f"Converting {input_file} into {output_file}")

    dataset = ad.read_h5ad(input_file)

    # Set schema version
    dataset.uns["schema_version"] = "3.1.0"

    # ONTOLOGY TERMS TO UPDATE ACROSS ALL DATASETS IN CORPUS
    # Initialization is AUTOMATED for newly deprecated terms that have 'Replaced By' terms in their ontology files

    # Curators should review the monthly 'Curator Report' and add deprecated term replacements to corresponding map if
    # 'Replaced By' is not available for a deprecated term.

    # If Curators have non-deprecated term changes to apply to all datasets in the corpus where applicable,
    # add them here.
    ontology_term_maps = {
        "assay": {
            "EFO:0030002 (BD Rhapsody)": "EFO:0700003",  # AUTOMATED
            "EFO:0010183 (BD Rhapsody)": "EFO:0700003",  # AUTOMATED
        },
        "cell_type": {
            "CL:0002609": "CL:0010012",  # AUTOMATED
            "CL:0011107": "CL:0000636",  # Curator-added (demonstrative example, curators do not need to annotate)
        },
        "development_stage": {},
        "disease": {
            "MONDO:0008345": "MONDO:0800029",  # Curator-added (demonstrative example, curators do not need to annotate)
        },
        "organism": {},
        "self_reported_ethnicity": {},
        "sex": {},
        "tissue": {},
    }

    # AUTOMATED, DO NOT CHANGE
    for ontology_name, deprecated_term_map in ontology_term_maps.items():
        utils.replace_ontology_term(dataset.obs, ontology_name, deprecated_term_map)

    # CURATOR-DEFINED, DATASET-SPECIFIC UPDATES
    # Use the template below to define dataset and collection specific ontology changes. Will only apply to dataset
    # if it matches a condition.
    # If no such changes are needed, leave blank

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

    # Example Curator Input
    df = dataset.obs
    if collection_id == "a48f5033-3438-4550-8574-cdff3263fdfd":
        utils.replace_ontology_term(df, "assay", {"EFO:0008913": "EFO:0700010"})
    elif dataset_id == "5cdbb2ea-c622-466d-9ead-7884ad8cb99f":
        utils.replace_ontology_term(df, "cell_type", {"CL:0000561": "CL:4030027"})
    elif dataset_id == "f8c77961-67a7-4161-b8c2-61c3f917b54f":
        df["cell_type_ontology_term_id"] = df["cell_type_ontology_term_id"].cat.add_categories("CL:0004219")
        df.loc[
            (df['author_cell_type'] == "1.0") & (df['cell_type_ontology_term_id'] == "CL:0000561"),
            "cell_type_ontology_term_id"
        ] = "CL:0004219"

        df.loc[
            (df['author_cell_type'] == "5.0") & (df['cell_type_ontology_term_id'] == "CL:0000561"),
            "cell_type_ontology_term_id"
        ] = "CL:4030028"

        df.loc[
            (df['author_cell_type'] == "11.0") & (df['cell_type_ontology_term_id'] == "CL:0000561"),
            "cell_type_ontology_term_id"
        ] = "CL:0004232"

    # AUTOMATED, DO NOT CHANGE -- IF GENCODE UPDATED, DEPRECATED FEATURE FILTERING ALGORITHM WILL GO HERE.
    # No Changes

    dataset.write(output_file, compression="gzip")
