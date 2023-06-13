import anndata as ad

from . import utils


def migrate(input_file, output_file, collection_id, dataset_id):
    print(f"Converting {input_file} into {output_file}")

    dataset = ad.read_h5ad(input_file)

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
    # fmt: off
    deprecated_features_ids = [
        "ENSG00000256374",
        "ENSG00000263464",
        "ENSG00000203812",
        "ENSG00000272196",
        "ENSG00000237133",
        "ENSG00000284741",
        "ENSG00000237838",
        "ENSG00000286699",
        "ENSG00000280250",
        "ENSG00000272370",
        "ENSG00000272354",
        "ENSG00000288639",
        "ENSG00000286228",
        "ENSG00000237513",
        "ENSG00000225932",
        "ENSG00000244693",
        "ENSG00000226403",
        "ENSG00000261534",
        "ENSG00000237548",
        "ENSG00000224745",
        "ENSG00000261438",
        "ENSG00000231575",
        "ENSG00000256863",
        "ENSG00000287388",
        "ENSG00000277077",
        "ENSG00000259444",
        "ENSG00000244952",
        "ENSG00000288630",
        "ENSG00000261963",
        "ENSG00000286065",
        "ENSG00000225178",
        "ENSG00000288593",
        "ENSG00000288546",
        "ENSG00000278198",
        "ENSG00000273496",
        "ENSG00000277666",
        "ENSG00000278782",
        "ENSG00000277761",
    ]
    # fmt: on
    dataset = utils.remove_deprecated_features(dataset, deprecated_features_ids)

    dataset.write(output_file, compression="gzip")