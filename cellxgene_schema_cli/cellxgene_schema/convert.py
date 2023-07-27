from warnings import warn

import anndata as ad
import numpy as np

from . import ontology
from .remove_labels import AnnDataLabelRemover


def convert(input_file, output_file):
    warn("This is deprecated, and will be removed in 4.0.0", DeprecationWarning, stacklevel=2)
    print(f"Converting {input_file} into {output_file}")

    dataset = ad.read_h5ad(input_file)

    # Rename ethnicity_ontology_term_id field
    dataset.obs.rename(
        columns={"ethnicity_ontology_term_id": "self_reported_ethnicity_ontology_term_id"},
        inplace=True,
    )

    # Set schema version to 3.0.0
    dataset.uns["schema_version"] = "3.0.0"

    # Remove Deprecated uns Fields
    deprecated_uns_fields = [
        "X_normalization",
        "default_field",
        "layer_descriptions",
        "tags",
        "version",
        "contributors",
        "preprint_doi",
        "project_description",
        "project_links",
        "project_name",
        "publication_doi",
    ]
    for field in deprecated_uns_fields:
        if field in dataset.uns:
            del dataset.uns[field]

    # remove labels that are meant to be added by data portal
    anndata_label_remover = AnnDataLabelRemover(dataset)
    anndata_label_remover.remove_labels()
    dataset = anndata_label_remover.adata

    # remove 'ethnicity' label (needs to be done separately since its no longer in the schema definition)
    if "ethnicity" in dataset.obs:
        del dataset.obs["ethnicity"]

    # Update deprecated ontologies with known direct replacements
    disease_ontology_update_map = {
        "MONDO:0008345": "MONDO:0800029",
        "MONDO:0004553": "MONDO:0017853",
    }
    cell_type_ontology_update_map = {
        "CL:0002609": "CL:0010012",
        "CL:0011107": "CL:0000636",
    }
    assay_ontology_update_map = {
        "EFO:0030002 (BD Rhapsody)": "EFO:0700003",
        "EFO:0010183 (BD Rhapsody)": "EFO:0700003",
    }

    def update_categorical_column_vals(dataframe, column_name, update_map):
        if dataframe[column_name].dtype != "category":
            dataframe[column_name] = dataframe[column_name].astype("category")
        for deprecated_term, new_term in update_map.items():
            if deprecated_term in dataframe[column_name].cat.categories:
                # add new one if not already in category, else continue
                if new_term not in dataframe[column_name].cat.categories:
                    dataframe[column_name] = dataframe[column_name].cat.add_categories(new_term)
                # replace in dataset
                dataframe.loc[dataframe[column_name] == deprecated_term, column_name] = new_term
                # remove deprecated_term from category
                dataframe[column_name] = dataframe[column_name].cat.remove_categories(deprecated_term)

    update_categorical_column_vals(dataset.obs, "disease_ontology_term_id", disease_ontology_update_map)
    update_categorical_column_vals(dataset.obs, "cell_type_ontology_term_id", cell_type_ontology_update_map)
    update_categorical_column_vals(dataset.obs, "assay_ontology_term_id", assay_ontology_update_map)

    # Set suspension type

    # mappings of assays (or assays + child term assays) to corresponding suspension_type
    # valid assays with multiple possible suspension_types shown but commented out
    match_assays = {
        # 'EFO:0010010': ['cell', 'nucleus'],
        "EFO:0008720": "nucleus",
        # 'EFO:0008722': ['cell', 'nucleus'], ,
        "EFO:0030002": "cell",
        "EFO:0008853": "cell",
        "EFO:0030026": "nucleus",
        # 'EFO:0010550': ['cell', 'nucleus'],
        "EFO:0008919": "cell",
        "EFO:0008939": "nucleus",
        "EFO:0030027": "nucleus",
    }

    match_assays_or_children = {
        # 'EFO:0030080': ['cell', 'nucleus'],
        "EFO:0007045": "nucleus",
        "EFO:0009294": "cell",
        # 'EFO:0010184': ['cell', 'nucleus'],
        "EFO:0009918": "na",
        "EFO:0700000": "na",
        "EFO:0030005": "na",
    }

    ontology_checker = ontology.OntologyChecker()

    def assign_suspension_type(item):
        if item in match_assays:
            return match_assays[item]
        else:
            for k, v in match_assays_or_children.items():
                try:
                    if k == item or ontology_checker.is_descendent_of("EFO", item, k):
                        return v
                except Exception:
                    return np.nan
        return np.nan

    if "suspension_type" not in dataset.obs:
        print("column 'suspension_type' not found in obs. Adding column and auto-assigning value where possible.")
        suspension_type_map = {}
        if dataset.obs["assay_ontology_term_id"].dtype != "category":
            dataset.obs.loc[:, ["assay_ontology_term_id"]] = dataset.obs.astype("category")
        for item in dataset.obs["assay_ontology_term_id"].cat.categories:
            suspension_type_map[item] = assign_suspension_type(item)
            if np.isnan(suspension_type_map[item]):
                print(
                    f"Dataset contains row(s) with assay_ontology_term_id {item}. These cannot be auto-assigned a "
                    f"suspension type, please assign a suspension_type manually and validate."
                )
            else:
                print(
                    f"Dataset contains row(s) with assay_ontology_term_id {item}. "
                    f"Automatically assigned suspension_type {suspension_type_map[item]}"
                )
        dataset.obs["suspension_type"] = dataset.obs.apply(
            lambda row: suspension_type_map.get(row.assay_ontology_term_id), axis=1
        )
        dataset.obs.loc[:, ["suspension_type"]] = dataset.obs.astype("category")
    else:
        if dataset.obs["suspension_type"].dtype != "category":
            dataset.obs.loc[:, ["suspension_type"]] = dataset.obs.astype("category")
        print(f"suspension_type already exists in obs, with categories {dataset.obs['suspension_type'].cat.categories}")

    print(
        f"Automatable conversions completed. Please run 'cellxgene-schema validate {output_file}' to check for "
        f"required manual changes, if any."
    )
    dataset.write(output_file, compression="gzip")
