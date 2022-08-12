import anndata as ad
import numpy as np
from . import ontology


def convert(input_file, output_file):
    print(f"converting {input_file} into {output_file}")
    
    dataset = ad.read_h5ad(input_file)

    # Rename ethnicity_ontology_term_id field
    dataset.obs.rename(columns={"ethnicity_ontology_term_id": "self_reported_ethnicity_ontology_term_id"}, inplace=True)

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

    # Set suspension type
    match_assays = {
        # 'EFO:0010010': ['cell', 'nucleus'], 
        'EFO:0008720': 'nucleus',
        'EFO:0008722': 'cell',
        'EFO:0030002': 'cell',
        'EFO:0008853': 'cell',
        'EFO:0030026': 'nucleus',
        # 'EFO:0010550': ['cell', 'nucleus'], 
        'EFO:0008919': 'cell',
        'EFO:0008939': 'nucleus',
        'EFO:0030027': 'nucleus',
    }

    match_assays_or_children = {
        # 'EFO:0030080': ['cell', 'nucleus'], 
        'EFO:0007045': 'nucleus',
        'EFO:0009294': 'cell',
        # 'EFO:0010184': ['cell', 'nucleus'], 
        'EFO:0009918': 'na',
        'EFO:0700000': 'na',
        'EFO:0030005': 'na',
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

    suspension_type_map = {}
    if dataset.obs["assay_ontology_term_id"].dtype != "category":
        dataset.obs.loc[:, ["assay_ontology_term_id"]] = dataset.obs.astype("category")
    for item in dataset.obs["assay_ontology_term_id"].cat.categories:
        suspension_type_map[item] = assign_suspension_type(item)
    dataset.obs["suspension_type"] = dataset.obs.apply(lambda row: suspension_type_map.get(row.assay_ontology_term_id),
                                                       axis=1)
    dataset.obs.loc[:, ["suspension_type"]] = dataset.obs.astype("category")
    dataset.write(output_file)

