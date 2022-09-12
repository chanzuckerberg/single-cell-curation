import anndata as ad


def remove_labels_inplace(dataset):
    labeled_var_fields = [
        "feature_name",
        "feature_reference",
        "feature_biotype",
    ]

    labeled_obs_fields = [
        "cell_type",
        "assay",
        "disease",
        "ethnicity",
        "organism",
        "sex",
        "tissue",
        "self_reported_ethnicity",
        "development_stage",
    ]

    for label in labeled_var_fields:
        if label in dataset.var:
            del dataset.var[label]
        try:  # raw may not be defined
            if label in dataset.raw.var:
                del dataset.raw.var[label]
        except Exception:
            pass

    for label in labeled_obs_fields:
        if label in dataset.obs:
            del dataset.obs[label]


def remove_labels(input_file, output_file):
    print(f"Creating copy of {input_file} with portal-added labels removed and saving as {output_file}")

    dataset = ad.read_h5ad(input_file)

    remove_labels_inplace(dataset)

    dataset.write(output_file)