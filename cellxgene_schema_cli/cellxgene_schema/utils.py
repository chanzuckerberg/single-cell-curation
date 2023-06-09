import os
from typing import List

import anndata as ad
import env


def replace_ontology_term(dataframe, ontology_name, update_map):
    column_name = f"{ontology_name}_ontology_term_id"
    if dataframe[column_name].dtype != "category":
        dataframe[column_name] = dataframe[column_name].astype("category")
    for old_term, new_term in update_map.items():
        if old_term in dataframe[column_name].cat.categories:
            # add new one if not already in category, else continue
            if new_term not in dataframe[column_name].cat.categories:
                dataframe[column_name] = dataframe[column_name].cat.add_categories(new_term)
            # replace in dataset
            dataframe.loc[dataframe[column_name] == old_term, column_name] = new_term
            # remove deprecated_term from category
            dataframe[column_name] = dataframe[column_name].cat.remove_categories(old_term)


def remove_deprecated_features(adata: ad.AnnData) -> ad.AnnData:
    deprecated = get_deprecated_features()

    # Filter out genes that don't appear in the approved annotation
    var_to_keep = adata.var.index[~adata.var.index.isin(deprecated)].tolist()
    adata = adata[:, var_to_keep]

    # Repeat much of the same steps for the raw.var, if it exists
    if adata.raw:
        raw_adata = ad.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs)
        var_to_keep = raw_adata.var.index[~raw_adata.var.index.isin(deprecated)].tolist()
        raw_adata = raw_adata[:, var_to_keep]
        adata.raw = raw_adata
    return adata


def get_deprecated_features() -> List[str]:
    # return a list of deprecated feature ids.
    diff_map = []
    suffix = "_diff.txt"
    files = os.listdir(env.ONTOLOGY_DIR)
    for file in files:
        if file.endswith(suffix):
            with open(f"{env.ONTOLOGY_DIR}/{file}") as fp:
                lines = fp.read().splitlines()
                diff_map.extend(lines)
    return diff_map
