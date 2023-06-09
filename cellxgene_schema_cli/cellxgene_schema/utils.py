import os
from typing import List

import env
import pandas


def replace_ontology_term(dataframe: pandas.DataFrame, ontology_name: str, update_map: dict):
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


def remove_deprecated_features(dataframe: pandas.DataFrame):
    deprecated = get_deprecated_features()
    var_in_deprecated = dataframe.var.index[~dataframe.var.index.isin(deprecated)].tolist()
    var_to_keep = dataframe.var.index.tolist()
    var_to_keep = [e for e in var_to_keep if e not in var_in_deprecated]
    dataset = dataframe[:, var_to_keep]
    return dataset


def get_deprecated_features() -> List[str]:
    # return a list of deprecated feature ids.
    diff_map = []
    suffix = "_diff.txt"
    files = os.listdir(env.ONTOLOGY_DIR)
    for file in files:
        if file.endswith(suffix):
            with open(f"{env.ONTOLOGY_DIR}/{file}") as fp:
                lines = fp.readlines()
                diff_map.extend(lines)
    return diff_map
