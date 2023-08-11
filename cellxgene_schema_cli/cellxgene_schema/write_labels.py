import logging
import traceback
from typing import Dict, List, Optional

import pandas as pd

from cellxgene_schema import ontology
from cellxgene_schema.validate import ONTOLOGY_CHECKER, Validator

from .utils import getattr_anndata

logger = logging.getLogger(__name__)


class AnnDataLabelAppender:
    """
    From valid h5ad, handles writing a new h5ad file with ontology/gene labels added
    to adata.obs and adata.var respectively as indicated in the schema definition
    """

    def __init__(self, validator: Validator):
        """
        From a list of ids and defined constraints, creates a mapping dictionary {id: label, ...}

        :param Validator validator: a Validator object, it's used to get adata and schema defintion for its validation,
        it's also used to make sure the validation on this adata was successful.
        """

        if not validator.is_valid:
            raise ValueError(
                "AnnData object is not valid or hasn't been run through validation. "
                "Validate AnnData first before attempting to write labels"
            )

        if validator.adata.isbacked:
            try:
                self.adata = validator.adata.to_memory()
            except ValueError:
                self.adata = validator.adata
        else:
            self.adata = validator.adata.copy()
        self.validator = validator
        self.schema_def = validator.schema_def
        self.errors = []
        self.was_writing_successful = False

    def _merge_dicts(self, dict1: dict, dict2: dict) -> dict:
        """
        Recursively merges two dicts, designed to be used to flatten a column definition.

        :params dict dict1: first dict
        :params dict dict2: second dict

        :rtype dict
        :return the merged dict
        """

        merged_dict = dict1.copy()

        for key, value_2 in dict2.items():
            if key == "rule":
                continue

            if key not in dict1:
                merged_dict[key] = value_2
            else:
                value_1 = dict1[key]

                if not type(value_1) == type(value_2):
                    raise ValueError("Inconsistent types, impossible to merge")

                if isinstance(value_2, str):
                    if key == "error_message_suffix":
                        merged_dict[key] = value_1 + " " + value_2
                        continue
                    if not value_2 == value_1:
                        raise ValueError(f"Strings types in dependencies cannot be different, {value_1} and {value_2}")

                elif isinstance(value_2, list):
                    merged_dict[key] = list(set(value_1 + value_2))
                elif isinstance(value_2, dict):
                    merged_dict[key] = self._merge_dicts(value_1, value_2)
                else:
                    raise ValueError(f"merging {type(value_2)} is not implemented")

        return merged_dict

    def _flatten_column_def_with_dependencies(self, column_def: dict) -> dict:
        """
        Flattens a column definition  that has dependencies, it essentially concatenates all the definitions of the
        dependencies into the on definition

        :params dict column_def: the column definition that has dependencies in it

        :rtype dict
        :return the flatten column definition
        """

        # Do nothing if ther are no dependencies
        if "dependencies" not in column_def:
            return column_def

        flatten = column_def.copy()
        del flatten["dependencies"]

        for dep in column_def["dependencies"]:
            flatten = self._merge_dicts(flatten, dep)

        return flatten

    def _get_mapping_dict_curie(self, ids: List[str], curie_constraints: dict) -> Dict[str, str]:
        """
        From defined constraints it creates a mapping dictionary of ontology IDs and labels.

        :param list[str] ids: Ontology IDs use for mapping
        :param list[str] curie_constraints: curie constraints e.g.
        schema_def["components"]["obs"]["columns"]["cell_type_ontology_term_id"]["curie_"]

        :return a mapping dictionary: {id: label, ...}
        :rtype dict
        """

        mapping_dict = {}
        allowed_ontologies = curie_constraints["ontologies"]

        # Remove any suffixes if any
        # original_ids will have untouched ids which will be used for mapping
        # id_suffixes will save suffixes if any, these will be used to append to labels
        # ids will have the ids without suffixes
        original_ids = ids.copy()
        id_suffixes = [""] * len(ids)

        if "suffixes" in curie_constraints:
            for i in range(len(ids)):
                ids[i], id_suffixes[i] = Validator._curie_remove_suffix(ids[i], curie_constraints["suffixes"])

        for original_id, id, id_suffix in zip(original_ids, ids, id_suffixes):
            # If there are exceptions the label should be the same as the id
            if "exceptions" in curie_constraints and original_id in curie_constraints["exceptions"]:
                mapping_dict[original_id] = original_id
                continue

            for ontology_name in allowed_ontologies:
                if ontology_name == "NA":
                    continue
                if ONTOLOGY_CHECKER.is_valid_term_id(ontology_name, id):
                    mapping_dict[original_id] = ONTOLOGY_CHECKER.get_term_label(ontology_name, id) + id_suffix

        # Check that all ids got a mapping. All ids should be found if adata was validated
        for id in original_ids:
            if id not in mapping_dict:
                raise ValueError(f"Add labels error: Unable to get label for '{id}'")

        return mapping_dict

    def _get_mapping_dict_feature_id(self, ids: List[str]) -> Dict[str, str]:
        """
        Creates a mapping dictionary of gene/feature IDs and labels.

        :param list[str] ids: Gene/feature IDs use for mapping

        :return a mapping dictionary: {id: label, ...}
        :rtype dict
        """

        mapping_dict = {}

        for i in ids:
            organism = ontology.get_organism_from_feature_id(i)
            mapping_dict[i] = self.validator.gene_checkers[organism].get_symbol(i)

        return mapping_dict

    def _get_mapping_dict_feature_reference(self, ids: List[str]) -> Dict[str, Optional[ontology.SupportedOrganisms]]:
        """
        Creates a mapping dictionary of gene/feature IDs and NCBITaxon curies

        :param list[str] ids: Gene/feature IDs use for mapping

        :return a mapping dictionary: {id: label, ...}
        :rtype dict
        """

        mapping_dict = {}

        for i in ids:
            organism = ontology.get_organism_from_feature_id(i)
            mapping_dict[i] = organism.value

        return mapping_dict

    def _get_mapping_dict_feature_biotype(self, ids: List[str]) -> Dict[str, str]:
        """
        Creates a mapping dictionary of feature IDs and biotype ("gene" or "spike-in")

        :param list[str] ids: feature IDs use for mapping

        :return a mapping dictionary: {id: "gene", id: "spike-in", ...}
        :rtype dict
        """
        mapping_dict = {}

        for i in ids:
            if i.startswith("ERCC"):
                mapping_dict[i] = "spike-in"
            elif i.startswith("ENS"):
                mapping_dict[i] = "gene"
            else:
                raise ValueError(f"{i} is not a recognized `feature_name` and cannot be assigned a `feature_type`")

        return mapping_dict

    def _get_labels(
        self,
        component: str,
        column: str,
        column_definition: dict,
        label_type: dict,
    ) -> pd.Categorical:
        """
        Retrieves a new column (pandas categorical) with labels based on the IDs in 'column' and the logic in the
        'column_definition'

        :param str component: what dataframe in self.adata to work with
        :param str column: Column in self.adata with IDs that will be used to retrieve values
        :param dict column_definition: schema definition of the column
        e.g. schema_def["obs"]["columns"]["cell_type_ontology_term_id"]
        :param dict label_type: the type of label

        :rtype pandas.Categorical
        :return new pandas column with labels corresponding to input column
        """

        # Set variables for readability
        current_df = getattr_anndata(self.adata, component)

        if column == "index":
            original_column = pd.Series(current_df.index)
            original_column.index = current_df.index
        else:
            original_column = getattr(current_df, column)

        ids = getattr(current_df, column).drop_duplicates().tolist()

        # Flatten column definition (will do so if there are dependencies in the definition
        column_definition = self._flatten_column_def_with_dependencies(column_definition)

        if label_type == "curie":
            if "curie_constraints" not in column_definition:
                raise ValueError(
                    f"Schema definition error: 'add_labels' with type 'curie' was found for '{column}' "
                    "but no curie constraints were found for the labels"
                )

            mapping_dict = self._get_mapping_dict_curie(ids, column_definition["curie_constraints"])

        elif label_type == "feature_id":
            mapping_dict = self._get_mapping_dict_feature_id(ids=ids)

        elif label_type == "feature_reference":
            mapping_dict = self._get_mapping_dict_feature_reference(ids=ids)

        elif label_type == "feature_biotype":
            mapping_dict = self._get_mapping_dict_feature_biotype(ids=ids)

        else:
            raise TypeError(f"'{label_type}' is not supported in 'add-labels' functionality")

        new_column = original_column.copy().replace(mapping_dict).astype("category")

        return new_column

    def _add_column(self, component: str, column: str, column_definition: dict):
        """
        Adds a new column (pandas categorical) to a component of adata with labels based on the IDs
        in 'column' and the logic in the 'column_def'

        :param str component: what dataframe in self.adata to work with
        :param str column: Column in self.adata with IDs that will be used to retrieve values
        :param dict column_definition: schema definition of the column
        e.g. schema_def["obs"]["columns"]["cell_type_ontology_term_id"]

        :rtype None
        """

        for label_def in column_definition["add_labels"]:
            new_column = self._get_labels(component, column, column_definition, label_def["type"])
            new_column_name = label_def["to_column"]

            # The sintax below is a programtic way to access obs and var in adata:
            # adata.__dict__["_obs"] is adata.obs
            # "raw.var" requires to levels of programtic access
            if "." in component:
                [first_elem, second_elem] = component.split(".")
                self.adata.__dict__["_" + first_elem].__dict__["_" + second_elem][new_column_name] = new_column
            else:
                self.adata.__dict__["_" + component][new_column_name] = new_column

    def _add_labels(self):
        """
        From a valid (per cellxgene's schema) adata, this function adds to self.adata ontology/gene labels
        to adata.obs, adata.var, and adata.raw.var respectively
        """
        for component in ["obs", "var", "raw.var"]:
            # If the component does not exist, skip (this is for raw.var)
            if getattr_anndata(self.adata, component) is None:
                continue

            # Doing it for columns
            if "columns" in self.schema_def["components"][component]:
                for column, column_def in self.schema_def["components"][component]["columns"].items():
                    if "add_labels" in column_def:
                        self._add_column(component, column, column_def)

            # Doing it for index
            index_def = self.schema_def["components"][component]["index"]
            if "add_labels" in index_def:
                self._add_column(component, "index", index_def)

    def _remove_categories_with_zero_values(self):
        df = self.adata.obs
        for column in df.columns:
            col = df[column]
            if col.dtype == "category":
                df[column] = col.cat.remove_unused_categories()

    def write_labels(self, add_labels_file: str):
        """
        From a valid (per cellxgene's schema) h5ad, this function writes a new h5ad file with ontology/gene labels added
        to adata.obs  and adata.var respectively

        :param str add_labels_file: Path to new h5ad file with ontology/gene labels added

        :rtype None
        """
        logger.info("Writing labels")
        # Add labels in obs
        self._add_labels()

        # Remove unused categories
        self._remove_categories_with_zero_values()

        # Update version
        self.adata.uns["schema_version"] = self.validator.schema_version

        # Write file
        try:
            self.adata.write_h5ad(add_labels_file, compression="gzip")
        except Exception as e:
            tb = traceback.format_exc()
            self.errors.append((f"Writing h5ad was unsuccessful, got exception '{e}'.", tb))

        # Print errors if any
        if self.errors:
            for e, tb in self.errors:
                logger.error(e, extra={"exec_info": tb})
            self.was_writing_successful = False
        else:
            self.was_writing_successful = True
