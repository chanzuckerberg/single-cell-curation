import sys
import anndata
import os
import yaml
import pandas as pd
import re
from typing import List, Dict, Tuple
from . import ontology
from . import env

ONTOLOGY_CHECKER = ontology.OntologyChecker()


def _is_null(v):

    """Return True if v is null, for one of the multiple ways a "null" value shows up in an h5ad."""
    return pd.isnull(v) or (hasattr(v, "__len__") and len(v) == 0)


def _get_schema_definition(version: str) -> dict:

    """
    Look up and read a schema definition based on a version number like "2.0.0".

    :param str version: Schema version

    :return The schema definition
    :rtype dict
    """

    path = os.path.join(env.SCHEMA_DEFINITIONS_DIR, version.replace(".", "_") + ".yaml")

    if not os.path.isfile(path):
        raise ValueError(f"No definition for version '{version}' found.")

    return yaml.load(open(path), Loader=yaml.FullLoader)


def _curie_remove_suffix(term_id: str, suffix_def: dict) -> Tuple[str, str]:
    """
    Remove suffix from a curie term id, if none present return it unmodified

    :param str term_id: the curie term id to validate
    :param dict{str: list[str], ...} suffix_def: dictionary whose keys are ontology term ids and values
    are list of allowed suffixes

    :rtype Tuple[str, str]
    :return the term_id with suffixed stripped, and the suffix
    """

    id_suffix = ""

    for ontology_name, suffixes in suffix_def.items():

        for suffix in suffixes:
            suffix = suffix.replace("(", "\(")
            suffix = suffix.replace(")", "\)")
            search_results = re.search(r"%s$" % suffix, term_id)
            if search_results:
                stripped_term_id = re.sub(r"%s$" % suffix, "", term_id)
                if ONTOLOGY_CHECKER.is_valid_term_id(ontology_name, stripped_term_id):
                    id_suffix = search_results.group(0)

                    return stripped_term_id, id_suffix

    return term_id, id_suffix


class Validator:
    """Handle validation of AnnData"""
    schema_definitions_dir = env.SCHEMA_DEFINITIONS_DIR

    def __init__(self):

        # Set initial state
        self.errors = []
        self.is_valid = False
        self.adata = anndata.AnnData()
        self.schema_def = dict()
        self.h5ad_path = ""



    def _validate_curie_allowed_terms(self, term_id: str, column_name: str, terms: Dict[str, List[str]]):

        """
        Validate a single curie term id is a valid children of any of allowed terms
        If there are any errors, it adds them to self.errors

        :param str term_id: the curie term id to validate
        :param str column_name: original column name in adata where the term_id comes from (used for error messages)
        :param dict{str: list[str]} terms: keys must be ontology names and values must lists of allowed terms

        :rtype None
        """

        checks = []
        for ontology, allowed_terms in terms.items():
            if ONTOLOGY_CHECKER.is_valid_term_id(ontology, term_id):
                checks.append(term_id in allowed_terms)

        if sum(checks) == 0 and len(checks) > 0:
            all_allowed = list(terms.values())
            self.errors.append(f"'{term_id}' in '{column_name}' is not an allowed term of '{all_allowed}'")

    def _validate_curie_ancestors(self, term_id: str, column_name: str, allowed_ancestors: Dict[str, List[str]]):

        """
        Validate a single curie term id is a valid children of any of allowed ancestors
        If there are any errors, it adds them to self.errors

        :param str term_id: the curie term id to validate
        :param str column_name: original column name in adata where the term_id comes from (used for error messages)
        :param dict{str: list[str]} allowed_ancestors, keys must be ontology names and values must lists of
        allowed ancestors

        :rtype None
        """

        checks = []

        for ontology, ancestors in allowed_ancestors.items():
            for ancestor in ancestors:
                if ONTOLOGY_CHECKER.is_valid_term_id(ontology, term_id) & ONTOLOGY_CHECKER.is_valid_term_id(ontology,
                                                                                                            ancestor):
                    checks.append(ONTOLOGY_CHECKER.is_descendent_of(ontology, term_id, ancestor))

        if sum(checks) == 0 and len(checks) > 0:
            all_ancestors = list(allowed_ancestors.values())
            self.errors.append(f"'{term_id}' in '{column_name}' is not a children term id of '{all_ancestors}'")

    def _validate_curie_ontology(self, term_id: str, column_name: str, allowed_ontologies: List[str]):

        """
        Validate a single curie term id belongs to specified ontologies
        If there are any errors, it adds them to self.errors

        :param str term_id: the curie term id to validate
        :param str column_name: original column name in adata where the term_id comes from (used for error messages)
        :param List[str] allowed_ontologies: allowed ontologies

        :rtype None
        """

        checks = []

        for ontology in allowed_ontologies:
            checks.append(ONTOLOGY_CHECKER.is_valid_term_id(ontology, term_id))

        if sum(checks) == 0:
            self.errors.append(
                f"'{term_id}' in '{column_name}' is not a valid ontology term id of '{', '.join(allowed_ontologies)}'"
            )


    def _validate_curie(self, term_id: str, column_name: str, curie_constraints: dict):

        """
        Validate a single curie term id based on some constraints.
        If there are any errors, it adds them to self.errors

        :param str term_id: the curie term id to validate
        :param str column_name: Name of the column in the dataframe
        :param dict curie_constraints: constraints for the curie term to be validated,
        this part of the schema definition

        :rtype None
        """

        # If there are exceptions and this is one then skip to end
        if "exceptions" in curie_constraints:
            if term_id in curie_constraints["exceptions"]:
                return None

        # Check if there are any allowed suffixes and remove them if needed
        if "suffixes" in curie_constraints:
            term_id, suffix = _curie_remove_suffix(term_id, curie_constraints["suffixes"])

        # Check that term id belongs to allowed ontologies
        self._validate_curie_ontology(term_id, column_name, curie_constraints["ontologies"])

        # If there are specified ancestors then make sure that this id is a valid child
        if "ancestors" in curie_constraints:
            self._validate_curie_ancestors(term_id, column_name, curie_constraints["ancestors"])

        # If there is a set of allowed terms check for it
        if "allowed_terms" in curie_constraints:
            self._validate_curie_allowed_terms(term_id, column_name, curie_constraints["allowed_terms"])

    def _validate_column(self, column: pd.Categorical, column_name: str, df_name: str, column_def: dict):

        """
        Given a schema definition and the column of a dataframe, verify that the column satisfies the schema.
        If there are any errors, it adds them to self.errors

        :param pandas.Categorical column: Column of a dataframe to validate
        :param str column_name: Name of the column in the dataframe
        :param str df_name: Name of the dataframe
        :param dict column_def: schema definition for this specific column,
        e.g. schema_def["obs"]["columns"]["cell_type_ontology_term_id"]

        :return A list of error messages. If that list is empty, the object passed validation.
        :rtype list
        """

        if column_def.get("unique"):
            if column.nunique() != len(column):
                self.errors.append(f"Column '{column_name}' in dataframe '{df_name}' is not unique.")

        if "enum" in column_def:
            bad_enums = [v for v in column if v not in column_def["enum"]]
            if bad_enums:
                self.errors.append(
                    f"Column '{column_name}' in dataframe '{df_name}' contains invalid values like "
                    f"'{bad_enums[0]}'. Values must be one of {column_def['enum']}."
                )

        if column_def.get("type") == "bool":
            if not column.dtype == bool:
                self.errors.append(
                    f"Column '{column_name}' in dataframe '{df_name}' must be boolean not {column.dtype.name}")

        if column_def.get("type") == "curie":
            if "curie_constraints" not in column_def:
                raise ValueError(f"Corrupt schema definition, no 'curie_constraints' were found for '{column_name}'")
            if "ontologies" not in column_def["curie_constraints"]:
                raise ValueError(f"allowed 'ontolgies' must be specified under 'curie constraints' for '{column_name}'")

            for curie in column.drop_duplicates():
                self._validate_curie(curie, column_name, column_def["curie_constraints"])

    def _validate_dataframe(self, df_name: str, component_def: dict):

        """
        Verifies the dataframe follows the schema. Adds errors to self.errors if any

        :param str df_name: Name of dataframe in the adata (e.g. "obs")
        :param dict component_def: component definition for this dataframe (from schema definition); e.g. schema_defintion["obs"]

        :rtype None
        """

        df = getattr(self.adata, df_name)

        if "index" in component_def:
            self._validate_column(df.index, "index", df_name, component_def["index"])

        for column in component_def.get("columns", []):
            if column not in df.columns:
                self.errors.append(f"Dataframe '{df_name}' is missing column '{column}'.")
            else:
                self._validate_column(
                    df[column], column, df_name, component_def["columns"][column]
                )

    def _deep_check(self):

        """
        Perform a "deep" check of the AnnData object using the schema definition. Adds errors to self.errors if any

        :rtype None
        """

        if "schema_version" not in self.adata.uns_keys():
            print("AnnData file is missing cellxgene's schema version")
            return

        for component, component_def in self.schema_def["components"].items():
            if component_def["type"] == "dataframe":
                self._validate_dataframe(component, component_def)
            elif component_def["type"] == "dict":
                # Placeholder
                self.errors.extend(
                    ""
                    # _validate_dict(getattr(adata, component), component, component_def)
                )
            else:
                raise ValueError(f"Unexpected component type '{component['type']}'")

    def _read_h5ad(self, h5ad_path: str):

        """Reads h5ad into self.adata"""
        try:
            self.adata = anndata.read_h5ad(h5ad_path, backed="r")
        except (OSError, TypeError):
            print(f"Unable to open '{h5ad_path}' with AnnData")
            sys.exit(1)

        # Set schema definition and path to h5ad
        self.schema_def = _get_schema_definition(self.adata.uns["schema_version"])
        self.h5ad_path = h5ad_path

    def validate_adata(self, h5ad_path) -> bool:

        """
        Validates adata

        :return True if successful validation, False otherwise
        :rtype bool
        """

        self._read_h5ad(h5ad_path)
        self._deep_check()

        if self.errors:
            print(*self.errors, sep="\n")
            self.is_valid = False
        else:
            self.is_valid = True

        return self.is_valid


class LabelWriter:
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
            raise ValueError("AnnData object is not valid or hasn't been run through validation. " 
                             "Validate AnnData first before attempting to write labels")

        # Read adata not in "full" mode (not read-only)
        self.adata = anndata.read_h5ad(validator.h5ad_path)

        self.schema_def = validator.schema_def
        self.errors = []
        self.was_writing_successful = False

    def _get_labels(self, component: str, column: str, column_definition: dict) -> pd.Categorical:

        """
        Retrieves a new column (pandas categorical) with labels based on the IDs in 'column' and the logic in the
        'column_definition'

        :param str component: what dataframe in self.adata to work with
        :param str column: Column in self.adata with IDs that will be used to retrieve values
        :param dict column_def: schema definition for this specific column,
        e.g. schema_def["obs"]["columns"]["cell_type_ontology_term_id"]

        :rtype pandas.Categorical
        :return new pandas column with labels corresponding to input column
        """

        type_labels = column_definition["add_labels"]["type"]
        current_df = getattr(self.adata, component)
        ids = getattr(current_df, column).drop_duplicates().tolist()

        if type_labels == "curie":

            if "curie_constraints" not in column_definition:
                raise ValueError(f"Schema definition error: 'add_lables' with type 'curie' was found for '{column}' "
                                 "but no curie constraints were found for the lables")

            mapping_dict = _get_mapping_dict_curie(ids=ids, curie_constraints=column_definition["curie_constraints"])

        else:
            raise TypeError(f"'{type_labels}' is not supported in 'add-labels' functionality")

        new_column = getattr(current_df, column).copy().replace(mapping_dict).astype("category")

        return new_column

    def _add_labels(self):

        """
        From a valid (per cellxgene's schema) adata, this function adds to self.adata ontology/gene labels
        to adata.obs and adata.var respectively
        """

        # First adata.obs
        for column, column_def in self.schema_def["components"]["obs"]["columns"].items():
            if "add_labels" in column_def:
                new_column = self._get_labels("obs", column, column_def)
                self.adata.obs[column_def["add_labels"]["to"]] = new_column

        # Second adata.var

    def _check_column_availability(self):

        """
        This method will add error messages to self.errors if reserved columns in self.adata.obs or
        self.adata.var already exist

        :rtype none
        """

        for component in ["obs"]:
            for column, columns_def in self.schema_def["components"]["obs"]["columns"].items():
                if "add_labels" in columns_def:
                    reserved_name = columns_def["add_labels"]["to"]
                    if reserved_name in getattr(self.adata, component):
                        self.errors.append(f"Add labels error: Column '{reserved_name}' is a reserved column name " 
                                           "of 'obs'. Remove it from h5ad and try again.")

    def write_labels(self, add_labels_file: str):

        """
        From a valid (per cellxgene's schema) h5ad, this function writes a new h5ad file with ontology/gene labels added
        to adata.obs  and adata.var respectively

        :param str add_labels_file: Path to new h5ad file with ontology/gene labels added

        :rtype None
        """

        # First check that columns to be created don't exist. Terminate process if errors found
        self._check_column_availability()
        if self.errors:
            print(*self.errors, sep="\n")
            self.was_writing_successful = False
            return

        # Add labels in obs
        self._add_labels()

        # Write file
        self.adata.write_h5ad(add_labels_file, compression="gzip")
        self.was_writing_successful = True


def _get_mapping_dict_curie(ids: List[str], curie_constraints: dict) -> Dict[str, str]:

    """
    From defined constraints it creates a mapping dictionary of ontology IDs and labels. Used internally
    by LabelWriter class

    :param list[str] ids: Ontology IDs use for mapping
    :param list[str] curie_constraints: curie constraints e.g. schema_def["components"]["obs"]["columns"]["cell_type_ontology_term_id"]["curie_"]

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
            ids[i], id_suffixes[i] = _curie_remove_suffix(ids[i], curie_constraints["suffixes"])

    for original_id, id, id_suffix in zip(original_ids, ids, id_suffixes):
        for ontology in allowed_ontologies:
            if ONTOLOGY_CHECKER.is_valid_term_id(ontology, id):
                mapping_dict[original_id] = ONTOLOGY_CHECKER.get_term_label(ontology, id) + id_suffix

        # If there are exceptions the label should be the same as the id
        if "exceptions" in curie_constraints:
            if original_id in curie_constraints["exceptions"]:
                mapping_dict[original_id] = original_id

    # Check that all ids got a mapping. All ids should be found if adata was validated
    for id in original_ids:
        if id not in mapping_dict:
            raise ValueError(f"Add labels error: Unable to get label for '{id}'")

    return mapping_dict


def validate(h5ad_path: str, add_labels_file: str = None):

    """
    Entry point for validation.

    :param str h5ad_path: Path to h5ad file to validate
    :param str add_labels_file: Path to new h5ad file with ontology/gene labels added

    :return True if successful validation, False otherwise
    :rtype bool
    """

    # Perform validation
    validator = Validator()
    validator.validate_adata(h5ad_path)

    # Stop if validation was unsuccessful
    if not validator.is_valid:
        return False

    success = False
    if add_labels_file:
        writer = LabelWriter(validator)
        writer.write_labels(add_labels_file)

        success = validator.is_valid & writer.was_writing_successful

    return success
