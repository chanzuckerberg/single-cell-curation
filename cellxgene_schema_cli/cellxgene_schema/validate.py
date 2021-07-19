import sys
import anndata
import os
import yaml
import pandas as pd
import re
from . import ontology

SCHEMA_DEFINITIONS_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "schema_definitions")
ONTOLOGY_CHECKER = ontology.ontologyChecker()

def _validate_curie_allowed_terms(term_id, column_name, terms):

    """
    Validate a single curie term id is a valid children of any of allowed terms

    :param str term_id: the curie term id to validate
    :param str column_name: original column name in adata where the term_id comes from (used for error messages)
    :param dict{str: list[str]} allowed_terms, keys must be ontology names and values must lists of allowed terms

    :return list with error messages, empty if none
    :rtype list[str|none]
    """

    errors = []
    checks = []

    for ontology, allowed_terms in terms.items():
        if ONTOLOGY_CHECKER.is_valid_term_id(ontology, term_id):
            checks.append(term_id in allowed_terms)

    if sum(checks) == 0 and len(checks) > 0:
        all_allowed = list(terms.values())
        errors.append(f"'{term_id}' in '{column_name}' is not an allowed term of '{all_allowed}'")

    return errors

def _validate_curie_ancestors(term_id, column_name, allowed_ancestors):
    """
    Validate a single curie term id is a valid children of any of allowed ancestors

    :param str term_id: the curie term id to validate
    :param str column_name: original column name in adata where the term_id comes from (used for error messages)
    :param dict{str: list[str]} allowed_ancestors, keys must be ontology names and values must lists of allowed ancestors

    :return list with error messages, empty if none
    :rtype list[str|none]
    """

    checks = []
    errors = []

    for ontology, ancestors in allowed_ancestors.items():
        for ancestor in ancestors:
            if ONTOLOGY_CHECKER.is_valid_term_id(ontology, term_id) & ONTOLOGY_CHECKER.is_valid_term_id(ontology, ancestor):
                checks.append(ONTOLOGY_CHECKER.is_descendent_of(ontology, term_id, ancestor))

    if sum(checks) == 0 and len(checks) > 0:
        all_ancestors = list(allowed_ancestors.values())
        errors.append(f"'{term_id}' in '{column_name}' is not a children term id of '{all_ancestors}'")

    return errors


def _validate_curie_ontology(term_id, column_name, allowed_ontologies):

    """
    Validate a single curie term id belongs to specified ontologies

    :param str term_id: the curie term id to validate
    :param str column_name: original column name in adata where the term_id comes from (used for error messages)
    :param list[str] allowed_ontologies: allowed ontologies

    :return list with error messages, empty if none
    :rtype list[str|none]
    """

    errors = []
    checks=[]

    for ontology in allowed_ontologies:
        checks.append(ONTOLOGY_CHECKER.is_valid_term_id(ontology, term_id))

    if sum(checks) == 0:
        errors.append(f"'{term_id}' in '{column_name}' is not a valid ontology term id of '{', '.join(allowed_ontologies)}'")

    return errors

def _validate_curie_remove_suffix(term_id, suffix_def):

    """
    Remove suffix from a term id, if none present return it unmodified

    :param str term_id: the curie term id to validate
    :param dict{str: list[str], ...} suffix_def: dictionary whose keys are ontologies and values are list of allowed
    suffixes

    :return the term_id with suffixed stripped
    :rtype str
    """

    id_suffix = ""

    for ontology, suffixes in suffix_def.items():

        for suffix in suffixes:
            suffix = suffix.replace("(", "\(")
            suffix = suffix.replace(")", "\)")
            search_results = re.search(r"%s$" % suffix, term_id)
            if search_results:
                stripped_term_id = re.sub(r"%s$" % suffix, "", term_id)
                if ONTOLOGY_CHECKER.is_valid_term_id(ontology, stripped_term_id):
                    id_suffix = search_results.group(0)

                    return stripped_term_id, id_suffix

    return term_id, id_suffix

def _validate_curie(term_id, column_name, curie_constraints):

    """
    Validate a single curie term id based on some constraints, if invalid it will return x other wise none

    :param str term_id: the curie term id to validate
    :param dict curie_constraints: constraints for the curie term to be validated, this part of the schema definition

    :return list with error messages, empty if none
    :rtype list[str|none]
    """

    errors = []

    term_id_original = term_id

    # If there are exceptions and this is one then skip to end
    if "exceptions" in curie_constraints:
        if term_id in curie_constraints["exceptions"]:
            return errors

    # Check if there are any allowed suffixes and remove them if needed
    if "suffixes" in curie_constraints:
        term_id, suffix = _validate_curie_remove_suffix(term_id, curie_constraints["suffixes"])

    # Check that term id belongs to allowed ontologies
    errors.extend(_validate_curie_ontology(term_id, column_name, curie_constraints["ontologies"]))

    # If there are specified ancestors then make sure that this id is a valid child
    if "ancestors" in curie_constraints:
        errors.extend(_validate_curie_ancestors(term_id, column_name, curie_constraints["ancestors"]))

    # If there is a set of allowed terms check for it
    if "allowed_terms" in curie_constraints:
        errors.extend(_validate_curie_allowed_terms(term_id, column_name, curie_constraints["allowed_terms"]))

    return errors



def _is_null(v):
    """Return True if v is null, for one of the multiple ways a "null" value shows up in an h5ad."""
    return pd.isnull(v) or (hasattr(v, "__len__") and len(v) == 0)


def _validate_column(column, column_name, df_name, column_def):
    """Given a schema definition and the column of a dataframe, verify that the column satifies the schema

    :param pandas.DataFrame column: Column of a dataframe to validate
    :param str column_name: Name of the column in the dataframe
    :param str df_name: Name of the dataframe
    :param dict column_def: schema definition for this specific column, e.g. schema_def["obs"]["columns"]["cell_type_ontology_term_id"]

    :return A list of error messages. If that list is empty, the object passed validation.
    :rtype list
    """

    errors = []

    if column_def.get("unique"):
        if column.nunique() != len(column):
            errors.append(f"Column '{column_name}' in dataframe '{df_name}' is not unique.")

    if "enum" in column_def:
        bad_enums = [v for v in column if v not in column_def["enum"]]
        if bad_enums:
            errors.append(
                f"Column '{column_name}' in dataframe '{df_name}' contains unpermitted values like "
                f"'{bad_enums[0]}'. Values must be one of {column_def['enum']}."
            )

    if column_def.get("type") == "curie":
        if "curie_constraints" not in column_def:
            raise ValueError(f"Corrupt schema definition, no 'curie_constraints' were found for '{column_name}'")
            if "ontologies" not in column_def[curie_constraints]:
                raise ValueError(f"allowed 'ontolgies' must be specified under 'curie constraints' for '{column_name}'")

        for curie in column.drop_duplicates():
            errors.extend(_validate_curie(curie, column_name, column_def["curie_constraints"]))

    return errors


def _validate_dataframe(df, df_name, component_def):

    """Verifies the dataframe follows the schema.
    :param pandas.DataFrame df: Dataframe to validate
    :param str df_name: Name that this data frame has in the adata (e.g. "obs")
    :param dict component_def: component definition read for this dataframe; e.g. schema_defintion["obs"]

    :return A list of error messages. If that list is empty, the object passed validation.
    :rtype list
    """

    errors = []

    if "index" in component_def:
        errors.extend(_validate_column(df.index, "index", df_name, component_def["index"]))

    for column in component_def.get("columns", []):
        if column not in df.columns:
            errors.append(f"Dataframe '{df_name}' is missing column '{column}'.")
        else:
            errors.extend(
                _validate_column(
                    df[column], column, df_name, component_def["columns"][column]
                )
            )

    return errors


def deep_check(adata, schema_def):
    """Perform a "deep" check of the AnnData object using the schema definition.

    :param anndata.AnnData adata: AnnData object

    :return A list of error messages. If that list is empty, the object passed validation.
    :rtype list
    """

    errors = []

    for component, component_def in schema_def["components"].items():
        if component_def["type"] == "dataframe":
            errors.extend(
                _validate_dataframe(getattr(adata, component), component, component_def)
            )
        elif component_def["type"] == "dict":
            errors.extend(
                ""
                #_validate_dict(getattr(adata, component), component, component_def)
            )
        else:
            raise ValueError(f"Unexpected component type '{component['type']}'")

    return errors


def get_schema_definition(version):
    """
    Look up and read a schema definition based on a version number like "2.0.0".

    :param str version: Schema version

    :return The schema definition
    :rtype dict
    """

    path = os.path.join(SCHEMA_DEFINITIONS_PATH, version.replace(".", "_") + ".yaml")

    if not os.path.isfile(path):
        raise ValueError(f"No definition for version '{version}' found.")

    return yaml.load(open(path), Loader=yaml.FullLoader)


def validate_adata(adata, schema_def):

    """
    Validates adata

    :param anndata.AnnData adata: AnnData object
    :param dict schema_def: schema definition read with "get_schema_definition"

    :return True if successful validation, False otherwise
    :rtype bool
    """

    if "schema_version" not in adata.uns_keys():
        print("AnnData file is missing cellxgene's schema version")
        return False

    errors = deep_check(adata, schema_def)

    if errors:
        print(*errors, sep="\n")
        return False
    else:
        return True


def _get_mapping_dict_curie(ids, curie_constraints):
    """
    From a list of ids and defined constraints, creates a mapping dictionary {id: label, ...}

    :param list[str] ids: Ontology IDs use for mapping
    :param list[str] allowed_ontologies: List of allowed ontologies for conversion

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
            ids[i], id_suffixes[i] = _validate_curie_remove_suffix(ids[i], curie_constraints["suffixes"])

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


def _get_labels(adata, component, column, column_definition):
    """
    Retrieves a new column (pandas categorical) with labels based on the IDs in 'column' and the 'type lables'

    :param anndata.AnnData adata: A valid (per cellxgene's schema) adata
    :param str component: what data frame in adata to work with
    :param str column: Column in adata with IDs that will be used to retrieve values
    :param dict column_def: schema definition for this specific column, e.g. schema_def["obs"]["columns"]["cell_type_ontology_term_id"]

    :rtype pandas.Categorical
    """

    type_labels = column_definition["add_labels"]["type"]
    current_df = getattr(adata, component)
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


def _check_column_availability(adata, schema_def):
    """
    From a valid (per cellxgene's schema) h5ad, it will return a list with error messages if reserved columns
    in adata.obs or adata.var already exist

    :param anndata.AnnData adata: A valid (per cellxgene's schema) adata
    :param dict schema_def: schema definition read with "get_schema_definition"

    :return A list with errors, empty if none
    :rtype List
    """

    errors = []
    for component in ["obs"]:
        for column, columns_def in schema_def["components"]["obs"]["columns"].items():
            if "add_labels" in columns_def:
                reserved_name = columns_def["add_labels"]["to"]
                if reserved_name in getattr(adata, component):
                    errors.append(f"Add labels error: Column '{reserved_name}' is a reserved column name of 'obs'. "
                                  "Remove it from h5ad and try again.")

    return errors


def _add_labels(adata, schema_def):

    """
    From a valid (per cellxgene's schema) adata, this function returns a new adata with ontology/gene labels added
    to adata.obs and adata.var respectively

    :param anndata.AnnData adata: A valid (per cellxgene's schema) adata
    :param dict schema_def: schema definition read with "get_schema_definition"

    :return The original adata with labels added to adata.obs and adata.var
    :rtype anndata.AnnData
    """

    for column, column_def in schema_def["components"]["obs"]["columns"].items():
        if "add_labels" in column_def:
            new_column = _get_labels(adata, "obs", column, column_def)
            adata.obs[column_def["add_labels"]["to"]] = new_column

    return adata


def write_labels(adata, add_labels_file, schema_def):

    """
    From a valid (per cellxgene's schema) h5ad, this function writes a new h5ad file with ontology/gene labels added
    to adata.obs  and adata.var respectively

    :param anndata.AnnData adata: A valid (per cellxgene's schema) adata
    :param str add_labels_file: Path to new h5ad file with ontology/gene labels added
    :param dict schema_def: schema definition read with "get_schema_definition"

    :rtype None
    """

    errors = []

    # First check that columns to be created don't exist. Terminate process if errors found
    errors.extend(_check_column_availability(adata, schema_def))
    if errors:
        print(*errors)
        return False

    # Add labels in obs
    adata = _add_labels(adata, schema_def)

    # Write file
    adata.write_h5ad(add_labels_file, compression="gzip")

    return True

def validate(h5ad_path, add_labels_file=None):
    """
    Entry point for validation.

    :param str h5ad_path: Path to h5ad file to validate
    :param str add_labels_file: Path to new h5ad file with ontology/gene labels added

    :return True if successful validation, False otherwise
    :rtype bool
    """

    try:
        adata = anndata.read_h5ad(h5ad_path, backed="r")
    except (OSError, TypeError):
        print(f"Unable to open '{h5ad_path}' with AnnData")
        sys.exit(1)

    schema_def = get_schema_definition(adata.uns["schema_version"])

    # Perform validation
    is_validation_successful = validate_adata(adata, schema_def)

    if not is_validation_successful:
        return is_validation_successful
    else:
        # Add labels if indicated
        is_add_labels_successful = True
        if add_labels_file:
            is_add_labels_successful = write_labels(adata, add_labels_file, schema_def)

        return is_validation_successful & is_add_labels_successful
