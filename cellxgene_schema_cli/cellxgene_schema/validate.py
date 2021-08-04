import sys
import anndata
import os
import yaml
import pandas as pd
from numpy import ndarray
from numpy import count_nonzero
from numpy import cumprod
from scipy import sparse
import re
from typing import List, Dict, Tuple, Union
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


def _get_organism_from_feature_id (feature_id: str) -> Union[ontology.SupportedOrganisms, None]:

    """
    Infers the organism of a feature id based on the prefix of a feature id, e.g. ENSG means Homo sapiens

    :param str feature_id: the feature id

    :rtype Union[ontology.SypportedOrganisms, None]
    :return: the organism the feature id is from
    """

    if feature_id.startswith("ENSG") or feature_id.startswith("ENST"):
        return ontology.SupportedOrganisms.HOMO_SAPIENS
    elif feature_id.startswith("ENSMUS"):
        return ontology.SupportedOrganisms.MUS_MUSCULUS
    elif feature_id.startswith("ENSSAS"):
        return ontology.SupportedOrganisms.HOMO_SAPIENS
    elif feature_id.startswith("ERCC"):
        return ontology.SupportedOrganisms.ERCC
    else:
        return None

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
            suffix = suffix.replace("(", r"\(")
            suffix = suffix.replace(")", r"\)")
            search_results = re.search(r"%s$" % suffix, term_id)
            if search_results:
                stripped_term_id = re.sub(r"%s$" % suffix, "", term_id)
                if ONTOLOGY_CHECKER.is_valid_term_id(ontology_name, stripped_term_id):
                    id_suffix = search_results.group(0)

                    return stripped_term_id, id_suffix

    return term_id, id_suffix




class Validator:
    """Handles validation of AnnData"""

    schema_definitions_dir = env.SCHEMA_DEFINITIONS_DIR

    def __init__(self):

        # Set initial state
        self.errors = []
        self.warnings = []
        self.is_valid = False
        self.adata = anndata.AnnData()
        self.schema_def = dict()
        self.h5ad_path = ""
        self.invalid_feature_ids = []
        self.gene_checkers = dict() # Values will be instances of ontology.GeneChecker,
                                    # keys will be one of ontology.SupportedOrganisms

    def _read_h5ad(self, h5ad_path: Union[str, bytes, os.PathLike]):

        """
        Reads h5ad into self.adata
        :params Union[str, bytes, os.PathLike] h5ad_path: path to h5ad to read

        :rtype None
        """
        try:
            # H5AD has to be loaded in memory mode. If not the types of X are not properly retrieved by anndata
            # see https://github.com/theislab/anndata/issues/326#issuecomment-892203924
            self.adata = anndata.read_h5ad(h5ad_path, backed=None)
        except (OSError, TypeError):
            print(f"Unable to open '{h5ad_path}' with AnnData")
            sys.exit(1)

        self.h5ad_path = h5ad_path

    def _set_schema_def(self):
        """
        Sets schema dictionary from using information in adata
        rtype: None
        """

        if "schema_version" not in self.adata.uns:
            raise ValueError("adata has no schema definition in 'adata.uns'")

        self.schema_def = _get_schema_definition(self.adata.uns["schema_version"])

    def _get_component_def(self, component: str) -> dict:
        """
        Gets the definition of an individual component in the schema (e.g. obs)

        :param component: the component name

        :rtype dict
        """

        if self.schema_def:
            return self.schema_def["components"][component]
        else:
            raise RuntimeError("Schema has not been set in this instance class")

    def _get_column_def(self, component: str, column_name: str) -> dict:
        """
        Gets the definition of a column from a component in the schema (e.g. obs)

        :param component: the component name
        :param str column_name: the column name

        :rtype dict
        """

        if column_name == "index":
            return self._get_component_def(component)["index"]
        else:
            return self._get_component_def(component)["columns"][column_name]

    def _validate_curie_allowed_terms(
        self, term_id: str, column_name: str, terms: Dict[str, List[str]]
    ):

        """
        Validate a single curie term id is a valid children of any of allowed terms
        If there are any errors, it adds them to self.errors

        :param str term_id: the curie term id to validate
        :param str column_name: original column name in adata where the term_id comes from (used for error messages)
        :param dict{str: list[str]} terms: keys must be ontology names and values must lists of allowed terms

        :rtype None
        """

        checks = []
        for ontology_name, allowed_terms in terms.items():
            if ONTOLOGY_CHECKER.is_valid_term_id(ontology_name, term_id):
                checks.append(term_id in allowed_terms)

        if sum(checks) == 0 and len(checks) > 0:
            all_allowed = list(terms.values())
            self.errors.append(
                f"'{term_id}' in '{column_name}' is not an allowed term of '{all_allowed}'"
            )

    def _validate_curie_ancestors(
        self, term_id: str, column_name: str, allowed_ancestors: Dict[str, List[str]]
    ):

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

        for ontology_name, ancestors in allowed_ancestors.items():
            for ancestor in ancestors:
                if ONTOLOGY_CHECKER.is_valid_term_id(
                    ontology_name, term_id
                ) & ONTOLOGY_CHECKER.is_valid_term_id(ontology_name, ancestor):
                    checks.append(
                        ONTOLOGY_CHECKER.is_descendent_of(
                            ontology_name, term_id, ancestor
                        )
                    )

        if True not in checks:
            all_ancestors = list(allowed_ancestors.values())
            self.errors.append(
                f"'{term_id}' in '{column_name}' is not a children term id of '{all_ancestors}'"
            )

    def _validate_curie_ontology(
        self, term_id: str, column_name: str, allowed_ontologies: List[str]
    ):

        """
        Validate a single curie term id belongs to specified ontologies
        If there are any errors, it adds them to self.errors

        :param str term_id: the curie term id to validate
        :param str column_name: original column name in adata where the term_id comes from (used for error messages)
        :param List[str] allowed_ontologies: allowed ontologies

        :rtype None
        """

        checks = []

        for ontology_name in allowed_ontologies:
            checks.append(ONTOLOGY_CHECKER.is_valid_term_id(ontology_name, term_id))

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
                return

        # If NA is found in allowed ontologies, it means only exceptions should be found. If no exceptions were found
        # then return error
        if curie_constraints["ontologies"] == ["NA"]:
            self.errors.append(
                f"'{term_id}' in '{column_name}' is not a valid value of '{column_name}'"
            )
            return

        # Check if there are any allowed suffixes and remove them if needed
        if "suffixes" in curie_constraints:
            term_id, suffix = _curie_remove_suffix(
                term_id, curie_constraints["suffixes"]
            )

        # Check that term id belongs to allowed ontologies
        self._validate_curie_ontology(
            term_id, column_name, curie_constraints["ontologies"]
        )

        # If there are specified ancestors then make sure that this id is a valid child
        if "ancestors" in curie_constraints:
            self._validate_curie_ancestors(
                term_id, column_name, curie_constraints["ancestors"]
            )

        # If there is a set of allowed terms check for it
        if "allowed_terms" in curie_constraints:
            self._validate_curie_allowed_terms(
                term_id, column_name, curie_constraints["allowed_terms"]
            )

    def _validate_feature_id(self, feature_id: str):

        """
        Validates a feature id, i.e. checks that it's present in the reference
        If there are any errors, it adds them to self.errors and adds it to the list of invalid features

        :param str feature_id: the feature id to be validated

        :rtype none
        """

        organism = _get_organism_from_feature_id(feature_id)

        if not organism:
            self.errors.append(f"Could not infer organism from '{feature_id}', make sure it is a valid ID")
            return

        if organism not in self.gene_checkers:
            self.gene_checkers[organism] = ontology.GeneChecker(organism)

        if not self.gene_checkers[organism].is_valid_id(feature_id):
            self.invalid_feature_ids.append(feature_id)
            self.errors.append(f"'{feature_id}' is not a valid feature ID")

        return

    def _validate_column(
        self, column: pd.Series, column_name: str, df_name: str, column_def: dict
    ):

        """
        Given a schema definition and the column of a dataframe, verify that the column satisfies the schema.
        If there are any errors, it adds them to self.errors

        :param pandas.Series column: Column of a dataframe to validate
        :param str column_name: Name of the column in the dataframe
        :param str df_name: Name of the dataframe
        :param dict column_def: schema definition for this specific column,
        e.g. schema_def["obs"]["columns"]["cell_type_ontology_term_id"]

        :rtype None
        """

        if column_def.get("unique"):
            if column.nunique() != len(column):
                self.errors.append(
                    f"Column '{column_name}' in dataframe '{df_name}' is not unique."
                )

        if column_def.get("type") == "bool":
            if not column.dtype == bool:
                self.errors.append(
                    f"Column '{column_name}' in dataframe '{df_name}' must be boolean not {column.dtype.name}"
                )

        if column_def.get("type") == "categorical":
            if not column.dtype.name == "category":
                self.errors.append(
                    f"Column '{column_name}' in dataframe '{df_name}' must be categorical not {column.dtype.name}"
                )

        if "enum" in column_def:
            bad_enums = [v for v in column.drop_duplicates() if v not in column_def["enum"]]
            if bad_enums:
                self.errors.append(
                    f"Column '{column_name}' in dataframe '{df_name}' contains invalid values "
                    f"'{bad_enums}'. Values must be one of {column_def['enum']}."
                )

        if column_def.get("type") == "feature_id":
            for feature_id in column:
                self._validate_feature_id(feature_id)

        if column_def.get("type") == "curie":
            if "curie_constraints" not in column_def:
                raise ValueError(
                    f"Corrupt schema definition, no 'curie_constraints' were found for '{column_name}'"
                )
            if "ontologies" not in column_def["curie_constraints"]:
                raise ValueError(
                    f"allowed 'ontolgies' must be specified under 'curie constraints' for '{column_name}'"
                )

            for term_id in column.drop_duplicates():
                self._validate_curie(
                    term_id, column_name, column_def["curie_constraints"]
                )

    def _validate_column_dependencies(
        self, df: pd.DataFrame, df_name: str, column_name: str, dependencies: List[dict]
    ) -> pd.Series:

        """
        Validates subset of columns based on dependecies, for instance development_stage_ontology_term_id has
        dependencies with organism_ontology_term_id -- the allowed values depend on whether organism is human, mouse
        or something else.

        After performing all validations, it will return the column that has been stripped of the already validated
        fields -- this has to still be validated.

        :param pd.DataFrame df: pandas dataframe containing the column to be validated
        :param str df_name: the name of dataframe in the adata object, e.g. "obs"
        :param str column_name: the name of the column to be validated
        :param list dependencies: a list of dependecy definitions, which is a list of column definitions with a "rule"
        """

        all_rules = []

        for dependency_def in dependencies:
            column = getattr(df.query(dependency_def["rule"]), column_name)
            all_rules.append(dependency_def["rule"])

            self._validate_column(column, column_name, df_name, dependency_def)

        # Set column with the data that's left
        all_rules = " | ".join(all_rules)
        column = getattr(df.query("not (" + all_rules + " )"), column_name)

        return column

    def _validate_list(
        self, list_name: str, current_list: List[str], element_type: str
    ):

        """
        Validates the elements of a list based on the type definition. Adds errors to self.errors if any

        :param str list_name: name of list to use for error messages (if any)
        :param str current_list: the list to be validated
        :param str element_type: type to be validated

        :rtype None
        """

        for i in current_list:

            if element_type == "match_obs_columns":
                if i not in self.adata.obs.columns:
                    self.errors.append(
                        f"Value '{i}' of list '{list_name}' is not a column in 'adata.obs'"
                    )

    def _validate_dict(self, dictionary: dict, dict_name: str, dict_def: dict):

        """
        Verifies the dictionary follows the schema. Adds errors to self.errors if any

        :param str dictionary: The dictionary to validate
        :param str dict_name: Name of dictionary in the adata (e.g. "uns")
        :param str dict_def: The schema definition for this specific dictionary

        :rtype None
        """

        for key, value_def in dict_def["keys"].items():

            if key not in dictionary:
                if "required" in value_def:
                    self.errors.append(f"'{key}' in '{dict_name}' is not present")
                continue

            value = dictionary[key]

            if value_def["type"] == "string":
                if not isinstance(value, str):
                    self.errors.append(
                        f"'{value}' in '{dict_name}['{key}']' is not valid, it must be a string."
                    )
                    continue

                if "enum" in value_def:
                    if value not in value_def["enum"]:
                        self.errors.append(
                            f"'{value}' in '{dict_name}['{key}']' is not valid. "
                            f"Allowed terms: {value_def['enum']}."
                        )

            if value_def["type"] == "match_obsm_keys":
                if value not in self.adata.obsm:
                    self.errors.append(
                        f"'{value}' in '{dict_name}['{key}']' is not valid, "
                        f"it must be a key of 'adata.obsm'."
                    )

            if value_def["type"] == "list":
                if not (isinstance(value, list) or isinstance(value, ndarray)):
                    self.errors.append(
                        f"'{value}' in '{dict_name}['{key}']' is not valid, "
                        f"it must be a list or numpy array"
                    )
                    continue

                self._validate_list(key, value, value_def["element_type"])

    def _validate_dataframe(self, df_name: str):

        """
        Verifies the dataframe follows the schema. Adds errors to self.errors if any

        :param str df_name: Name of dataframe in the adata (e.g. "obs")

        :rtype None
        """

        df = getattr(self.adata, df_name)

        if "index" in self._get_component_def(df_name):
            self._validate_column(
                pd.Series(df.index), "index", df_name, self._get_column_def(df_name, "index")
            )

        for column_name in self._get_component_def(df_name)["columns"].keys():

            if column_name not in df.columns:
                self.errors.append(
                    f"Dataframe '{df_name}' is missing column '{column_name}'."
                )
                continue

            column_def = self._get_column_def(df_name, column_name)
            column = getattr(df, column_name)

            # First check if there are dependencies with other columns and wokr with a subset of the data if so
            if "dependencies" in column_def:
                column = self._validate_column_dependencies(
                    df, df_name, column_name, column_def["dependencies"]
                )

            # If after validating dependencies there's still values in the column, validate them.
            if len(column) > 0:
                self._validate_column(column, column_name, df_name, column_def)

    def _validate_sparsity(self):

        """
        calculates sparsity of x and raw.x, if bigger than indicated in the schema and not a scipy sparse matrix, then
        adds to warnings

        :rtype none
        """

        max_sparsity = self.schema_def["sparsity"]

        to_validate = [self.adata.X]
        to_validate_name = ["X"]

        # check if there's raw data
        if self.adata.raw:
            to_validate.append(self.adata.raw.X)
            to_validate_name.append("raw")

        # check if there's other expression matrices under layers
        if self.adata.layers:
            for key, value in self.adata.layers.items():
                to_validate.append(value)
                to_validate_name.append(f"layers['{key}']")

        # Check sparsity
        for x, x_name in zip(to_validate, to_validate_name):
            if not isinstance(x, sparse.csr_matrix):

                if isinstance(x, sparse.csc_matrix):
                    sparsity = 1 - x.count_nonzero() / float(cumprod(x.shape)[-1])
                elif isinstance(x, ndarray):
                    sparsity = 1 - count_nonzero(x) / float(cumprod(x.shape)[-1])
                else:
                    self.warnings.append(f"{x_name} matrix is of type {type(x)}, sparsity calculation hasn't been "
                                         f"implemented")
                    continue

                if sparsity > max_sparsity:
                    self.warnings.append(f"Sparsity of '{x_name}' is {sparsity} which is greater than {max_sparsity}, "
                                         f"and it is not a 'scipy.sparse.csr_matrix'. It is STRONGLY RECOMMENDED "
                                         f"to use this type of matrix for the given sparsity.")

    def _validate_embedding_dict(self):

        """
        Validates the embedding dictionary -- it checks that all values of adata.obms are numpy arrays with the correct
        dimension. Adds errors to self.errors if any.

        :rtype none
        """

        if not self.adata.obsm:
            self.errors.append("No embeddings found in 'adata.obsm'")
            return

        for key, value in self.adata.obsm.items():

            if not isinstance(value, ndarray):
                self.errors.append(f"All embeddings have to be of 'numpy.ndarray' type, "
                                   f"'adata.obsm['{key}']' is {type(value)}')")
                continue

            # obsm
            if value.shape[0] != self.adata.n_obs or value.shape[1] < 2:
                self.errors.append(f"All embeddings must have as many rows as cells, and at least two columns."
                                   f"'adata.obsm['{key}']' has shape of '{value.shape}'")
        return

    def _deep_check(self):

        """
        Perform a "deep" check of the AnnData object using the schema definition. Adds errors to self.errors if any

        :rtype None
        """

        # Checks that adata is fully loaded
        if self.adata.isbacked:
            raise RuntimeError("adata is loaded in backed mode, validation requires anndata to be loaded in memory")

        # Checks sparsity
        self._validate_sparsity()

        # Checks each component
        for component, component_def in self.schema_def["components"].items():
            if component_def["type"] == "dataframe":
                self._validate_dataframe(component)
            elif component_def["type"] == "dict":
                dictionary = getattr(self.adata, component)
                self._validate_dict(dictionary, component, component_def)
            elif component_def["type"] == "embedding_dict":
                self._validate_embedding_dict()
            else:
                raise ValueError(f"Unexpected component type '{component['type']}'")

    def validate_adata(self, h5ad_path: Union[str, bytes, os.PathLike] = None) -> bool:

        """
        Validates adata

        :params Union[str, bytes, os.PathLike] h5ad_path: path to h5ad to validate, if None it will try to validate
        from self.adata

        :return True if successful validation, False otherwise
        :rtype bool
        """

        # Re-start errors in case a new h5ad is being validated
        self.errors = []

        if h5ad_path:
            self._read_h5ad(h5ad_path)

        # Fetches schema def from anndata if schema version is not found in AnnData, this fails
        self._set_schema_def()

        self._deep_check()

        # Print warnings if any
        if self.warnings:
            self.warnings = ["WARNING: " + i for i in self.warnings]
            print(*self.warnings, sep="\n")

        # Print errors if any
        if self.errors:
            self.errors = ["ERROR: " + i for i in self.errors]
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
            raise ValueError(
                "AnnData object is not valid or hasn't been run through validation. "
                "Validate AnnData first before attempting to write labels"
            )

        if validator.adata.isbacked:
            self.adata = validator.adata.to_memory()
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
                    if not value_2 == value_1:
                        raise ValueError(
                            f"Strings types in dependencies cannot be different, {value_1} and {value_2}"
                        )
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
                ids[i], id_suffixes[i] = _curie_remove_suffix(
                    ids[i], curie_constraints["suffixes"]
                )

        for original_id, id, id_suffix in zip(original_ids, ids, id_suffixes):
            # If there are exceptions the label should be the same as the id
            if "exceptions" in curie_constraints:
                if original_id in curie_constraints["exceptions"]:
                    mapping_dict[original_id] = original_id
                    continue

            for ontology_name in allowed_ontologies:
                if ontology_name == "NA":
                    continue
                if ONTOLOGY_CHECKER.is_valid_term_id(ontology_name, id):
                    mapping_dict[original_id] = (
                            ONTOLOGY_CHECKER.get_term_label(ontology_name, id) + id_suffix
                    )

        # Check that all ids got a mapping. All ids should be found if adata was validated
        for id in original_ids:
            if id not in mapping_dict:
                raise ValueError(f"Add labels error: Unable to get label for '{id}'")

        return mapping_dict

    def _get_mapping_dict_features(self, ids: List[str]) -> Dict[str, str]:

        """
        Creates a mapping dictionary of gene/feature IDs and labels.

        :param list[str] ids: Gene/feature IDs use for mapping

        :return a mapping dictionary: {id: label, ...}
        :rtype dict
        """

        mapping_dict = {}

        for i in ids:
            organism = _get_organism_from_feature_id(i)
            print(self.validator.gene_checkers[organism].get_symbol(i))
            mapping_dict[i] = self.validator.gene_checkers[organism].get_symbol(i)

        return mapping_dict

    def _get_labels(
        self, component: str, column: str, column_definition: dict
    ) -> pd.Categorical:

        """
        Retrieves a new column (pandas categorical) with labels based on the IDs in 'column' and the logic in the
        'column_definition'

        :param str component: what dataframe in self.adata to work with
        :param str column: Column in self.adata with IDs that will be used to retrieve values
        :param dict column_definition: schema definition of the column
        e.g. schema_def["obs"]["columns"]["cell_type_ontology_term_id"]

        :rtype pandas.Categorical
        :return new pandas column with labels corresponding to input column
        """

        # Set variables for readability
        type_labels = column_definition["add_labels"]["type"]
        current_df = getattr(self.adata, component)

        if column == "index":
            original_column = pd.Series(current_df.index)
        else:
            original_column = getattr(current_df, column)

        ids = getattr(current_df, column).drop_duplicates().tolist()

        # Flatten column definition (will do so if there are dependencies in the definition
        column_definition = self._flatten_column_def_with_dependencies(
            column_definition
        )

        if type_labels == "curie":

            if "curie_constraints" not in column_definition:
                raise ValueError(
                    f"Schema definition error: 'add_lables' with type 'curie' was found for '{column}' "
                    "but no curie constraints were found for the lables"
                )

            mapping_dict = self._get_mapping_dict_curie(ids, column_definition["curie_constraints"])

        elif type_labels == "feature_id":
            mapping_dict = self._get_mapping_dict_features(ids=ids)

        else:
            raise TypeError(
                f"'{type_labels}' is not supported in 'add-labels' functionality"
            )

        new_column = (
            original_column.copy().replace(mapping_dict).astype("category")
        )

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

        new_column = self._get_labels(component, column, column_definition)

        # The sintax below is a programtic way to access obs and var in adata:
        # adata.__dict__["_obs"] is adata.obs
        self.adata.__dict__["_" + component][column_definition["add_labels"]["to"]] = new_column

    def _add_labels(self):

        """
        From a valid (per cellxgene's schema) adata, this function adds to self.adata ontology/gene labels
        to adata.obs and adata.var respectively
        """

        for component in ["obs", "var"]:

            # Doing it for columns
            for column, column_def in self.schema_def["components"][component]["columns"].items():

                if "add_labels" in column_def:
                    self._add_column(component, column, column_def)

            # Doing it for index
            index_def = self.schema_def["components"][component]["index"]
            if "add_labels" in index_def:
                self._add_column(component, "index", index_def)

    def _check_column_availability(self):

        """
        This method will add error messages to self.errors if reserved columns in self.adata.obs or
        self.adata.var already exist

        :rtype none
        """

        for component in ["obs"]:
            for column, columns_def in self.schema_def["components"]["obs"][
                "columns"
            ].items():
                if "add_labels" in columns_def:
                    reserved_name = columns_def["add_labels"]["to"]
                    if reserved_name in getattr(self.adata, component):
                        self.errors.append(
                            f"Add labels error: Column '{reserved_name}' is a reserved column name "
                            "of 'obs'. Remove it from h5ad and try again."
                        )

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


def validate(h5ad_path: Union[str, bytes, os.PathLike], add_labels_file: str = None):

    """
    Entry point for validation.

    :param Union[str, bytes, os.PathLike] h5ad_path: Path to h5ad file to validate
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
