import re
import sys
import anndata
import os
import pandas as pd
from pandas.core.computation.ops import UndefinedVariableError
from numpy import ndarray
from numpy import count_nonzero
from numpy import cumprod
from numpy import isnan
from numpy import nditer
from scipy import sparse
from typing import List, Dict, Union, Optional, Tuple
from . import ontology
from . import schema
from . import env

ONTOLOGY_CHECKER = ontology.OntologyChecker()


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
        self._raw_layer_exists = None

        # Values will be instances of ontology.GeneChecker,
        # keys will be one of ontology.SupportedOrganisms
        self.gene_checkers = dict()

    @staticmethod
    def getattr_anndata(adata: anndata.AnnData, attr: str = None):

        """
        same as getattr but handles the special case of "raw.var" for an anndata.AndData object

        :param anndata.AnnData adata: the anndata.AnnData object from which to extract an attribute
        :param str attr: name of the attribute to extract

        :return the attribute or none if it does not exist
        """

        if attr == "raw.var":
            if adata.raw:
                return getattr(getattr(adata, "raw"), "var")
            else:
                return None
        else:
            return getattr(adata, attr)

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
        Sets schema dictionary from using information in adata. If there are any errors, it adds them to self.errors
        """

        # Check if schema version slot is present
        if "schema_version" not in self.adata.uns:
            self.errors.append(
                "adata has no schema definition in 'adata.uns'. Validation cannot be performed."
            )
        else:

            # Check if schema version is supported
            version = self.adata.uns["schema_version"]
            path = schema.get_schema_file_path(version)

            if not os.path.isfile(path):
                supported_versions = schema.get_schema_versions_supported()
                self.errors.append(
                    f"Schema version '{version}' is not supported. Current supported versions: '{supported_versions}'. "
                    f"Validation cannot be performed."
                )
            else:
                self.schema_def = schema.get_schema_definition(version)

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
                f"'{term_id}' in '{column_name}' is not an allowed term: '{all_allowed}'."
            )

    def _validate_curie_ancestors(
        self,
        term_id: str,
        column_name: str,
        allowed_ancestors: Dict[str, List[str]],
        inclusive: bool,
    ):

        """
        Validate a single curie term id is a valid children of any of allowed ancestors
        If there are any errors, it adds them to self.errors

        :param str term_id: the curie term id to validate
        :param str column_name: original column name in adata where the term_id comes from (used for error messages)
        :param dict{str: list[str]} allowed_ancestors: keys must be ontology names and values must lists of
        :param bool inclusive:  if True then the ancestors themselves are allowed
        allowed ancestors

        :rtype None
        """

        checks = []

        for ontology_name, ancestors in allowed_ancestors.items():
            for ancestor in ancestors:

                if inclusive and term_id == ancestor:
                    checks.append(True)

                is_valid_term_id = ONTOLOGY_CHECKER.is_valid_term_id(
                    ontology_name, term_id
                )
                is_valid_ancestor_id = ONTOLOGY_CHECKER.is_valid_term_id(
                    ontology_name, ancestor
                )
                if is_valid_term_id & is_valid_ancestor_id:
                    is_child = ONTOLOGY_CHECKER.is_descendent_of(
                        ontology_name, term_id, ancestor
                    )
                    checks.append(is_child)

        if True not in checks:
            all_ancestors = list(allowed_ancestors.values())
            self.errors.append(
                f"'{term_id}' in '{column_name}' is not a child term id of '{all_ancestors}'."
            )

    def _validate_curie_ontology(
        self, term_id: str, column_name: str, allowed_ontologies: List[str]
    ):

        """
        Validate a single curie term id belongs to specified ontologies. If it does belong to an allowed ontology
        verifies that it is not deprecated (obsolete).
        If there are any errors, it adds them to self.errors

        :param str term_id: the curie term id to validate
        :param str column_name: original column name in adata where the term_id comes from (used for error messages)
        :param List[str] allowed_ontologies: allowed ontologies

        :rtype None
        """

        checks = []

        for ontology_name in allowed_ontologies:

            is_valid = ONTOLOGY_CHECKER.is_valid_term_id(ontology_name, term_id)
            checks.append(is_valid)

            if is_valid:
                if ONTOLOGY_CHECKER.is_term_id_deprecated(ontology_name, term_id):
                    self.errors.append(
                        f"'{term_id}' in '{column_name}' is a deprecated term id of '{ontology_name}'."
                    )

        if sum(checks) == 0:
            self.errors.append(
                f"'{term_id}' in '{column_name}' is not a valid ontology term id of '{', '.join(allowed_ontologies)}'."
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

        # If there are forbidden terms
        if "forbidden" in curie_constraints:
            if term_id in curie_constraints["forbidden"]:
                self.errors.append(f"'{term_id}' in '{column_name}' is not allowed'.")
                return

        # If NA is found in allowed ontologies, it means only exceptions should be found. If no exceptions were found
        # then return error
        if curie_constraints["ontologies"] == ["NA"]:
            self.errors.append(
                f"'{term_id}' in '{column_name}' is not a valid value of '{column_name}'."
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
                term_id, column_name, curie_constraints["ancestors"], False
            )

        # Ancestors themselves are allowed
        if "ancestors_inclusive" in curie_constraints:
            self._validate_curie_ancestors(
                term_id, column_name, curie_constraints["ancestors_inclusive"], True
            )

        # If there is a set of allowed terms check for it
        if "allowed_terms" in curie_constraints:
            self._validate_curie_allowed_terms(
                term_id, column_name, curie_constraints["allowed_terms"]
            )

    def _validate_feature_id(self, feature_id: str, df_name: str):

        """
        Validates a feature id, i.e. checks that it's present in the reference
        If there are any errors, it adds them to self.errors and adds it to the list of invalid features

        :param str feature_id: the feature id to be validated
        :param str df_name: name of dataframe the feauter id comes from (var or raw.var)

        :rtype none
        """

        organism = ontology.get_organism_from_feature_id(feature_id)

        if not organism:
            self.errors.append(
                f"Could not infer organism from feature ID '{feature_id}' in '{df_name}', "
                f"make sure it is a valid ID."
            )
            return

        if organism not in self.gene_checkers:
            self.gene_checkers[organism] = ontology.GeneChecker(organism)

        if not self.gene_checkers[organism].is_valid_id(feature_id):
            self.invalid_feature_ids.append(feature_id)
            self.errors.append(
                f"'{feature_id}' is not a valid feature ID in '{df_name}'."
            )

        return

    def _validate_column_feature_is_filtered(
        self, column: pd.Series, column_name: str, df_name: str
    ):

        """
        Validates the "is_feature_filtered" in adata.var. This column must be bool, and for genes that are set to
        True, their expression values in X must be 0.
        If there are any errors, it adds them to self.errors.

        :rtype none
        """

        if not column.dtype == bool:
            self.errors.append(
                f"Column '{column_name}' in dataframe '{df_name}' must be boolean, not '{column.dtype.name}'."
            )

        if sum(column) > 0:

            n_nonzero = 0

            if (
                isinstance(self.adata.X, sparse.csc_matrix)
                or isinstance(self.adata.X, sparse.csr_matrix)
                or isinstance(self.adata.X, sparse.coo_matrix)
            ):

                n_nonzero = self.adata.X[:, column].count_nonzero()

            elif isinstance(self.adata.X, ndarray):
                n_nonzero = count_nonzero(self.adata.X[:, column])

            else:
                self.errors.append(
                    f"X matrix is of type {type(self.adata.X)}, validation of 'feature_is_filtered' "
                    f"cannot be completed."
                )

            if n_nonzero > 0:
                self.errors.append(
                    f"Some features are 'True' in '{column_name}' of dataframe '{df_name}', but there are "
                    f"{n_nonzero} non-zero values in the corresponding columns of the matrix 'X'. All values for "
                    f"these features must be 0."
                )

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

        # error_original_count will count the number of error messages prior to validating the column, this
        # will be useful in case there's an error prefix to be added to errors found here
        error_original_count = len(self.errors)

        if column_def.get("unique"):
            if column.nunique() != len(column):
                self.errors.append(
                    f"Column '{column_name}' in dataframe '{df_name}' is not unique."
                )

        if column_def.get("type") == "bool":
            if not column.dtype == bool:
                self.errors.append(
                    f"Column '{column_name}' in dataframe '{df_name}' must be boolean, not '{column.dtype.name}'."
                )

        if column_def.get("type") == "categorical":
            if not column.dtype.name == "category":
                self.errors.append(
                    f"Column '{column_name}' in dataframe '{df_name}' must be categorical, not {column.dtype.name}."
                )

        if column_def.get("type") == "feature_is_filtered":
            self._validate_column_feature_is_filtered(column, column_name, df_name)

        if "enum" in column_def:
            bad_enums = [
                v for v in column.drop_duplicates() if v not in column_def["enum"]
            ]
            if bad_enums:
                self.errors.append(
                    f"Column '{column_name}' in dataframe '{df_name}' contains invalid values "
                    f"'{bad_enums}'. Values must be one of {column_def['enum']}."
                )

        if column_def.get("type") == "feature_id":

            # Validates each id
            for feature_id in column:
                self._validate_feature_id(feature_id, df_name)

        if column_def.get("type") == "curie":
            if "curie_constraints" not in column_def:
                raise ValueError(
                    f"Corrupt schema definition, no 'curie_constraints' were found for '{column_name}'"
                )
            if "ontologies" not in column_def["curie_constraints"]:
                raise ValueError(
                    f"allowed 'ontologies' must be specified under 'curie constraints' for '{column_name}'"
                )

            for term_id in column.drop_duplicates():
                self._validate_curie(
                    term_id, column_name, column_def["curie_constraints"]
                )

        # Add error suffix to errors found here
        if "error_message_suffix" in column_def:
            error_total_count = len(self.errors)
            for i in range(error_original_count, error_total_count):
                self.errors[i] = (
                    self.errors[i] + " " + column_def["error_message_suffix"]
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
            try:
                column = getattr(
                    df.query(dependency_def["rule"], engine="python"), column_name
                )
            except UndefinedVariableError:
                self.errors.append(
                    f"Checking values with dependencies failed for adata.{df_name}['{column_name}'], "
                    f"this is likely due to missing dependent column in adata.{df_name}."
                )
                return pd.Series()

            all_rules.append(dependency_def["rule"])

            self._validate_column(column, column_name, df_name, dependency_def)

        # Set column with the data that's left
        all_rules = " | ".join(all_rules)
        column = getattr(
            df.query("not (" + all_rules + " )", engine="python"), column_name
        )

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
                        f"Value '{i}' of list '{list_name}' is not a column in 'adata.obs'."
                    )

    def _validate_str_in_dict(self, value, dict_name: str, key: str):

        """
        Validates that a value from a dictionary is a string and it does not have leading, trailing or double spaces.
        Adds errors to self.errors if any

        :param str value: The dictionary to validate
        :param str dict_name: Name of dictionary in the adata (e.g. "uns")
        :param str key: The key in the dictionary

        :rtype bool
        :return True if passed, False otherwise
        """

        is_valid = True

        if not isinstance(value, str):
            self.errors.append(
                f"'{value}' in '{dict_name}['{key}']' is not valid, it must be a string."
            )

            is_valid = False
        else:

            if value != value.rstrip():
                self.errors.append(
                    f"'{value}' in '{dict_name}['{key}']' is not valid, it contains trailing spaces."
                )

                is_valid = False

            if value != value.lstrip():
                self.errors.append(
                    f"'{value}' in '{dict_name}['{key}']' is not valid, it contains leading spaces."
                )

                is_valid = False

            if value.strip() != " ".join(value.split()):
                self.errors.append(
                    f"'{value}' in '{dict_name}['{key}']' is not valid, it contains double spaces."
                )

                is_valid = False

        return is_valid

    def _validate_enum_in_dict(self, value, enum: List[str], dict_name: str, key: str):

        """
        Validates that a value from a dictionary is part of a list. Adds errors to self.errors if any

        :param str value: The dictionary to validate
        :param List[str] enum: The allowed values
        :param str dict_name: Name of dictionary in the adata (e.g. "uns")
        :param str key: The key in the dictionary

        :rtype  None
        """

        if value not in enum:
            self.errors.append(
                f"'{value}' in '{dict_name}['{key}']' is not valid. "
                f"Allowed terms: {enum}."
            )

    def _validate_X_normalization(self, value):

        """
        Verifies the adata.uns["X_normalization"] folllows the schema. Adds errors to self.errors if any

        :param value: value in adata.uns["X_normalization"]

        :rtype None
        """

        if self._validate_str_in_dict(value, "uns", "X_normalization"):

            # If value is "none" then there must NOT be a raw layer
            if value == "none":

                if self._get_raw_x_loc() == "raw.X":
                    self.errors.append(
                        "uns['X_normalization'] is 'none' but 'raw.X' is present. Please indicate "
                        "the normalization used for 'X'."
                    )

                if self._get_raw_x_loc() == "raw.X" and not self._is_raw():
                    self.warnings.append(
                        "uns['X_normalization'] is 'none' but 'raw.X' doesn't appear to have "
                        "raw counts (integers)"
                    )

                if self._get_raw_x_loc() == "X" and not self._is_raw():
                    self.warnings.append(
                        "uns['X_normalization'] is 'none', there is no 'raw.X' and 'X' doesn't appear to have "
                        "raw counts (integers)"
                    )

            else:
                if self._get_raw_x_loc() == "X" and self._is_raw():
                    self.warnings.append(
                        f"uns['X_normalization'] is '{value}' but 'X' seems to have mostly "
                        f"raw values (integers)."
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
                    self.errors.append(f"'{key}' in '{dict_name}' is not present.")
                continue

            value = dictionary[key]

            if value_def["type"] == "string":
                if not self._validate_str_in_dict(value, dict_name, key):
                    continue

                if "enum" in value_def:
                    self._validate_enum_in_dict(
                        value, value_def["enum"], dict_name, key
                    )

            if value_def["type"] == "X_normalization":
                self._validate_X_normalization(value)

            if value_def["type"] == "match_obsm_keys":

                if not self._validate_str_in_dict(value, dict_name, key):
                    continue

                if value not in self.adata.obsm:
                    self.errors.append(
                        f"'{value}' in '{dict_name}['{key}']' is not valid, "
                        f"it must be a key of 'adata.obsm'."
                    )

            if value_def["type"] == "list":
                if not (isinstance(value, list) or isinstance(value, ndarray)):
                    self.errors.append(
                        f"'{value}' in '{dict_name}['{key}']' is not valid, "
                        f"it must be a list or numpy array."
                    )
                    continue

                self._validate_list(key, value, value_def["element_type"])

    def _validate_dataframe(self, df_name: str):

        """
        Verifies the dataframe follows the schema. Adds errors to self.errors if any

        :param str df_name: Name of dataframe in the adata (e.g. "obs")

        :rtype None
        """

        df = self.getattr_anndata(self.adata, df_name)
        df_definition = self._get_component_def(df_name)

        # Validate index if needed
        if "index" in self._get_component_def(df_name):
            self._validate_column(
                pd.Series(df.index),
                "index",
                df_name,
                self._get_column_def(df_name, "index"),
            )

        # Validate columns
        for column_name in df_definition["columns"].keys():
            if column_name not in df.columns:
                self.errors.append(
                    f"Dataframe '{df_name}' is missing column '{column_name}'."
                )
                continue

            column_def = self._get_column_def(df_name, column_name)
            column = getattr(df, column_name)

            # First check if there are dependencies with other columns and work with a subset of the data if so
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

        max_sparsity = float(self.schema_def["sparsity"])

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

                if isinstance(x, sparse.csc_matrix) or isinstance(x, sparse.coo_matrix):
                    sparsity = 1 - x.count_nonzero() / float(cumprod(x.shape)[-1])
                elif isinstance(x, ndarray):
                    sparsity = 1 - count_nonzero(x) / float(cumprod(x.shape)[-1])
                else:
                    self.warnings.append(
                        f"{x_name} matrix is of type {type(x)}, sparsity calculation hasn't been "
                        f"implemented"
                    )
                    continue

                if sparsity > max_sparsity:
                    self.warnings.append(
                        f"Sparsity of '{x_name}' is {sparsity} which is greater than {max_sparsity}, "
                        f"and it is not a 'scipy.sparse.csr_matrix'. It is STRONGLY RECOMMENDED "
                        f"to use this type of matrix for the given sparsity."
                    )

    def _validate_embedding_dict(self):

        """
        Validates the embedding dictionary -- it checks that all values of adata.obms are numpy arrays with the correct
        dimension. Adds errors to self.errors if any. Checks that the keys start with "X_"

        :rtype none
        """

        if not self.adata.obsm:
            self.errors.append("No embeddings found in 'adata.obsm'.")
            return

        obsm_with_x_prefix = 0
        for key, value in self.adata.obsm.items():

            if not isinstance(value, ndarray):
                self.errors.append(
                    f"All embeddings have to be of 'numpy.ndarray' type, "
                    f"'adata.obsm['{key}']' is {type(value)}')."
                )
                continue

            # Embeddings to be shown in cellxgene explorer
            if key.startswith("X_"):
                obsm_with_x_prefix += 1

                if len(value.shape) < 2:
                    self.errors.append(
                        f"All embeddings must have as many rows as cells, and at least two columns."
                        f"'adata.obsm['{key}']' has shape of '{value.shape}'."
                    )
                else:
                    if value.shape[0] != self.adata.n_obs or value.shape[1] < 2:
                        self.errors.append(
                            f"All embeddings must have as many rows as cells, and at least two columns."
                            f"'adata.obsm['{key}']' has shape of '{value.shape}'."
                        )

        if obsm_with_x_prefix == 0:
            self.errors.append(
                "At least one embedding in 'obsm' has to have a key with an 'X_' prefix."
            )

    def _are_children_of(
        self, component: str, column: str, ontology_name: str, ancestors: List[str]
    ) -> bool:

        """
        Checks if elements in the specified column of the component (e.g. 'assay_ontology_term_id' of 'adata.obs') are
        children of the given ancestors.

        Ancestors checks are inclusive, meaning that a value is its own ancestor as well.

        :param str component: the name of the component that's been checked.
        :param str column: Column in the component to check
        :param str ontology_name: Name of the ontology (e.g. "EFO")
        :param List[str] ancestors: List of ancestors

        :rtype bool
        :return True if any value in column is children of any ancestor.
        """

        curies = getattr(getattr(self.adata, component), column)
        curies = curies.drop_duplicates()

        for curie in curies:
            if ONTOLOGY_CHECKER.is_valid_term_id(ontology_name, curie):
                curie_ancestors = ONTOLOGY_CHECKER.get_term_ancestors(
                    ontology_name, curie
                )
                curie_ancestors.add(curie)
                if bool(curie_ancestors & set(ancestors)):
                    return True

        return False

    def _get_raw_x(self) -> Union[ndarray, sparse.csc_matrix, sparse.csr_matrix]:

        """
        gets raw x (best guess, i.e. not guarantee it's actually raw)
        """

        if self.adata.raw:
            return self.adata.raw.X
        else:
            return self.adata.X

    def _get_raw_x_loc(self) -> str:

        """
        gets raw x location (best guess, i.e. not guarantee it's actually raw)
        """

        if self.adata.raw:
            return "raw.X"
        else:
            return "X"

    def _is_raw(self, max_values_to_check: int = 5000, force: bool = False) -> bool:

        """
        Checks if the first non-zero "max_values_to_check" in the best guess for the raw matrix (adata.X or adata.raw.X)
        are integers. Returns False if at least one value is not an integer,

        True otherwise.

        Since this process is memory intensive, it will return a cache value if this function has been called before.
        If calculation needs to be repeated use `force = True`

        :param int max_values_to_check: total values to check, default set to 5000 due to performance concerns.

        :rtype bool
        :return False if at least one value is not an integer, True otherwise
        """

        if force:
            self._raw_layer_exists = None

        if self._raw_layer_exists is None:

            # Get potential raw_X
            raw_loc = self._get_raw_x_loc()
            if raw_loc == "raw.X":
                x = self.adata.raw.X
            else:
                x = self.adata.X

            # Get array without zeros
            non_zeroes_index = x.nonzero()

            if max_values_to_check > len(non_zeroes_index[0]):
                max_values_to_check = len(non_zeroes_index[0])

            x_non_zeroes = x[non_zeroes_index[0][:max_values_to_check], non_zeroes_index[1][:max_values_to_check]]

            # If all values are zeros then is raw, otherwise if a single value is not an int then return is not raw
            if x_non_zeroes.size < 1:
                self._raw_layer_exists = True
            else:
                self._raw_layer_exists = True
                for i in nditer(x_non_zeroes):

                    if isnan(i):
                        continue

                    if i % 1 != 0:
                        self._raw_layer_exists = False
                        break

        return self._raw_layer_exists

    def _validate_x_raw_x_dimensions(self):

        """
        Validates that X and raw.X have the same shape. Adds errors to self.errors if any.
        """

        if self._get_raw_x_loc() == "raw.X":
            if self.adata.n_vars != self.adata.raw.n_vars:
                self.errors.append(
                    f"Number of genes in X ({self.adata.n_vars}) is different "
                    f"than raw.X ({self.adata.raw.n_vars})."
                )
            else:
                if not (self.adata.var.index == self.adata.raw.var.index).all():
                    self.errors.append(
                        "Index of 'raw.var' is not identical to index of 'var'."
                    )
            if self.adata.n_obs != self.adata.raw.n_obs:
                self.errors.append(
                    f"Number of cells in X ({self.adata.n_obs}) is different "
                    f"than raw.X ({self.adata.raw.n_obs})."
                )
            else:
                if not (self.adata.obs_names == self.adata.raw.obs_names).all():
                    self.errors.append("Cells in X and raw.X are different.")

    def _validate_raw(self):

        """
        Validates raw only if the rules in the schema definition are fulfilled and that X and raw.X have the same shape
        The validation entails checking that:
         1. X and raw.X have the same column and and raw indeces
         2. there's an expression matrix containing raw (integer) values, first in adata.raw.X and then adata.X if
         the former does not exist.

        Adds errors to self.errors if any.

        :rtype None
        """

        # Check that raw and raw.X have the same shape
        self._validate_x_raw_x_dimensions()

        # Asses if we actually need to perform validation of raw based on the rules in the schema
        # As of now, this means that we only do validation of raw if it's RNA data
        checks = []
        for component, component_rules in self.schema_def["raw"].items():
            for column, column_rules in component_rules.items():
                for rule, rule_def in column_rules.items():
                    if rule == "not_children_of":
                        for ontology_name, ancestors in rule_def.items():
                            checks.append(
                                not self._are_children_of(
                                    component, column, ontology_name, ancestors
                                )
                            )
                    else:
                        raise ValueError(
                            f"'{rule}' rule in raw definition of the schema is not implemented "
                        )

        # If all checks passed then proceed with validation
        if all(checks):

            normalization = self.adata.uns["X_normalization"]

            # If raw data is missing: no "raw.X", and "X_normalization" is not "none"
            if (
                not self._is_raw()
                and self._get_raw_x_loc() == "X"
                and normalization != "none"
            ):

                self.errors.append(
                    "Raw data is missing: there is no 'raw.X' and 'X_normalization' is not 'none'."
                )

            # If raw data is in X but X_normalization is NOT none raise an error
            if (
                self._is_raw()
                and self._get_raw_x_loc() == "X"
                and normalization != "none"
            ):

                self.errors.append(
                    f"uns['X_normalization'] is '{normalization}' but raw data seems to be in X, "
                    f"if X is raw then uns['X_normalization'] MUST be 'none'."
                )

            # If raw data is in X and there is nothing in raw.X (i.e. normalized values are not provided), then
            # add a warning because normalized data for RNA data is STRONGLY RECOMMENDED
            if (
                self._is_raw()
                and self._get_raw_x_loc() == "X"
                and normalization == "none"
            ):

                self.warnings.append(
                    "Only raw data was found, i.e. there is no 'raw.X' and 'uns['X_normalization']' is 'none'. "
                    "It is STRONGLY RECOMMENDED that 'final' (normalized) data is provided."
                )

    def _check_single_column_availability(
        self, component: str, add_labels_def: List[dict]
    ):

        """
        This method checks a single reserved column in adata.obs or adata.var and adds a message to self.error if
        it already exists

        :param str component: the name of the component that's been checked.
        :param List[dict] add_labels_def: the "add_labels" definition, contains the information of
        the reserved column(s)

        :rtype none
        """

        for label_def in add_labels_def:
            reserved_name = label_def["to_column"]

            if reserved_name in self.getattr_anndata(self.adata, component):
                self.errors.append(
                    f"Add labels error: Column '{reserved_name}' is a reserved column name "
                    f"of '{component}'. Remove it from h5ad and try again."
                )

    def _check_column_availability(self):

        """
        This method will check for columns that are reserved in self.adata.obs or
        self.adata.var already exist

        :rtype none
        """

        for component, component_def in self.schema_def["components"].items():

            # If not a dataframe we don't need to check for columns
            if component_def["type"] != "dataframe":
                continue

            # Skip if component does not exist
            if self.getattr_anndata(self.adata, component) is None:
                continue

            # Do it for columns that are forbidden
            if "forbidden_columns" in component_def:
                for column in component_def["forbidden_columns"]:
                    if column in self.getattr_anndata(self.adata, component):
                        self.errors.append(
                            f"Column '{column}' must not be present in '{component}'."
                        )

            # Do it for columns that map to columns
            for column, columns_def in component_def["columns"].items():

                if "add_labels" in columns_def:
                    self._check_single_column_availability(
                        component, columns_def["add_labels"]
                    )

            # Do it for index that map to columns
            if "index" in component_def:

                index_def = component_def["index"]
                if "add_labels" in index_def:
                    self._check_single_column_availability(
                        component, index_def["add_labels"]
                    )

    def _deep_check(self):

        """
        Perform a "deep" check of the AnnData object using the schema definition. Adds errors to self.errors if any

        :rtype None
        """

        # Checks that adata is fully loaded
        if self.adata.isbacked:
            raise RuntimeError(
                "adata is loaded in backed mode, validation requires anndata to be loaded in memory"
            )

        # Checks that reserved columns are not used
        self._check_column_availability()

        # Checks sparsity
        self._validate_sparsity()

        # Checks each component
        for component, component_def in self.schema_def["components"].items():

            # Skip if component does not exist: only useful for adata.raw.var
            if self.getattr_anndata(self.adata, component) is None:
                if "required" in component_def:
                    self.errors.append(
                        f"'{component}' is missing from adata and it's required."
                    )
                continue

            if component_def["type"] == "dataframe":
                self._validate_dataframe(component)
            elif component_def["type"] == "dict":
                dictionary = getattr(self.adata, component)
                self._validate_dict(dictionary, component, component_def)
            elif component_def["type"] == "embedding_dict":
                self._validate_embedding_dict()
            else:
                raise ValueError(f"Unexpected component type '{component['type']}'")

        # Checks for raw only if there are no errors, because it depends on the
        # existence of adata.obs["assay_ontology_term_id"] and adata.obs["X_normalization"]
        if not self.errors and "raw" in self.schema_def:
            self._validate_raw()
        else:
            self.warnings.append(
                "Validation of raw layer was not performed due to current errors, try again after "
                "fixing current errors."
            )

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

        if not self.errors:
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
                    if key == "error_message_suffix":
                        merged_dict[key] = value_1 + " " + value_2
                        continue
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

    def _get_mapping_dict_curie(
        self, ids: List[str], curie_constraints: dict
    ) -> Dict[str, str]:

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

    def _get_mapping_dict_feature_reference(
        self, ids: List[str]
    ) -> Dict[str, Optional[ontology.SupportedOrganisms]]:

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
        current_df = Validator.getattr_anndata(self.adata, component)

        if column == "index":
            original_column = pd.Series(current_df.index)
            original_column.index = current_df.index
        else:
            original_column = getattr(current_df, column)

        ids = getattr(current_df, column).drop_duplicates().tolist()

        # Flatten column definition (will do so if there are dependencies in the definition
        column_definition = self._flatten_column_def_with_dependencies(
            column_definition
        )

        if label_type == "curie":

            if "curie_constraints" not in column_definition:
                raise ValueError(
                    f"Schema definition error: 'add_lables' with type 'curie' was found for '{column}' "
                    "but no curie constraints were found for the lables"
                )

            mapping_dict = self._get_mapping_dict_curie(
                ids, column_definition["curie_constraints"]
            )

        elif label_type == "feature_id":
            mapping_dict = self._get_mapping_dict_feature_id(ids=ids)

        elif label_type == "feature_reference":
            mapping_dict = self._get_mapping_dict_feature_reference(ids=ids)

        else:
            raise TypeError(
                f"'{label_type}' is not supported in 'add-labels' functionality"
            )

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

            new_column = self._get_labels(
                component, column, column_definition, label_def["type"]
            )
            new_column_name = label_def["to_column"]

            # The sintax below is a programtic way to access obs and var in adata:
            # adata.__dict__["_obs"] is adata.obs
            self.adata.__dict__["_" + component][new_column_name] = new_column

    def _add_labels(self):

        """
        From a valid (per cellxgene's schema) adata, this function adds to self.adata ontology/gene labels
        to adata.obs and adata.var respectively
        """

        for component in ["obs", "var"]:

            # Doing it for columns
            for column, column_def in self.schema_def["components"][component][
                "columns"
            ].items():

                if "add_labels" in column_def:
                    self._add_column(component, column, column_def)

            # Doing it for index
            index_def = self.schema_def["components"][component]["index"]
            if "add_labels" in index_def:
                self._add_column(component, "index", index_def)

    def write_labels(self, add_labels_file: str):

        """
        From a valid (per cellxgene's schema) h5ad, this function writes a new h5ad file with ontology/gene labels added
        to adata.obs  and adata.var respectively

        :param str add_labels_file: Path to new h5ad file with ontology/gene labels added

        :rtype None
        """

        # Add labels in obs
        self._add_labels()

        # Write file
        try:
            self.adata.write_h5ad(add_labels_file, compression="gzip")
        except Exception as e:
            self.errors.append(f"Writing h5ad was unsuccessful, got exception '{e}'.")

        # Print errors if any
        if self.errors:
            self.errors = ["ERROR: " + i for i in self.errors]
            print(*self.errors, sep="\n")
            self.was_writing_successful = False
        else:
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

    if add_labels_file:
        writer = AnnDataLabelAppender(validator)
        writer.write_labels(add_labels_file)

        return validator.is_valid & writer.was_writing_successful

    return True
