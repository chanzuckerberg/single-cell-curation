import logging
import math
import os
import re
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Union

import anndata
import numpy as np
import pandas as pd
import semver
from pandas.core.computation.ops import UndefinedVariableError
from scipy import sparse

from . import env, ontology, schema
from .utils import get_matrix_format, getattr_anndata, read_h5ad

logger = logging.getLogger(__name__)

ONTOLOGY_CHECKER = ontology.OntologyChecker()


class Validator:
    """Handles validation of AnnData"""

    schema_definitions_dir = env.SCHEMA_DEFINITIONS_DIR

    def __init__(self, ignore_labels=False):
        # Set initial state
        self.errors = []
        self.warnings = []
        self.is_valid = False
        self.adata = anndata.AnnData()
        self.schema_def = dict()
        self.schema_version: str = None
        self.h5ad_path = ""
        self.invalid_feature_ids = []
        self._raw_layer_exists = None
        self.is_seurat_convertible: bool = True
        self.ignore_labels = ignore_labels

        # Values will be instances of ontology.GeneChecker,
        # keys will be one of ontology.SupportedOrganisms
        self.gene_checkers = dict()

        # Matrix (e.g., X, raw.X, ...) number non-zero cache
        self.number_non_zero = dict()

    @staticmethod
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

    def _validate_encoding_version(self):
        import h5py

        with h5py.File(self.h5ad_path, "r") as f:
            encoding_dict = dict(f.attrs)
            encoding_version = encoding_dict.get("encoding-version")
            if encoding_version != "0.1.0":
                self.errors.append("The h5ad artifact was generated with an AnnData version different from 0.8.0.")

    def _set_schema_def(self):
        """
        Sets schema dictionary from using information in adata. If there are any errors, it adds them to self.errors
        """
        version = self.adata.uns.get("schema_version")
        supported_major_versions = schema.get_schema_versions_supported()
        supported_major_versions.sort()
        latest_major_version: semver.Version = semver.Version.parse(supported_major_versions[-1])
        if version and not semver.Version.parse(version).is_compatible(latest_major_version):
            logger.warning(
                f"Schema version '{version}' is not supported. Current supported versions: '{latest_major_version}'. "
                f"Validating with latest version '{latest_major_version}'."
            )
        if not self.schema_version:
            self.schema_version = schema.get_current_schema_version()
            self.schema_def = schema.get_schema_definition(str(latest_major_version))

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
        for ontology_name, allowed_terms in terms.items():
            if ONTOLOGY_CHECKER.is_valid_term_id(ontology_name, term_id):
                checks.append(term_id in allowed_terms)

        if sum(checks) == 0 and len(checks) > 0:
            all_allowed = list(terms.values())
            self.errors.append(f"'{term_id}' in '{column_name}' is not an allowed term: '{all_allowed}'.")

    def _validate_curie_ancestors(
        self,
        term_id: str,
        column_name: str,
        allowed_ancestors: Dict[str, List[str]],
        inclusive: bool,
        errors=True,
    ):
        """
        Validate a single curie term id is a valid children of any of allowed ancestors
        If there are any errors, it adds them to self.errors if the errors flag is True.

        :param str term_id: the curie term id to validate
        :param str column_name: original column name in adata where the term_id comes from (used for error messages)
        :param dict{str: list[str]} allowed_ancestors: keys must be ontology names and values must lists of
        allowed ancestors
        :param bool inclusive:  if True then the ancestors themselves are allowed
        :param bool errors: if True then errors are appended to self.errors

        :rtype Bool
        """

        checks = []

        for ontology_name, ancestors in allowed_ancestors.items():
            for ancestor in ancestors:
                if inclusive and term_id == ancestor:
                    checks.append(True)

                is_valid_term_id = ONTOLOGY_CHECKER.is_valid_term_id(ontology_name, term_id)
                is_valid_ancestor_id = ONTOLOGY_CHECKER.is_valid_term_id(ontology_name, ancestor)
                if is_valid_term_id & is_valid_ancestor_id:
                    is_child = ONTOLOGY_CHECKER.is_descendent_of(ontology_name, term_id, ancestor)
                    checks.append(is_child)

        if True not in checks:
            if errors:
                all_ancestors = list(allowed_ancestors.values())
                self.errors.append(f"'{term_id}' in '{column_name}' is not a child term id of '{all_ancestors}'.")
            return False
        return True

    def _validate_curie_ontology(self, term_id: str, column_name: str, allowed_ontologies: List[str]):
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

            if is_valid and ONTOLOGY_CHECKER.is_term_id_deprecated(ontology_name, term_id):
                self.errors.append(f"'{term_id}' in '{column_name}' is a deprecated term id of '{ontology_name}'.")

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
        if "exceptions" in curie_constraints and term_id in curie_constraints["exceptions"]:
            return

        # If there are forbidden terms
        if "forbidden" in curie_constraints and term_id in curie_constraints["forbidden"]:
            self.errors.append(f"'{term_id}' in '{column_name}' is not allowed'.")
            return

        # If NA is found in allowed ontologies, it means only exceptions should be found. If no exceptions were found
        # then return error
        if curie_constraints["ontologies"] == ["NA"]:
            self.errors.append(f"'{term_id}' in '{column_name}' is not a valid value of '{column_name}'.")
            return

        # Check if there are any allowed suffixes and remove them if needed
        if "suffixes" in curie_constraints:
            term_id, suffix = self._curie_remove_suffix(term_id, curie_constraints["suffixes"])

        # Check that term id belongs to allowed ontologies
        self._validate_curie_ontology(term_id, column_name, curie_constraints["ontologies"])

        # If the term id does not belong to an allowed ontology, the subsequent checks are redundant
        if self.errors:
            return

        # If there are specified ancestors then make sure that this id is a valid child
        if "ancestors" in curie_constraints:
            self._validate_curie_ancestors(term_id, column_name, curie_constraints["ancestors"], False)

        # Ancestors themselves are allowed
        if "ancestors_inclusive" in curie_constraints:
            self._validate_curie_ancestors(term_id, column_name, curie_constraints["ancestors_inclusive"], True)

        # If there is a set of allowed terms check for it
        if "allowed_terms" in curie_constraints:
            self._validate_curie_allowed_terms(term_id, column_name, curie_constraints["allowed_terms"])

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
            self.errors.append(f"'{feature_id}' is not a valid feature ID in '{df_name}'.")

        return

    def _chunk_matrix(
        self,
        matrix: Union[np.ndarray, sparse.spmatrix],
        obs_chunk_size: Optional[int] = 10_000,
    ):
        """
        Iterator which chunks the _named_ or _specified_ matrix by the
        first (obs) dimension

        The parameter type restrictions are strictly for ensuring that the
        AnnData read fast-path is used (as of AnnData 0.8.0).

        Iterator produces a sequence of tuples, each containing
        (chunk, start, end)
        """
        start = 0
        n = matrix.shape[0]
        for i in range(int(n // obs_chunk_size)):
            logger.debug(f"_chunk_matrix [{i} of {math.ceil(n/obs_chunk_size)}]")
            end = start + obs_chunk_size
            yield (matrix[start:end], start, end)
            start = end
        if start < n:
            yield (matrix[start:n], start, n)

    def _count_matrix_nonzero(self, matrix_name: str, matrix: Union[np.ndarray, sparse.spmatrix]) -> int:
        if matrix_name in self.number_non_zero:
            return self.number_non_zero[matrix_name]

        logger.debug(f"Counting non-zero values in {matrix_name}")

        nnz = 0
        format = get_matrix_format(self.adata, matrix)
        for matrix_chunk, _, _ in self._chunk_matrix(matrix):
            nnz += matrix_chunk.count_nonzero() if format != "dense" else np.count_nonzero(matrix_chunk)

        self.number_non_zero[matrix_name] = nnz
        return nnz

    def _validate_column_feature_is_filtered(self, column: pd.Series, column_name: str, df_name: str):
        """
        Validates the "is_feature_filtered" in adata.var. This column must be bool, and for genes that are set to
        True, their expression values in X must be 0.
        If there are any errors, it adds them to self.errors.

        :rtype none
        """

        if column.dtype != bool:
            self.errors.append(
                f"Column '{column_name}' in dataframe '{df_name}' must be boolean, not '{column.dtype.name}'."
            )
            return

        if sum(column) > 0:
            n_nonzero = 0

            X_format = get_matrix_format(self.adata, self.adata.X)
            if X_format in ["csc", "csr", "coo"]:
                n_nonzero = self.adata.X[:, column].count_nonzero()

            elif X_format == "dense":
                n_nonzero = np.count_nonzero(self.adata.X[:, column])

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

    def _validate_column(self, column: pd.Series, column_name: str, df_name: str, column_def: dict):
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

        if column_def.get("unique") and column.nunique() != len(column):
            self.errors.append(f"Column '{column_name}' in dataframe '{df_name}' is not unique.")

        if column_def.get("type") == "bool" and column.dtype != bool:
            self.errors.append(
                f"Column '{column_name}' in dataframe '{df_name}' must be boolean, not '{column.dtype.name}'."
            )

        if column_def.get("type") == "categorical":
            if column.dtype.name != "category":
                self.errors.append(
                    f"Column '{column_name}' in dataframe '{df_name}' must be categorical, not {column.dtype.name}."
                )
            else:
                if column_def.get("subtype") == "str":
                    if column.dtype.categories.dtype != "object" and column.dtype.categories.dtype != "string":
                        self.errors.append(
                            f"Column '{column_name}' in dataframe '{df_name}' must be object or string, not"
                            f" {column.dtype.categories.dtype}."
                        )
                    else:
                        if any(len(cat.strip()) == 0 for cat in column.dtype.categories):
                            self.errors.append(
                                f"Column '{column_name}' in dataframe '{df_name}' must not contain empty values."
                            )

                # check for null values--skip on column defs with enums, since it will already be part of that check
                if not column_def.get("enum") and column.isnull().any():
                    self.errors.append(f"Column '{column_name}' in dataframe '{df_name}' must not contain NaN values.")

        if column_def.get("type") == "feature_is_filtered":
            self._validate_column_feature_is_filtered(column, column_name, df_name)

        if "enum" in column_def:
            bad_enums = [v for v in column.drop_duplicates() if v not in column_def["enum"]]
            if bad_enums:
                self.errors.append(
                    f"Column '{column_name}' in dataframe '{df_name}' contains invalid values "
                    f"'{bad_enums}'. Values must be one of {column_def['enum']}"
                )

        if column_def.get("type") == "feature_id":
            # Validates each id
            for feature_id in column:
                self._validate_feature_id(feature_id, df_name)

        if column_def.get("type") == "curie":
            # Check for NaN values
            if column.isnull().any():
                self.errors.append(f"Column '{column_name}' in dataframe '{df_name}' must not contain NaN values.")
                return

            if "curie_constraints" not in column_def:
                raise ValueError(f"Corrupt schema definition, no 'curie_constraints' were found for '{column_name}'")
            if "ontologies" not in column_def["curie_constraints"]:
                raise ValueError(
                    f"allowed 'ontologies' must be specified under 'curie constraints' for '{column_name}'"
                )

            for term_id in column.drop_duplicates():
                self._validate_curie(term_id, column_name, column_def["curie_constraints"])

        # Add error suffix to errors found here
        if "error_message_suffix" in column_def:
            error_total_count = len(self.errors)
            for i in range(error_original_count, error_total_count):
                self.errors[i] = self.errors[i] + " " + column_def["error_message_suffix"]

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
            if "complex_rule" in dependency_def:
                if "match_ancestors" in dependency_def["complex_rule"]:
                    query_fn, args = self._generate_match_ancestors_query_fn(
                        dependency_def["complex_rule"]["match_ancestors"]
                    )
                    term_id, ontologies, ancestors, ancestor_inclusive = args
                    query_exp = f"@query_fn({term_id}, {ontologies}, {ancestors}, {ancestor_inclusive})"
            elif "rule" in dependency_def:
                query_exp = dependency_def["rule"]
            else:
                continue

            try:
                column = getattr(df.query(query_exp, engine="python"), column_name)
            except UndefinedVariableError:
                self.errors.append(
                    f"Checking values with dependencies failed for adata.{df_name}['{column_name}'], "
                    f"this is likely due to missing dependent column in adata.{df_name}."
                )
                return pd.Series(dtype=np.float64)

            all_rules.append(query_exp)

            self._validate_column(column, column_name, df_name, dependency_def)

        # Set column with the data that's left
        all_rules = " | ".join(all_rules)
        column = getattr(df.query("not (" + all_rules + " )", engine="python"), column_name)

        return column

    def _generate_match_ancestors_query_fn(self, rule_def: Dict):
        """
        Generates vectorized function and args to query a pandas dataframe. Function will determine whether values from
        a specified column is a child term to a group of specified ancestors, returning a Bool.
        :param rule_def: defines arguments to pass into vectorized ancestor match validation function
        :return: Tuple(function, Tuple(str, List[str], List[str]))
        """
        validate_curie_ancestors_vectorized = np.vectorize(self._validate_curie_ancestors)
        ancestor_map = rule_def["ancestors"]
        inclusive = rule_def["inclusive"]

        # hack: pandas dataframe query doesn't support Dict inputs
        ontology_keys = []
        ancestor_list = []
        for key, val in ancestor_map.items():
            ontology_keys.append(key)
            ancestor_list.append(val)

        def is_ancestor_match(
            term_id: str,
            ontologies: List[str],
            ancestors: List[str],
            ancestor_inclusive: bool,
        ) -> bool:
            allowed_ancestors = dict(zip(ontologies, ancestors))
            return validate_curie_ancestors_vectorized(term_id, "", allowed_ancestors, ancestor_inclusive, False)

        return is_ancestor_match, (
            rule_def["column"],
            ontology_keys,
            ancestor_list,
            inclusive,
        )

    def _validate_list(self, list_name: str, current_list: List[str], element_type: str):
        """
        Validates the elements of a list based on the type definition. Adds errors to self.errors if any

        :param str list_name: name of list to use for error messages (if any)
        :param str current_list: the list to be validated
        :param str element_type: type to be validated

        :rtype None
        """

        for i in current_list:
            if element_type == "match_obs_columns" and i not in self.adata.obs.columns:
                self.errors.append(f"Value '{i}' of list '{list_name}' is not a column in 'adata.obs'.")

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
            self.errors.append(f"'{value}' in '{dict_name}['{key}']' is not valid, it must be a string.")

            is_valid = False
        else:
            if value != value.rstrip():
                self.errors.append(f"'{value}' in '{dict_name}['{key}']' is not valid, it contains trailing spaces.")

                is_valid = False

            if value != value.lstrip():
                self.errors.append(f"'{value}' in '{dict_name}['{key}']' is not valid, it contains leading spaces.")

                is_valid = False

            if value.strip() != " ".join(value.split()):
                self.errors.append(f"'{value}' in '{dict_name}['{key}']' is not valid, it contains double spaces.")

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
            self.errors.append(f"'{value}' in '{dict_name}['{key}']' is not valid. " f"Allowed terms: {enum}.")

    def _validate_dict(self, dictionary: dict, dict_name: str, dict_def: dict):
        """
        Verifies the dictionary follows the schema. Adds errors to self.errors if any

        :param str dictionary: The dictionary to validate
        :param str dict_name: Name of dictionary in the adata (e.g. "uns")
        :param str dict_def: The schema definition for this specific dictionary

        :rtype None
        """

        for key, value_def in dict_def["keys"].items():
            logger.debug(f"Validating uns dict for key: {key}")
            if key not in dictionary:
                if "required" in value_def:
                    self.errors.append(f"'{key}' in '{dict_name}' is not present.")
                continue

            value = dictionary[key]

            if value_def["type"] == "string":
                if not self._validate_str_in_dict(value, dict_name, key):
                    continue

                if "enum" in value_def:
                    self._validate_enum_in_dict(value, value_def["enum"], dict_name, key)

            if value_def["type"] == "match_obsm_keys":
                if not self._validate_str_in_dict(value, dict_name, key):
                    continue

                if value not in self.adata.obsm:
                    self.errors.append(
                        f"'{value}' in '{dict_name}['{key}']' is not valid, " f"it must be a key of 'adata.obsm'."
                    )

            if value_def["type"] == "list":
                if not (isinstance(value, (list, np.ndarray))):
                    self.errors.append(
                        f"'{value}' in '{dict_name}['{key}']' is not valid, " f"it must be a list or numpy array."
                    )
                    continue

                self._validate_list(key, value, value_def["element_type"])

    def _validate_dataframe(self, df_name: str):
        """
        Verifies the dataframe follows the schema. Adds errors to self.errors if any

        :param str df_name: Name of dataframe in the adata (e.g. "obs")

        :rtype None
        """

        df = getattr_anndata(self.adata, df_name)
        df_definition = self._get_component_def(df_name)

        # Validate index if needed
        if "index" in self._get_component_def(df_name):
            logger.debug("Validating index...")
            self._validate_column(
                pd.Series(df.index),
                "index",
                df_name,
                self._get_column_def(df_name, "index"),
            )

        low_rows_threshold = self._get_component_def(df_name).get("warn_if_less_than_rows")
        if low_rows_threshold is not None:
            num_rows = df.shape[0]
            if num_rows < low_rows_threshold:
                self.warnings.append(
                    f"Dataframe '{df_name}' only has {num_rows} rows. "
                    f"Features SHOULD NOT be filtered from expression matrix."
                )

        # Check for columns that have a category defined 0 times (obs only)
        if df_name == "obs":
            for column in df.columns:
                col = df[column]
                if col.dtype != "category":
                    continue
                for category in col.dtype.categories:
                    if category not in col.values:
                        self.warnings.append(
                            f"Column '{column}' in dataframe '{df_name}' contains a category '{category}' with "
                            f"zero observations. These categories will be removed when `--add-labels` flag is present."
                        )

        # Validate columns
        if "columns" in df_definition:
            for column_name in df_definition["columns"]:
                logger.debug(f"Validating column: {column_name}...")
                if column_name not in df.columns:
                    self.errors.append(f"Dataframe '{df_name}' is missing column '{column_name}'.")
                    continue

                column_def = self._get_column_def(df_name, column_name)
                column = getattr(df, column_name)

                # First check if there are dependencies with other columns and work with a subset of the data if so
                if "dependencies" in column_def:
                    column = self._validate_column_dependencies(df, df_name, column_name, column_def["dependencies"])

                # If after validating dependencies there's still values in the column, validate them.
                if len(column) > 0:
                    if "warning_message" in column_def:
                        self.warnings.append(column_def["warning_message"])
                    self._validate_column(column, column_name, df_name, column_def)

    def _validate_sparsity(self):
        """
        calculates sparsity of x and raw.x, if bigger than indicated in the schema and not a scipy sparse matrix, then
        adds to warnings

        :rtype none
        """
        max_sparsity = float(self.schema_def["sparsity"])

        to_validate = [(self.adata.X, "X")]

        # check if there's raw data
        if self.adata.raw:
            to_validate.append((self.adata.raw.X, "raw.X"))

        # check if there's other expression matrices under layers
        if self.adata.layers:
            for key, value in self.adata.layers.items():
                to_validate.append((value, f"layers['{key}']"))

        # Check sparsity
        for x, x_name in to_validate:
            matrix_format = get_matrix_format(self.adata, x)
            if matrix_format == "csr":
                continue
            assert format != "unknown"

            # It seems silly to perform this test for 'coo' and 'csc' formats,
            # which are, by definition, already sparse. But the old code
            # performs test, and so we continue the tradition. It is possible
            # that the prolog comment is incorrect, and the purpose of this
            # function is to recommend CSR for _any_ matrix with sparsity beyond
            # a given limit.

            nnz = self._count_matrix_nonzero(x_name, x)
            sparsity = 1 - nnz / np.prod(x.shape)
            if sparsity > max_sparsity:
                self.warnings.append(
                    f"Sparsity of '{x_name}' is {sparsity} which is greater than {max_sparsity}, "
                    f"and it is not a 'scipy.sparse.csr_matrix'. It is STRONGLY RECOMMENDED "
                    f"to use this type of matrix for the given sparsity."
                )

    def _validate_seurat_convertibility(self):
        """
        Use length of component matrices to determine if the anndata object will be unable to be converted to Seurat by
        virtue of the R language's array size limit (4-byte signed int length). Add warning for each matrix which is
        too large.
        rtype: None
        """
        to_validate = [(self.adata.X, "X")]
        # check if there's raw data
        if self.adata.raw:
            to_validate.append((self.adata.raw.X, "raw.X"))
        # Check length of component arrays
        for matrix, matrix_name in to_validate:
            format = get_matrix_format(self.adata, matrix)
            if format in ["csc", "csr", "coo"]:
                effective_r_array_size = self._count_matrix_nonzero(matrix_name, matrix)
                is_sparse = True
            elif format == "dense":
                effective_r_array_size = max(matrix.shape)
                is_sparse = False
            else:
                self.warnings.append(
                    f"Unable to verify seurat convertibility for matrix {matrix_name} " f"of type {type(matrix)}"
                )
                continue

            if effective_r_array_size > self.schema_def["max_size_for_seurat"]:
                if is_sparse:
                    self.warnings.append(
                        f"This dataset cannot be converted to the .rds (Seurat v4) format. "
                        f"{effective_r_array_size} nonzero elements in matrix {matrix_name} exceed the "
                        f"limitations in the R dgCMatrix sparse matrix class (2^31 - 1 nonzero "
                        f"elements)."
                    )
                else:
                    self.warnings.append(
                        f"This dataset cannot be converted to the .rds (Seurat v4) format. "
                        f"{effective_r_array_size} elements in at least one dimension of matrix "
                        f"{matrix_name} exceed the limitations in the R dgCMatrix sparse matrix class "
                        f"(2^31 - 1 nonzero elements)."
                    )

                self.is_seurat_convertible = False

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
            if not isinstance(value, np.ndarray):
                self.errors.append(
                    f"All embeddings have to be of 'numpy.ndarray' type, " f"'adata.obsm['{key}']' is {type(value)}')."
                )
                continue

            # Embeddings to be shown in cellxgene explorer
            if key.startswith("X_"):
                obsm_with_x_prefix += 1

                if len(value.shape) < 2 or value.shape[0] != self.adata.n_obs or value.shape[1] < 2:
                    self.errors.append(
                        f"All embeddings must have as many rows as cells, and at least two columns."
                        f"'adata.obsm['{key}']' has shape of '{value.shape}'."
                    )

        if obsm_with_x_prefix == 0:
            self.errors.append("At least one embedding in 'obsm' has to have a key with an 'X_' prefix.")

    def _are_children_of(self, component: str, column: str, ontology_name: str, ancestors: List[str]) -> bool:
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
                curie_ancestors = ONTOLOGY_CHECKER.get_term_ancestors(ontology_name, curie)
                curie_ancestors.add(curie)
                if bool(curie_ancestors & set(ancestors)):
                    return True

        return False

    def _get_raw_x(self) -> Union[np.ndarray, sparse.csc_matrix, sparse.csr_matrix]:
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
            x = self.adata.raw.X if raw_loc == "raw.X" else self.adata.X

            num_values_checked = 0
            format = get_matrix_format(self.adata, x)
            assert format != "unknown"
            self._raw_layer_exists = True
            for matrix_chunk, _, _ in self._chunk_matrix(x):
                data = matrix_chunk if isinstance(matrix_chunk, np.ndarray) else matrix_chunk.data
                if (data % 1 > 0).any():
                    self._raw_layer_exists = False
                    break

                num_values_checked += matrix_chunk.nnz if format != "dense" else np.count_nonzero(matrix_chunk)
                if num_values_checked > max_values_to_check:
                    break

        return self._raw_layer_exists

    def _validate_x_raw_x_dimensions(self):
        """
        Validates that X and raw.X have the same shape. Adds errors to self.errors if any.
        """

        if self._get_raw_x_loc() == "raw.X":
            if self.adata.n_vars != self.adata.raw.n_vars:
                self.errors.append(
                    f"Number of genes in X ({self.adata.n_vars}) is different " f"than raw.X ({self.adata.raw.n_vars})."
                )
            else:
                if not (self.adata.var.index == self.adata.raw.var.index).all():
                    self.errors.append("Index of 'raw.var' is not identical to index of 'var'.")
            if self.adata.n_obs != self.adata.raw.n_obs:
                self.errors.append(
                    f"Number of cells in X ({self.adata.n_obs}) is different " f"than raw.X ({self.adata.raw.n_obs})."
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
                            checks.append(not self._are_children_of(component, column, ontology_name, ancestors))
                    else:
                        raise ValueError(f"'{rule}' rule in raw definition of the schema is not implemented ")

        # If all checks passed then proceed with validation
        if all(checks):
            # If both "raw.X" and "X" exist but neither are raw
            # This is testing for when sometimes data contributors put a normalized matrix in both "X" and "raw.X".
            if not self._is_raw() and self._get_raw_x_loc() == "raw.X":
                self.errors.append("Raw data may be missing: data in 'raw.X' contains non-integer values.")

            # Only "X" exists but it's not raw
            # This is testing for when there is only a normalized matrix in "X" and there is no "raw.X".
            if not self._is_raw() and self._get_raw_x_loc() == "X":
                self.errors.append("Raw data is missing: there is only a normalized matrix in X and no raw.X")

            # If raw data is in X and there is nothing in raw.X (i.e. normalized values are not provided), then
            # add a warning because normalized data for RNA data is STRONGLY RECOMMENDED
            if self._is_raw() and self._get_raw_x_loc() == "X":
                self.warnings.append(
                    "Only raw data was found, i.e. there is no 'raw.X'. "
                    "It is STRONGLY RECOMMENDED that 'final' (normalized) data is provided."
                )

    def _check_single_column_availability(self, component: str, add_labels_def: List[dict]):
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

            if reserved_name in getattr_anndata(self.adata, component):
                self.errors.append(
                    f"Add labels error: Column '{reserved_name}' is a reserved column name "
                    f"of '{component}'. Remove it from h5ad and try again."
                )

    def _check_deprecated_columns(self):
        """
        This method will check for columns or keys that have been deprecated

        :rtype none
        """
        for component, component_def in self.schema_def["components"].items():
            if "deprecated_columns" in component_def:
                for column in component_def["deprecated_columns"]:
                    if column in getattr_anndata(self.adata, component):
                        self.errors.append(f"The field '{column}' is present in '{component}', but it is deprecated.")

    def _check_invalid_columns(self):
        """
        This method will check for columns or keys that are invalid
        - Columns that start with '__' (will fail the cxg conversion)

        :rtype none
        """
        for component, _ in self.schema_def["components"].items():
            df = getattr_anndata(self.adata, component)
            if df is None:
                continue
            for column in df:
                if column.startswith("__"):
                    self.errors.append(
                        f"The field '{column}' in '{component}' is invalid. Fields that start with '__' are reserved."
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
            if getattr_anndata(self.adata, component) is None:
                continue

            # Do it for columns that are forbidden
            if "forbidden_columns" in component_def:
                for column in component_def["forbidden_columns"]:
                    if column in getattr_anndata(self.adata, component):
                        self.errors.append(f"Column '{column}' must not be present in '{component}'.")

            # If ignore_labels is set, we will skip all the subsequent label checks
            if self.ignore_labels:
                continue

            # Do it for columns that map to columns
            if "columns" in component_def:
                for _column, columns_def in component_def["columns"].items():
                    if "add_labels" in columns_def:
                        self._check_single_column_availability(component, columns_def["add_labels"])

            # Do it for index that map to columns
            if "index" in component_def:
                index_def = component_def["index"]
                if "add_labels" in index_def:
                    self._check_single_column_availability(component, index_def["add_labels"])

    def _deep_check(self):
        """
        Perform a "deep" check of the AnnData object using the schema definition. Adds errors to self.errors if any

        :rtype None
        """

        # Checks for deprecated columns
        self._check_deprecated_columns()

        # Checks for invalid columns
        self._check_invalid_columns()

        # Checks that reserved columns are not used
        self._check_column_availability()

        # Checks sparsity
        logger.debug("Validating sparsity...")
        self._validate_sparsity()

        # Checks Seurat convertibility
        logger.debug("Validating Seurat convertibility...")
        self._validate_seurat_convertibility()

        # Checks each component
        for component, component_def in self.schema_def["components"].items():
            logger.debug(f"Validating component: {component}")
            # Skip if component does not exist: only useful for adata.raw.var
            if getattr_anndata(self.adata, component) is None:
                if "required" in component_def:
                    self.errors.append(f"'{component}' is missing from adata and it's required.")
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
        # existence of adata.obs["assay_ontology_term_id"]
        if not self.errors and "raw" in self.schema_def:
            logger.debug("Validating raw layer...")
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
        logger.info("Starting validation...")
        # Re-start errors in case a new h5ad is being validated
        self.errors = []

        if h5ad_path:
            logger.debug("Reading the h5ad file...")
            self.adata = read_h5ad(h5ad_path)
            self.h5ad_path = h5ad_path
            self._validate_encoding_version()
            logger.debug("Successfully read the h5ad file")

        # Fetches schema def from anndata if schema version is not found in AnnData, this fails
        self._set_schema_def()

        if not self.errors:
            self._deep_check()

        # Print warnings if any
        if self.warnings:
            self.warnings = ["WARNING: " + i for i in self.warnings]
            for w in self.warnings:
                logger.warning(w)

        # Print errors if any
        if self.errors:
            self.errors = ["ERROR: " + i for i in self.errors]
            for e in self.errors:
                logger.error(e)
            self.is_valid = False
        else:
            self.is_valid = True

        return self.is_valid


def validate(
    h5ad_path: Union[str, bytes, os.PathLike],
    add_labels_file: str = None,
    ignore_labels: bool = False,
    verbose: bool = False,
) -> (bool, list, bool):
    from .write_labels import AnnDataLabelAppender

    """
    Entry point for validation.

    :param Union[str, bytes, os.PathLike] h5ad_path: Path to h5ad file to validate
    :param str add_labels_file: Path to new h5ad file with ontology/gene labels added

    :return (True, [], <bool>) if successful validation, (False, [list_of_errors], <bool>) otherwise; last bool is for
    seurat convertibility
    :rtype tuple
    """

    # Perform validation
    start = datetime.now()
    if verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO, format="%(message)s")
    validator = Validator(ignore_labels=ignore_labels)
    validator.validate_adata(h5ad_path)
    logger.info(f"Validation complete in {datetime.now() - start} with status is_valid={validator.is_valid}")

    # Stop if validation was unsuccessful
    if not validator.is_valid:
        return False, validator.errors, validator.is_seurat_convertible

    if add_labels_file:
        label_start = datetime.now()
        writer = AnnDataLabelAppender(validator)
        writer.write_labels(add_labels_file)
        logger.info(
            f"H5AD label writing complete in {datetime.now() - label_start}, was_writing_successful: "
            f"{writer.was_writing_successful}"
        )

        return (
            validator.is_valid and writer.was_writing_successful,
            validator.errors + writer.errors,
            validator.is_seurat_convertible,
        )

    return True, validator.errors, validator.is_seurat_convertible
