import logging
import math
import numbers
import os
import re
from datetime import datetime
from typing import Dict, List, Mapping, Optional, Union

import anndata
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import scipy
from anndata._core.sparse_dataset import SparseDataset
from cellxgene_ontology_guide.ontology_parser import OntologyParser
from pandas.errors import UndefinedVariableError
from scipy import sparse

from . import gencode, schema
from .utils import SPARSE_MATRIX_TYPES, get_matrix_format, getattr_anndata, read_h5ad

logger = logging.getLogger(__name__)

ONTOLOGY_PARSER = OntologyParser(schema_version=f"v{schema.get_current_schema_version()}")

ASSAY_VISIUM = "EFO:0010961"
ASSAY_SLIDE_SEQV2 = "EFO:0030062"

VISIUM_AND_IS_SINGLE_TRUE_MATRIX_SIZE = 4992
SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE = 2000

ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE = "obs['assay_ontology_term_id'] 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True"
ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_FORBIDDEN = f"is only allowed for {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE}"
ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_REQUIRED = f"is required for {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE}"
ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_IN_TISSUE_0 = f"{ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE} and in_tissue is 0"


class Validator:
    """Handles validation of AnnData"""

    def __init__(self, ignore_labels=False):
        self.schema_def = dict()
        self.schema_version: str = None
        self.ignore_labels = ignore_labels
        self.visium_and_is_single_true_matrix_size = VISIUM_AND_IS_SINGLE_TRUE_MATRIX_SIZE

        # Values will be instances of gencode.GeneChecker,
        # keys will be one of gencode.SupportedOrganisms
        self.gene_checkers = dict()

    def reset(self):
        self.errors = []
        self.warnings = []
        self.is_valid = False
        self.h5ad_path = ""
        self._raw_layer_exists = None
        self.is_seurat_convertible: bool = True
        self.is_spatial = None
        self.is_visium = None
        self.is_visium_and_is_single_true = None

        # Matrix (e.g., X, raw.X, ...) number non-zero cache
        self.number_non_zero = dict()

    @property
    def adata(self) -> anndata.AnnData:
        return self._adata

    @adata.setter
    def adata(self, adata: anndata.AnnData):
        self.reset()
        self._adata = adata

    def _is_single(self) -> bool | None:
        """
        Determine value of uns.spatial.is_single. None if non-spatial.

        :return Value of uns.spatial.is_single if specified, None otherwise.
        :rtype bool | None
        """
        return (
            self.adata.uns["spatial"]["is_single"]
            if hasattr(self.adata, "uns")
            and "spatial" in self.adata.uns
            and isinstance(self.adata.uns["spatial"], dict)
            and "is_single" in self.adata.uns["spatial"]
            else None
        )

    def _is_supported_spatial_assay(self) -> bool:
        """
        Determine if the assay_ontology_term_id is either Visium (EFO:0010961) or Slide-seqV2 (EFO:0030062).

        :return True if assay_ontology_term_id is Visium or Slide-seqV2, False otherwise.
        :rtype bool
        """
        if self.is_spatial is None:
            try:
                self.is_spatial = False
                if self.adata.obs.assay_ontology_term_id.isin([ASSAY_VISIUM, ASSAY_SLIDE_SEQV2]).any():
                    self.is_spatial = True
            except AttributeError:
                # specific error reporting will occur downstream in the validation
                self.is_spatial = False
        return self.is_spatial

    def _is_visium_and_is_single_true(self) -> bool:
        """
        Determine if the assay_ontology_term_id is Visium (EFO:0010961) and uns.spatial.is_single is True.

        :return True if assay_ontology_term_id is Visium and is_single_cell is True, False otherwise.
        :rtype bool
        """
        if self.is_visium_and_is_single_true is None:
            self.is_visium_and_is_single_true = bool(self._is_visium() and self._is_single())
        return self.is_visium_and_is_single_true

    def _validate_encoding_version(self):
        import h5py

        with h5py.File(self.h5ad_path, "r") as f:
            encoding_dict = dict(f.attrs)
            encoding_version = encoding_dict.get("encoding-version")
            if encoding_version != "0.1.0":
                self.errors.append("The h5ad artifact was generated with an AnnData version different from 0.8.0.")

    def _set_schema_def(self):
        """
        Sets schema dictionary
        """
        if not self.schema_version:
            self.schema_version = schema.get_current_schema_version()
        if not self.schema_def:
            self.schema_def = schema.get_schema_definition()

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

    def _has_forbidden_curie_ancestor(
        self, term_id: str, column_name: str, forbidden_def: Dict[str, List[str]]
    ) -> bool:
        """
        Validate if a single curie term id is a descendant term of any forbidden ancestors.
        If there is a forbidden ancestor detected, it adds it to self.errors.

        :param str term_id: the curie term id to validate
        :param str column_name: original column name in adata where the term_id comes from (used for error messages)
        :param Dict[str, List[str] forbidden_def: mapping of ontologies to list of ancestor terms to validate against

        :returns bool
        """
        for ontology_name in forbidden_def:
            for ancestor in forbidden_def[ontology_name]:
                if ancestor in ONTOLOGY_PARSER.get_term_ancestors(term_id):
                    self.errors.append(
                        f"'{term_id}' in '{column_name}' is not allowed. Descendant terms of "
                        f"'{ancestor}' are not allowed."
                    )
                    return True
        return False

    def _validate_curie_ancestors(
        self,
        term_id: str,
        allowed_ancestors: Dict[str, List[str]],
        inclusive: bool = False,
    ) -> bool:
        """
        Validate a single curie term id is a valid descendant of any allowed ancestors

        :param str term_id: the curie term id to validate
        :param dict{str: list[str]} allowed_ancestors: keys must be ontology names and values must lists of
        allowed ancestors
        :param bool inclusive:  if True then the ancestors themselves are allowed

        :rtype Bool
        """

        checks = []

        for _, ancestors in allowed_ancestors.items():
            for ancestor in ancestors:
                if inclusive and term_id == ancestor:
                    checks.append(True)

                is_valid_term_id = ONTOLOGY_PARSER.is_valid_term_id(term_id)
                is_valid_ancestor_id = ONTOLOGY_PARSER.is_valid_term_id(ancestor)
                if is_valid_term_id & is_valid_ancestor_id:
                    is_descendant = ancestor in ONTOLOGY_PARSER.get_term_ancestors(term_id)
                    checks.append(is_descendant)

        if True not in checks:
            return False
        return True

    def _validate_curie_ontology(self, term_id: str, column_name: str, allowed_ontologies: List[str]) -> bool:
        """
        Validate a single curie term id belongs to specified ontologies. If it does belong to an allowed ontology
        verifies that it is not deprecated (obsolete).
        If there are any errors, it adds them to self.errors

        :param str term_id: the curie term id to validate
        :param str column_name: original column name in adata where the term_id comes from (used for error messages)
        :param List[str] allowed_ontologies: allowed ontologies

        :rtype bool
        """

        checks = []

        for ontology_name in allowed_ontologies:
            try:
                is_valid = ONTOLOGY_PARSER.is_valid_term_id(term_id, ontology_name)
            except ValueError:
                is_valid = False
            checks.append(is_valid)

            if is_valid and ONTOLOGY_PARSER.is_term_deprecated(term_id):
                self.errors.append(f"'{term_id}' in '{column_name}' is a deprecated term id of '{ontology_name}'.")
                return False

        if sum(checks) == 0:
            self.errors.append(
                f"'{term_id}' in '{column_name}' is not a valid ontology term id of '{', '.join(allowed_ontologies)}'."
            )
            return False
        return True

    def _validate_curie_str(self, term_str: str, column_name, curie_constraints: dict) -> None:
        """
        Validate a curie str based on some constraints. If there are any errors, it adds them to self.errors

        :param str term_str: the curie term str to validate
        :param str column_name: Name of the column in the dataframe
        :param dict curie_constraints: constraints for the curie term to be validated,
        this part of the schema definition

        :rtype None
        """
        if "exceptions" in curie_constraints and term_str in curie_constraints["exceptions"]:
            return

        # If NA is found in allowed ontologies, it means only exceptions should be found. If no exceptions were found
        # then return error
        if curie_constraints["ontologies"] == ["NA"]:
            self.errors.append(f"'{term_str}' in '{column_name}' is not a valid value of '{column_name}'.")
            return

        if not isinstance(term_str, str):
            self.errors.append(
                f"'{term_str}' in '{column_name}' is not a valid ontology term value, it must be a string."
            )
            return

        # if multi_term is defined, split str into individual ontology terms and validate each
        if "multi_term" in curie_constraints:
            delimiter = curie_constraints["multi_term"]["delimiter"]
            term_ids = term_str.split(delimiter)
            for term_id in term_ids:
                self._validate_curie(term_id, column_name, curie_constraints)
            if (
                curie_constraints["multi_term"].get("sorted", False)
                and len(term_ids) > 1
                and not all(term_ids[i].strip() <= term_ids[i + 1].strip() for i in range(len(term_ids) - 1))
            ):
                self.errors.append(f"'{term_str}' in '{column_name}' is not in ascending lexical order.")
            if len(set(term_ids)) != len(term_ids):
                self.errors.append(f"'{term_str}' in '{column_name}' contains duplicates.")
        else:
            self._validate_curie(term_str, column_name, curie_constraints)

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
        # If the term id does not belong to an allowed ontology, the subsequent checks are redundant
        if not self._validate_curie_ontology(term_id, column_name, curie_constraints["ontologies"]):
            return

        # Check if term_id is forbidden by schema definition despite being a valid ontology term
        if "forbidden" in curie_constraints:
            if "terms" in curie_constraints["forbidden"] and term_id in curie_constraints["forbidden"]["terms"]:
                self.errors.append(f"'{term_id}' in '{column_name}' is not allowed.")
                return
            if "ancestors" in curie_constraints["forbidden"] and self._has_forbidden_curie_ancestor(
                term_id, column_name, curie_constraints["forbidden"]["ancestors"]
            ):
                return

        # If there are allow-lists, validate against them
        if "allowed" in curie_constraints:
            is_allowed = False
            if "terms" in curie_constraints["allowed"]:
                for ontology_name, allow_list in curie_constraints["allowed"]["terms"].items():
                    if ONTOLOGY_PARSER.is_valid_term_id(term_id, ontology_name) and (
                        allow_list == ["all"] or term_id in set(allow_list)
                    ):
                        is_allowed = True
                        break
            if (
                not is_allowed
                and "ancestors" in curie_constraints["allowed"]
                and self._validate_curie_ancestors(term_id, curie_constraints["allowed"]["ancestors"])
            ):
                is_allowed = True

            if not is_allowed:
                self.errors.append(f"'{term_id}' in '{column_name}' is not an allowed term id.")

    def _validate_feature_id(self, feature_id: str, df_name: str):
        """
        Validates a feature id, i.e. checks that it's present in the reference
        If there are any errors, it adds them to self.errors and adds it to the list of invalid features

        :param str feature_id: the feature id to be validated
        :param str df_name: name of dataframe the feauter id comes from (var or raw.var)

        :rtype none
        """

        organism = gencode.get_organism_from_feature_id(feature_id)

        if not organism:
            self.errors.append(
                f"Could not infer organism from feature ID '{feature_id}' in '{df_name}', "
                f"make sure it is a valid ID."
            )
            return

        if organism not in self.gene_checkers:
            self.gene_checkers[organism] = gencode.GeneChecker(organism)

        if not self.gene_checkers[organism].is_valid_id(feature_id):
            self.errors.append(f"'{feature_id}' is not a valid feature ID in '{df_name}'.")

        return

    @staticmethod
    def _chunk_matrix(
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
            logger.debug(f"_chunk_matrix [{i} of {math.ceil(n / obs_chunk_size)}]")
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
        matrix_format = get_matrix_format(self.adata, matrix)
        for matrix_chunk, _, _ in self._chunk_matrix(matrix):
            nnz += matrix_chunk.count_nonzero() if matrix_format != "dense" else np.count_nonzero(matrix_chunk)

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
            if X_format in SPARSE_MATRIX_TYPES:
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

            if "curie_constraints" in column_def:
                for term_str in column.drop_duplicates():
                    self._validate_curie_str(term_str, column_name, column_def["curie_constraints"])

        # Add error suffix to errors found here
        if "error_message_suffix" in column_def:
            error_total_count = len(self.errors)
            for i in range(error_original_count, error_total_count):
                self.errors[i] = self.errors[i] + " " + column_def["error_message_suffix"]

    def _validate_column_dependencies(
        self, df: pd.DataFrame, df_name: str, column_name: str, dependencies: List[dict]
    ) -> pd.Series:
        """
        Validates subset of columns based on dependencies, for instance development_stage_ontology_term_id has
        dependencies with organism_ontology_term_id -- the allowed values depend on whether organism is human, mouse
        or something else.

        After performing all validations, it will return the column that has been stripped of the already validated
        fields -- this has to still be validated.

        :param pd.DataFrame df: pandas dataframe containing the column to be validated
        :param str df_name: the name of dataframe in the adata object, e.g. "obs"
        :param str column_name: the name of the column to be validated
        :param list dependencies: a list of dependency definitions, which is a list of column definitions with a "rule"
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
        a specified column is a descendant term to a group of specified ancestors, returning a Bool.
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
            return validate_curie_ancestors_vectorized(term_id, allowed_ancestors, inclusive=ancestor_inclusive)

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
                if value_def.get("required", False):
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

        for column_name in df.columns:
            column = df[column_name]
            if column.dtype.name != "category":
                # Check for columns with mixed values, which is not supported by anndata 0.8.0
                # TODO: check if this can be removed after upgading to anndata 0.10.0
                value_types = {type(x) for x in column.values}
                if len(value_types) != 1:
                    self.errors.append(
                        f"Column '{column_name}' in dataframe '{df_name}' cannot contain mixed types. Found {value_types}."
                    )
            else:
                # Check for columns that have a category defined 0 times (obs only)
                if df_name == "obs":
                    for category in column.dtype.categories:
                        if category not in column.values:
                            self.warnings.append(
                                f"Column '{column_name}' in dataframe '{df_name}' contains a category '{category}' with "
                                f"zero observations. These categories will be removed when `--add-labels` flag is present."
                            )
                categorical_types = {type(x) for x in column.dtype.categories.values}
                # Check for columns that have illegal categories, which are not supported by anndata 0.8.0
                # TODO: check if this can be removed after upgading to anndata 0.10.0
                blocked_categorical_types = {bool}
                illegal_categorical_types = categorical_types & blocked_categorical_types
                if illegal_categorical_types:
                    self.errors.append(
                        f"Column '{column_name}' in dataframe '{df_name}' contains {illegal_categorical_types=}."
                    )
                # Check for categorical column has mixed types, which is not supported by anndata 0.8.0
                # TODO: check if this can be removed after upgading to anndata 0.10.0
                categorical_types = {type(x) for x in column.dtype.categories.values}
                if len(categorical_types) > 1:
                    self.errors.append(
                        f"Column '{column_name}' in dataframe '{df_name}' contains {len(categorical_types)} categorical types. "
                        f"Only one type is allowed."
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

    def _validate_uns_dict(self, uns_dict: dict) -> None:
        df = getattr_anndata(self.adata, "obs")

        # Mapping from obs column name to number of unique categorical values
        category_mapping = {}

        # Check for categorical dtypes in the dataframe directly
        for column_name in df.columns:
            column = df[column_name]
            if column.dtype.name == "category":
                category_mapping[column_name] = column.nunique()

        for key, value in uns_dict.items():
            if any(
                isinstance(value, sparse_class)
                for sparse_class in (scipy.sparse.csr_matrix, scipy.sparse.csc_matrix, scipy.sparse.coo_matrix)
            ):
                if value.nnz == 0:  # number non-zero
                    self.errors.append(f"uns['{key}'] cannot be an empty value.")
            elif value is not None and not isinstance(value, (np.bool_, bool, numbers.Number)) and len(value) == 0:
                self.errors.append(f"uns['{key}'] cannot be an empty value.")
            if key.endswith("_colors"):
                # 1. Verify that the corresponding categorical field exists in obs
                column_name = key.replace("_colors", "")
                obs_unique_values = category_mapping.get(column_name)
                if not obs_unique_values:
                    error_message = f"Colors field uns[{key}] does not have a corresponding categorical field in obs"
                    if column_name in df.columns:
                        error_message += f". {column_name} is present but is dtype {df[column_name].dtype.name}"
                    portal_column_names = [
                        "assay",
                        "cell_type",
                        "development_stage",
                        "disease",
                        "organism",
                        "self_reported_ethnicity",
                        "sex",
                        "tissue",
                    ]
                    if column_name in portal_column_names:
                        error_message += f". Annotate {column_name}_ontology_term_id_colors instead"
                    self.errors.append(error_message)
                    continue
                # 2. Verify that the value is a numpy array
                if value is None or not isinstance(value, np.ndarray):
                    self.errors.append(
                        f"Colors field uns['{key}'] must be of 'numpy.ndarray' type, it is {type(value)}"
                    )
                    # Skip over all subsequent validations which expect a numpy array
                    continue
                # 3. Verify that we have strings in the array
                all_strings = all(isinstance(color, str) for color in value)
                if not all_strings:
                    self.errors.append(
                        f"Colors in uns[{key}] must be strings. Found: {value} which are {value.dtype.name}"
                    )
                    continue
                # 4. Verify that we have at least as many colors as unique values in the corresponding categorical field
                if len(value) < obs_unique_values:
                    self.errors.append(
                        f"Annotated categorical field {key.replace('_colors', '')} must have at least {obs_unique_values} color options "
                        f"in uns[{key}]. Found: {value}"
                    )
                # 5. Verify that either all colors are hex OR all colors are CSS4 named colors strings
                all_hex_colors = all(re.match(r"^#([0-9a-fA-F]{6})$", color) for color in value)
                all_css4_colors = all(color in mcolors.CSS4_COLORS for color in value)
                if not (all_hex_colors or all_css4_colors):
                    self.errors.append(
                        f"Colors in uns[{key}] must be either all hex colors or all CSS4 named colors. Found: {value}"
                    )

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
            assert matrix_format != "unknown"

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
        # Seurat conversion is not supported for Visium datasets.
        if self._is_visium():
            self.warnings.append(
                "Datasets with assay_ontology_term_id 'EFO:0010961' (Visium Spatial Gene Expression) are not compatible with Seurat."
            )
            self.is_seurat_convertible = False
            return

        to_validate = [(self.adata.X, "X")]
        # check if there's raw data
        if self.adata.raw:
            to_validate.append((self.adata.raw.X, "raw.X"))
        # Check length of component arrays
        for matrix, matrix_name in to_validate:
            matrix_format = get_matrix_format(self.adata, matrix)
            if matrix_format in SPARSE_MATRIX_TYPES:
                effective_r_array_size = self._count_matrix_nonzero(matrix_name, matrix)
                is_sparse = True
            elif matrix_format == "dense":
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

        if self.adata.raw and self.adata.raw.X.shape[1] != self.adata.raw.var.shape[0]:
            self.errors.append(
                "This dataset has a mismatch between 1) the number of features in raw.X and 2) the number of features "
                "in raw.var. These counts must be identical."
            )
            self.is_seurat_convertible = False

    def _validate_obsm(self):
        """
        Validates the embedding dictionary -- it checks that all values of adata.obsm are numpy arrays with the correct
        dimension. Adds errors to self.errors if any. Checks that the keys start with "X_", have no whitespace, and have
        a suffix at least 1 character long. For keys that don't start with "X_", we will run them through the same
        validation checks, but raise warnings instead of errors.

        :rtype none
        """

        if not self.adata.obsm:
            self.errors.append("No embeddings found in 'adata.obsm'.")
            return

        obsm_with_x_prefix = 0
        for key, value in self.adata.obsm.items():
            issue_list = self.errors

            regex_pattern = r"^[a-zA-Z][a-zA-Z0-9_.-]*$"

            unknown_key = False  # an unknown key does not match 'spatial' or 'X_{suffix}'
            if key.startswith("X_"):
                obsm_with_x_prefix += 1
                if key.lower() == "x_spatial":  # TODO undo after 5.0 patch release
                    self.errors.append(f"Embedding key in 'adata.obsm' {key} cannot be used.")
                elif not re.match(regex_pattern, key[2:]):
                    self.errors.append(
                        f"Suffix for embedding key in 'adata.obsm' {key} does not match the regex pattern {regex_pattern}."
                    )
            elif key.lower() != "spatial":
                if not re.match(regex_pattern, key):
                    self.errors.append(
                        f"Embedding key in 'adata.obsm' {key} does not match the regex pattern {regex_pattern}."
                    )
                self.warnings.append(
                    f"Embedding key in 'adata.obsm' {key} is not 'spatial' nor does it start with 'X_'. Thus, it will "
                    f"not be available in Explorer"
                )
                issue_list = self.warnings
                unknown_key = True

            if not isinstance(value, np.ndarray):
                self.errors.append(
                    f"All embeddings have to be of 'numpy.ndarray' type, " f"'adata.obsm['{key}']' is {type(value)}')."
                )
                # Skip over the subsequent checks that require the value to be an array
                continue

            if len(value.shape) < 2:
                self.errors.append(
                    f"All embeddings must at least two dimensions. 'adata.obsm['{key}']' has a shape length of '{len(value.shape)}'."
                )
            else:
                if value.shape[0] != self.adata.n_obs:
                    self.errors.append(
                        f"All embeddings must have as many rows as cells. 'adata.obsm['{key}']' has rows='{value.shape[0]}'."
                    )

                if unknown_key and value.shape[1] < 1:
                    self.errors.append(
                        f"All unspecified embeddings must have at least one column. 'adata.obsm['{key}']' has columns='{value.shape[1]}'."
                    )

                if not unknown_key and value.shape[1] < 2:
                    self.errors.append(
                        f"All 'X_' and 'spatial' embeddings must have at least two columns. 'adata.obsm['{key}']' has columns='{value.shape[1]}'."
                    )

            if not (np.issubdtype(value.dtype, np.integer) or np.issubdtype(value.dtype, np.floating)):
                issue_list.append(
                    f"adata.obsm['{key}'] has an invalid data type. It should be "
                    "float, integer, or unsigned integer of any precision (8, 16, 32, or 64 bits)."
                )
            else:
                # Check for inf/NaN values only if the dtype is numeric
                if np.isinf(value).any():
                    issue_list.append(f"adata.obsm['{key}'] contains positive infinity or negative infinity values.")
                if np.all(np.isnan(value)):
                    issue_list.append(f"adata.obsm['{key}'] contains all NaN values.")

        if self._is_supported_spatial_assay() is False and obsm_with_x_prefix == 0:
            self.errors.append("At least one embedding in 'obsm' has to have a key with an 'X_' prefix.")

        is_single = self._is_single()
        has_spatial_embedding = "spatial" in self.adata.obsm
        if is_single and not has_spatial_embedding:
            self.errors.append(
                "'spatial' embedding is required in 'adata.obsm' if " "adata.uns['spatial']['is_single'] is True."
            )
        elif is_single is None and has_spatial_embedding:
            self.errors.append(
                "'spatial' embedding is forbidden in 'adata.obsm' if " "adata.uns['spatial']['is_single'] is not set."
            )

    def _validate_annotation_mapping(self, component_name: str, component: Mapping):
        for key, value in component.items():
            # Check for empty ndarrays
            if isinstance(value, np.ndarray) and not value.size:
                self.errors.append(
                    f"The size of the ndarray stored for a 'adata.{component_name}['{key}']' MUST NOT be zero."
                )

    def _are_descendants_of(self, component: str, column: str, ontology_name: str, ancestors: List[str]) -> bool:
        """
        Checks if elements in the specified column of the component (e.g. 'assay_ontology_term_id' of 'adata.obs') are
        descendants of the given ancestors.

        Ancestors checks are inclusive, meaning that a value is its own ancestor as well.

        :param str component: the name of the component that's been checked.
        :param str column: Column in the component to check
        :param str ontology_name: Name of the ontology (e.g. "EFO")
        :param List[str] ancestors: List of ancestors

        :rtype bool
        :return True if any value in column is a descendant of any ancestor.
        """

        curies = getattr(getattr(self.adata, component), column)
        curies = curies.drop_duplicates()

        for curie in curies:
            if ONTOLOGY_PARSER.is_valid_term_id(curie):
                curie_ancestors = ONTOLOGY_PARSER.get_term_ancestors(curie, include_self=True)
                if bool(set(curie_ancestors) & set(ancestors)):
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

    def _has_valid_raw(self, force: bool = False) -> bool:
        """
        Checks if the non-zero values for the raw matrix (adata.X or adata.raw.X)
        are positive integers stored as numpy.float32. Also validates that every row contains at least one non-zero
        value.

        Returns False if at least one value / row does not meet requirements. True otherwise.

        Since this process is memory intensive, it will return a cache value if this function has been called before.
        If calculation needs to be repeated use `force = True`

        :rtype bool
        """
        if force:
            self._raw_layer_exists = None

        if self._raw_layer_exists is None:
            # Get potential raw_X
            x = self._get_raw_x()
            if x.dtype != np.float32:
                self._raw_layer_exists = False
                self.errors.append("Raw matrix values must have type numpy.float32.")
                return self._raw_layer_exists

            matrix_format = get_matrix_format(self.adata, x)
            assert matrix_format != "unknown"
            self._raw_layer_exists = True
            is_sparse_matrix = matrix_format in SPARSE_MATRIX_TYPES

            is_visium_and_is_single_true = self._is_visium_and_is_single_true()
            if is_visium_and_is_single_true and x.shape[0] != self.visium_and_is_single_true_matrix_size:
                self._raw_layer_exists = False
                self.errors.append(
                    f"When {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE}, the raw matrix must be the "
                    f"unfiltered feature-barcode matrix 'raw_feature_bc_matrix'. It must have exactly "
                    f"{self.visium_and_is_single_true_matrix_size} rows. Raw matrix row count is "
                    f"{x.shape[0]}."
                )

            if (
                is_visium_and_is_single_true
                and "in_tissue" in self.adata.obs
                and 0 in self.adata.obs["in_tissue"].values
            ):
                self._validate_raw_data_with_in_tissue_0(x, is_sparse_matrix)
            else:
                self._validate_raw_data(x, is_sparse_matrix)

        return self._raw_layer_exists

    def _validate_raw_data(self, x: Union[np.ndarray, sparse.spmatrix], is_sparse_matrix: bool):
        """
        Validates the data values in the raw matrix. Matrix size is chunked for large matrices.

        :param x: raw matrix
        :param is_sparse_matrix: bool indicating if the matrix is sparse {csc, csr, coo}
        """
        has_row_of_zeros = False
        has_invalid_nonzero_value = False
        for matrix_chunk, _, _ in self._chunk_matrix(x):
            if not has_row_of_zeros:
                if is_sparse_matrix:
                    row_indices, _ = matrix_chunk.nonzero()
                    if len(set(row_indices)) != matrix_chunk.shape[0]:
                        has_row_of_zeros = True
                # else, must be dense matrix, confirm that all rows have at least 1 nonzero value
                elif not all(np.apply_along_axis(np.any, axis=1, arr=matrix_chunk)):
                    has_row_of_zeros = True

            if not has_invalid_nonzero_value and self._matrix_has_invalid_nonzero_values(matrix_chunk):
                has_invalid_nonzero_value = True

            if has_row_of_zeros and has_invalid_nonzero_value:
                # Fail fast, exit loop and report
                break

        if has_row_of_zeros:
            self._raw_layer_exists = False
            self.errors.append("Each cell must have at least one non-zero value in its row in the raw matrix.")
        if has_invalid_nonzero_value:
            self._raw_layer_exists = False
            self.errors.append("All non-zero values in raw matrix must be positive integers of type numpy.float32.")

    def _validate_raw_data_with_in_tissue_0(
        self, x: Union[np.ndarray, sparse.spmatrix, SparseDataset], is_sparse_matrix: bool
    ):
        """
        Special case validation checks for Visium data with is_single = True and in_tissue column in obs where in_tissue
        has at least one value 0. Static matrix size of 4992 rows, so chunking is not required.

        :param x: raw matrix
        :param is_sparse_matrix: bool indicating if the matrix is sparse {csc, csr, coo}
        """
        has_tissue_0_non_zero_row = False
        has_tissue_1_zero_row = False
        if isinstance(x, SparseDataset):
            x = x.to_memory()
        if is_sparse_matrix:
            nonzero_row_indices, _ = x.nonzero()
        else:  # must be dense matrix
            nonzero_row_indices = np.where(np.any(x != 0, axis=1))[0]
        for i in range(x.shape[0]):
            if not has_tissue_0_non_zero_row and i in nonzero_row_indices and self.adata.obs["in_tissue"].iloc[i] == 0:
                has_tissue_0_non_zero_row = True
            elif (
                not has_tissue_1_zero_row and i not in nonzero_row_indices and self.adata.obs["in_tissue"].iloc[i] == 1
            ):
                has_tissue_1_zero_row = True
            if has_tissue_0_non_zero_row and has_tissue_1_zero_row:
                # exit early and report
                break

        if not has_tissue_0_non_zero_row:
            self._raw_layer_exists = False
            self.errors.append(
                "If obs['in_tissue'] contains at least one value 0, then there must be at least "
                "one row with obs['in_tissue'] == 0 that has a non-zero value in the raw matrix."
            )
        if has_tissue_1_zero_row:
            self._raw_layer_exists = False
            self.errors.append(
                "Each observation with obs['in_tissue'] == 1 must have at least one "
                "non-zero value in its row in the raw matrix."
            )
        if self._matrix_has_invalid_nonzero_values(x):
            self._raw_layer_exists = False
            self.errors.append("All non-zero values in raw matrix must be positive integers of type numpy.float32.")

    def _matrix_has_invalid_nonzero_values(self, x: Union[np.ndarray, sparse.spmatrix]) -> bool:
        """
        Checks whether the matrix has invalid non-zero values. The matrix must have all non-zero values as positive
        integers (type is numpy.float32).

        :param x: The matrix to validate

        :rtype bool
        :return True if any non-zero values are invalid given validation rules, False otherwise
        """
        data = x if isinstance(x, np.ndarray) else x.data
        return np.any((data % 1 > 0) | (data < 0))

    def _validate_x_raw_x_dimensions(self):
        """
        Validates that X and raw.X have the same shape, if raw.X exists. Adds errors to self.errors if any.
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
         1. X and raw.X have the same column and row indices
         2. there's an expression matrix containing raw (integer) values, first in adata.raw.X and then adata.X if
         the former does not exist.
         3. For applicable assays, checks that each row has at least one non-zero value

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
                    if rule == "not_descendants_of":
                        for ontology_name, ancestors in rule_def.items():
                            checks.append(not self._are_descendants_of(component, column, ontology_name, ancestors))
                    else:
                        raise ValueError(f"'{rule}' rule in raw definition of the schema is not implemented ")

        # If all checks passed then proceed with validation
        if all(checks):
            # If both "raw.X" and "X" exist but neither are raw
            # This is testing for when sometimes data contributors put a normalized matrix in both "X" and "raw.X".
            if not self._has_valid_raw() and self._get_raw_x_loc() == "raw.X":
                self.errors.append("Raw data may be missing: data in 'raw.X' does not meet schema requirements.")

            # Only "X" exists but it's not raw
            # This is testing for when there is only a normalized matrix in "X" and there is no "raw.X".
            if not self._has_valid_raw() and self._get_raw_x_loc() == "X":
                self.errors.append("Raw data is missing: there is only a normalized matrix in X and no raw.X")

            # If raw data is in X and there is nothing in raw.X (i.e. normalized values are not provided), then
            # add a warning because normalized data for RNA data is STRONGLY RECOMMENDED
            if self._has_valid_raw() and self._get_raw_x_loc() == "X":
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

    def _check_var_and_obs_column_name_uniqueness(self):
        """
        This method checks that all column names in the 'var' and 'obs' DataFrames are unique

        :rtype none
        """
        dataframe_components = ["obs", "var", "raw.var"]
        for df_component in dataframe_components:
            adata_component = getattr_anndata(self.adata, df_component)
            if adata_component is None:
                continue
            component_columns = set()
            for column in adata_component.columns:
                if column in component_columns:
                    raise ValueError(
                        f"Duplicate column name '{column}' detected in 'adata.{df_component}' DataFrame. All DataFrame column names must be unique."
                    )
                component_columns.add(column)

    def _check_column_availability(self):
        """
        This method will check for columns that are reserved in components and validate that they are
         available as expected

        :rtype none
        """

        for component, component_def in self.schema_def["components"].items():
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

            # Do it for metadata columns that are reserved for annotation after data portal upload
            if "reserved_columns" in component_def:
                for column in component_def["reserved_columns"]:
                    if column in getattr_anndata(self.adata, component):
                        self.errors.append(
                            f"Column '{column}' is a reserved column name "
                            f"of '{component}'. Remove it from h5ad and try again."
                        )

            # Do it for columns that map to other columns, for post-upload annotation
            if "columns" in component_def:
                for column_def in component_def["columns"].values():
                    if "add_labels" in column_def:
                        self._check_single_column_availability(component, column_def["add_labels"])

            # Do it for index that map to columns
            if "index" in component_def:
                index_def = component_def["index"]
                if "add_labels" in index_def:
                    self._check_single_column_availability(component, index_def["add_labels"])

    def _check_spatial(self):
        """
        Sequence validation of spatial-related values of the AnnData object.

        :rtype none
        """
        self._check_spatial_uns()
        self._check_spatial_obs()

    def _check_spatial_obs(self):
        """
        Validate obs spatial-related values of the AnnData object. Validation is not defined in schema definition yaml.
        Errors are added to self.errors.

        :rtype none
        """

        # Exit if obs is not specified. Error is reported in core validate functionality.
        obs_component = getattr_anndata(self.adata, "obs")
        if obs_component is None:
            return

        # Validate assay ontology term id.
        self._validate_spatial_assay_ontology_term_id()

        # Validate tissue positions.
        self._validate_spatial_tissue_positions()

        # Validate cell type.
        self._validate_spatial_cell_type_ontology_term_id()

        self._validate_spatial_is_primary_data()

    def _validate_spatial_is_primary_data(self):
        """
        Validate is_primary_data for spatial datasets.
        """
        obs = getattr_anndata(self.adata, "obs")
        if obs is None or "is_primary_data" not in obs:
            return
        if self._is_single() is False and obs["is_primary_data"].any():
            self.errors.append(
                "When uns['spatial']['is_single'] is False, obs['is_primary_data'] must be False for all rows."
            )

    def _validate_spatial_assay_ontology_term_id(self):
        """
        If assay is spatial, all assay ontology term ids should be identical.

        :rtype none
        """
        # Identical check requires assay_ontology_term_id.
        obs = getattr_anndata(self.adata, "obs")
        if obs is None or "assay_ontology_term_id" not in obs:
            return

        # Identical check is only applicable to spatial datasets.
        if not self._is_supported_spatial_assay():
            return

        # Validate assay ontology term ids are identical.
        term_count = obs["assay_ontology_term_id"].nunique()
        if term_count > 1:
            self.errors.append(
                "When obs['assay_ontology_term_id'] is either 'EFO:0010961' (Visium Spatial Gene Expression) or "
                "'EFO:0030062' (Slide-seqV2), all observations must contain the same value."
            )

    def _validate_spatial_cell_type_ontology_term_id(self):
        """
        Validate cell type ontology term id is "unknown" if Visium, is_single is True and in_tissue is 0.

        :rtype none
        """
        # Exit if:
        # - not Visium and is_single is True as no further checks are necessary
        # - in_tissue is not specified as checks are dependent on this value
        if not self._is_visium_and_is_single_true() or "in_tissue" not in self.adata.obs:
            return

        # Validate cell type: must be "unknown" if Visium and is_single is True and in_tissue is 0.
        if (
            (self.adata.obs["assay_ontology_term_id"] == ASSAY_VISIUM)
            & (self.adata.obs["in_tissue"] == 0)
            & (self.adata.obs["cell_type_ontology_term_id"] != "unknown")
        ).any():
            self.errors.append(
                f"obs['cell_type_ontology_term_id'] must be 'unknown' when {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_IN_TISSUE_0}."
            )

    def _validate_spatial_tissue_position(self, tissue_position_name: str, min: int, max: int):
        """
        Validate tissue position is allowed and required, and are integers within the given range. Validation is not defined in
        schema definition yaml.

        :rtype none
        """
        # Tissue position is foribidden if assay is not Visium and is_single is True.
        if tissue_position_name in self.adata.obs and (
            not self._is_visium_and_is_single_true()
            or (
                ~(self.adata.obs["assay_ontology_term_id"] == ASSAY_VISIUM)
                & (self.adata.obs[tissue_position_name].notnull())
            ).any()
        ):
            self.errors.append(f"obs['{tissue_position_name}'] {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_FORBIDDEN}.")
            return

        # Exit if we're not dealing with Visium and _is_single True as no further checks are necessary.
        if not self._is_visium_and_is_single_true():
            return

        # At this point, is_single is True and:
        # - there's at least one row with Visum, tissue position column is required
        # - for any Visium row, tissue position is required.
        if (
            tissue_position_name not in self.adata.obs
            or (
                (self.adata.obs["assay_ontology_term_id"] == ASSAY_VISIUM)
                & (self.adata.obs[tissue_position_name].isnull())
            ).any()
        ):
            self.errors.append(f"obs['{tissue_position_name}'] {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_REQUIRED}.")
            return

        # Tissue position must be an int.
        obs_tissue_position = self.adata.obs.get(tissue_position_name)
        if not np.issubdtype(obs_tissue_position.dtype, np.integer):
            self.errors.append(f"obs['{tissue_position_name}'] must be of int type, it is {obs_tissue_position.dtype}.")
            return

        # Tissue position must be within the given range.
        if not ((obs_tissue_position >= min) & (obs_tissue_position <= max)).all():
            if tissue_position_name == "in_tissue":
                error_message_token = f"{min} or {max}"
            else:
                error_message_token = f"between {min} and {max}"
            self.errors.append(
                f"obs['{tissue_position_name}'] must be {error_message_token}, the min and max are {obs_tissue_position.min()} and {obs_tissue_position.max()}. "
                f"This must be the value of the column tissue_positions_in_tissue from the tissue_positions_list.csv or tissue_positions.csv."
            )

    def _validate_spatial_tissue_positions(self):
        """
        Validate tissue positions of spatial datasets.

        :rtype none
        """
        self._validate_spatial_tissue_position("array_col", 0, 127)
        self._validate_spatial_tissue_position("array_row", 0, 77)
        self._validate_spatial_tissue_position("in_tissue", 0, 1)

    def _check_spatial_uns(self):
        """
        Validate uns spatial-related values of the AnnData object. Validation is not defined in schema definition yaml.
        Errors are added to self.errors.

        :rtype none
        """

        # Exit if uns is not specified. Error is reported in core validate functionality.
        uns_component = getattr_anndata(self.adata, "uns")
        if uns_component is None:
            return

        # uns spatial validation is dependent on obs.assay_ontology_term_id; exit if not specified. Error is
        # reported in core validate functionality.
        obs_component = getattr_anndata(self.adata, "obs")
        if obs_component is None or "assay_ontology_term_id" not in obs_component:
            return

        # spatial is forbidden if assay it not a supported spatial assay.
        uns_spatial = self.adata.uns.get("spatial")
        is_supported_spatial_assay = self._is_supported_spatial_assay()
        if uns_spatial is not None and not is_supported_spatial_assay:
            self.errors.append(
                "uns['spatial'] is only allowed for obs['assay_ontology_term_id'] values "
                "'EFO:0010961' (Visium Spatial Gene Expression) and 'EFO:0030062' (Slide-seqV2)."
            )
            return

        # Exit if we aren't dealing with a supported spatial assay as no further checks are necessary.
        if not is_supported_spatial_assay:
            return

        # spatial is required for supported spatial assays.
        if not isinstance(uns_spatial, dict):
            self.errors.append(
                "A dict in uns['spatial'] is required for obs['assay_ontology_term_id'] values "
                "'EFO:0010961' (Visium Spatial Gene Expression) and 'EFO:0030062' (Slide-seqV2)."
            )
            return

        # is_single is required.
        if "is_single" not in uns_spatial:
            self.errors.append("uns['spatial'] must contain the key 'is_single'.")
            # Exit if is_single is missing as all further checks are dependent on its value.
            return

        # is_single must be a boolean.
        uns_is_single = uns_spatial["is_single"]
        if not isinstance(uns_is_single, (np.bool_, bool)):
            self.errors.append(f"uns['spatial']['is_single'] must be of boolean type, it is {type(uns_is_single)}.")
            # Exit if is_single is not valid as all further checks are dependent on its value.
            return

        # Check there is at most one library_id.
        uns_spatial_keys = list(uns_spatial.keys())
        library_ids = list(filter(lambda x: x != "is_single", uns_spatial_keys))
        if len(library_ids) > 1:
            self.errors.append(
                "uns['spatial'] must contain only two top-level keys: 'is_single' and a library_id. "
                f"More than two top-level keys detected: {library_ids}."
            )
            # Exit if there is more than one library_id as we don't know which library_id to validate.
            return

        # library_id is forbidden if assay is not Visium or is_single is false.
        is_visium_and_uns_is_single = self._is_visium_and_is_single_true()
        if len(library_ids) > 0 and not is_visium_and_uns_is_single:
            self.errors.append(f"uns['spatial'][library_id] {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_FORBIDDEN}.")
            # Exit as library_id is not allowed.
            return

        # Exit if we're not dealing with Visium and _is_single True as no further checks are necessary.
        if not is_visium_and_uns_is_single:
            return

        # library_id is required if assay is Visium and is_single is True.
        if len(library_ids) == 0:
            self.errors.append(
                f"uns['spatial'] must contain at least one key representing the library_id when {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE}."
            )
            # Exit as library_id is missing.
            return

        # Confirm shape of library_id is valid: allowed keys are images and scalefactors.
        library_id_key = library_ids[0]
        uns_library_id = uns_spatial[library_id_key]
        if not isinstance(uns_library_id, dict):
            self.errors.append("uns['spatial'][library_id] must be a dictionary.")
            return
        elif not self._has_no_extra_keys(uns_library_id, ["images", "scalefactors"]):
            self.errors.append(
                "uns['spatial'][library_id] can only contain the keys 'images' and 'scalefactors'."
                f"Detected keys: {list(uns_library_id.keys())}."
            )

        # images is required.
        if "images" not in uns_library_id:
            self.errors.append("uns['spatial'][library_id] must contain the key 'images'.")
        # images is specified: proceed with validation of images.
        elif not isinstance(uns_library_id["images"], dict):
            self.errors.append("uns['spatial'][library_id]['images'] must be a dictionary.")
        else:
            # Confirm shape of images is valid: allowed keys are fullres and hires.
            uns_images = uns_library_id["images"]

            if not self._has_no_extra_keys(uns_images, ["fullres", "hires"]):
                self.errors.append(
                    "uns['spatial'][library_id]['images'] can only contain the keys 'fullres' and 'hires'."
                    f"Detected keys: {list(uns_images.keys())}."
                )

            # hires is required.
            if "hires" not in uns_images:
                self.errors.append("uns['spatial'][library_id]['images'] must contain the key 'hires'.")
            # hires is specified: proceed with validation of hires.
            else:
                self._validate_spatial_image_shape("hires", uns_images["hires"], SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE)

            # fullres is optional.
            uns_fullres = uns_images.get("fullres")
            if uns_fullres is None:
                # Warn if no fullres is specified as it is strongly recommended.
                self.warnings.append(
                    "No uns['spatial'][library_id]['images']['fullres'] was found. "
                    "It is STRONGLY RECOMMENDED that uns['spatial'][library_id]['images']['fullres'] is provided."
                )
            else:
                self._validate_spatial_image_shape("fullres", uns_fullres)

        # scalefactors is required.
        if "scalefactors" not in uns_library_id:
            self.errors.append("uns['spatial'][library_id] must contain the key 'scalefactors'.")
        # scalefactors is specified: proceed with validation of scalefactors.
        elif not isinstance(uns_library_id["scalefactors"], dict):
            self.errors.append("uns['spatial'][library_id]['scalefactors'] must be a dictionary.")
        else:
            # Confirm shape of scalefactors is valid: allowed keys are spot_diameter_fullres and tissue_hires_scalef.
            uns_scalefactors = uns_library_id["scalefactors"]
            if not self._has_no_extra_keys(uns_scalefactors, ["spot_diameter_fullres", "tissue_hires_scalef"]):
                self.errors.append(
                    "uns['spatial'][library_id]['scalefactors'] can only contain the keys "
                    "'spot_diameter_fullres' and 'tissue_hires_scalef'."
                    f"Detected keys: {list(uns_scalefactors.keys())}."
                )

            # spot_diameter_fullres is required.
            if "spot_diameter_fullres" not in uns_scalefactors:
                self.errors.append(
                    "uns['spatial'][library_id]['scalefactors'] must contain the key 'spot_diameter_fullres'."
                )
            # spot_diameter_fullres is specified: proceed with validation.
            else:
                spot_diameter_fullres = uns_scalefactors["spot_diameter_fullres"]
                if not isinstance(spot_diameter_fullres, float):
                    self.errors.append(
                        "uns['spatial'][library_id]['scalefactors']['spot_diameter_fullres'] must be of type float, it is "
                        f"{type(spot_diameter_fullres)}. This must be the value of the spot_diameter_fullres field from scalefactors_json.json"
                    )

            # tissue_hires_scalef is required.
            if "tissue_hires_scalef" not in uns_scalefactors:
                self.errors.append(
                    "uns['spatial'][library_id]['scalefactors'] must contain the key 'tissue_hires_scalef'."
                )
            # tissue_hires_scalef is specified: proceed with validation.
            else:
                tissue_hires_scalef = uns_scalefactors["tissue_hires_scalef"]
                if not isinstance(tissue_hires_scalef, float):
                    self.errors.append(
                        "uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef'] must be of type float, it is "
                        f"{type(tissue_hires_scalef)}. This must be the value of the tissue_hires_scalef field from scalefactors_json.json"
                    )

    def _has_no_extra_keys(self, dictionary: dict, allowed_keys: List[str]) -> bool:
        """
        Determine if the dictionary has only the given allowed keys. Keys can be missing (required
        checks are executed separately) but no additional keys are allowed.

        :param dict dictionary: the dictionary to check.
        :param List[str] allowed_keys: the list of allowed keys.

        :rtype bool
        """
        return set(dictionary.keys()).issubset(allowed_keys)

    def _is_valid_visium_image_shape(self, image: np.ndarray) -> bool:
        """
        Determine if the image has shape (,,3 or 4); image is expected to be a 3D numpy array
        with the size of the last dimension being three or four.

        :param np.ndarray image: the image to check the shape of.

        :return True if image has shape (,,3 or 4), False otherwise.
        :rtype bool
        """
        return len(image.shape) == 3 and image.shape[2] in [3, 4]

    def _is_visium(self) -> bool:
        """
        Determine if the assay_ontology_term_id is Visium (EFO:0010961).

        :return True if assay_ontology_term_id is Visium, False otherwise.
        :rtype bool
        """
        if self.is_visium is None:
            assay_ontology_term_id = self.adata.obs.get("assay_ontology_term_id")
            self.is_visium = assay_ontology_term_id is not None and (assay_ontology_term_id == ASSAY_VISIUM).any()
        return self.is_visium

    def _validate_spatial_image_shape(self, image_name: str, image: np.ndarray, max_dimension: int = None):
        """
        Validate the spatial image is of shape (,,3 or 4) and has a max dimension, if specified. A spatial image
        is either spatial[library_id]['images']['hires'] or spatial[library_id]['images']['fullres']. Errors
        are added to self.errors if any.

        :param str image_name: the name of the image, either "hires" or "fullres".
        :param np.ndarray image: the image to validate.
        :param int max_dimension: the largest allowed dimension of the image, optional.

        :rtype None
        """
        # Image must be an ndarray.
        if not isinstance(image, np.ndarray):
            self.errors.append(
                f"uns['spatial'][library_id]['images']['{image_name}'] must be of numpy.ndarray type, "
                f"it is {type(image)}."
            )
            return

        # Confirm type of ndarray is uint8.
        if image.dtype != np.uint8:
            self.errors.append(
                f"uns['spatial'][library_id]['images']['{image_name}'] must be of type numpy.uint8, "
                f"it is {image.dtype}."
            )

        # Confirm shape of image is valid: allowed shape is (,,3 or 4).
        if not self._is_valid_visium_image_shape(image):
            self.errors.append(
                f"uns['spatial'][library_id]['images']['{image_name}'] must have a length of 3 and "
                "either 3 (RGB color model for example) or 4 (RGBA color model for example) for its last dimension, "
                f"it has shape {image.shape}."
            )

        # Confirm max dimension of image, if specified, is valid.
        if max_dimension is not None and max(image.shape) != max_dimension:
            self.errors.append(
                f"The largest dimension of uns['spatial'][library_id]['images']['{image_name}'] must be "
                f"{max_dimension} pixels, it has a largest dimension of {max(image.shape)} pixels."
            )

    def _deep_check(self):
        """
        Perform a "deep" check of the AnnData object using the schema definition. Adds errors to self.errors if any

        :rtype None
        """
        # Checks DataFrame column name uniqueness
        self._check_var_and_obs_column_name_uniqueness()

        # Checks for deprecated columns
        self._check_deprecated_columns()

        # Checks for invalid columns
        self._check_invalid_columns()

        # Checks that reserved columns are not used
        self._check_column_availability()

        # Checks sparsity
        logger.debug("Validating sparsity...")
        self._validate_sparsity()

        # Checks spatial
        self._check_spatial()

        # Checks Seurat convertibility
        logger.debug("Validating Seurat convertibility...")
        self._validate_seurat_convertibility()

        # Checks each component
        for component_name, component_def in self.schema_def["components"].items():
            logger.debug(f"Validating component: {component_name}")
            component = getattr_anndata(self.adata, component_name)

            # Skip if component does not exist: only useful for adata.raw.var
            if component is None:
                # Check for required components
                if component_def.get("required", False):
                    self.errors.append(f"'{component}' is missing from adata and it's required.")
                continue
            elif component_def["type"] == "dataframe":
                self._validate_dataframe(component_name)
            elif component_def["type"] == "dict":
                self._validate_dict(component, component_name, component_def)
                if component_name == "uns":
                    self._validate_uns_dict(component)
            elif component_def["type"] == "annotation_mapping":
                self._validate_annotation_mapping(component_name, component)
                if component_name == "obsm":
                    self._validate_obsm()
            else:
                raise ValueError(f"Unexpected component type '{component_def['type']}'")

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
        self.reset()

        if h5ad_path:
            logger.debug("Reading the h5ad file...")
            self.adata = read_h5ad(h5ad_path)
            self.h5ad_path = h5ad_path
            self._validate_encoding_version()
            logger.debug("Successfully read the h5ad file")

        # Fetches schema def for latest major schema version
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
    validator = Validator(
        ignore_labels=ignore_labels,
    )
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
