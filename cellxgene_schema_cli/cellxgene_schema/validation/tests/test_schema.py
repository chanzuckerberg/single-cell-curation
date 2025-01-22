import numbers
import re
import warnings
from typing import Dict, List, Mapping, Optional, Tuple, Union

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest
import scipy
from anndata import AnnData
from anndata.compat import DaskArray
from dask.array import map_blocks
from pytest_check import check
from scipy import sparse

from ... import gencode
from ...utils import (
    SPARSE_MATRIX_TYPES,
    count_matrix_nonzero,
    get_matrix_format,
    getattr_anndata,
    is_ontological_descendant_of,
    read_h5ad,
)
from ..constants import (
    ASSAY_SLIDE_SEQV2,
    ASSAY_VISIUM,
    ASSAY_VISIUM_11M,
    ERROR_SUFFIX_IS_SINGLE,
    ERROR_SUFFIX_SPATIAL,
    ERROR_SUFFIX_VISIUM,
    ERROR_SUFFIX_VISIUM_11M,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_FORBIDDEN,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_IN_TISSUE_0,
    ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_REQUIRED,
    GENE_CHECKERS,
    ONTOLOGY_PARSER,
    SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE,
    SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE_VISIUM_11MM,
    VISIUM_11MM_AND_IS_SINGLE_TRUE_MATRIX_SIZE,
    VISIUM_11MM_TISSUE_POSITION_MAX_COL,
    VISIUM_11MM_TISSUE_POSITION_MAX_ROW,
    VISIUM_AND_IS_SINGLE_TRUE_MATRIX_SIZE,
    VISIUM_TISSUE_POSITION_MAX_COL,
    VISIUM_TISSUE_POSITION_MAX_ROW,
)
from ..marks import ignoring_labels


def test_anndata_encoding_version(h5ad_path):
    import h5py

    with h5py.File(h5ad_path, "r") as f:
        encoding_dict = dict(f.attrs)
        encoding_version = encoding_dict.get("encoding-version")
        assert encoding_version == "0.1.0", (
            "The h5ad artifact was generated with an AnnData version different from " "0.8.0."
        )


@pytest.mark.parametrize("component_name", ["obs", "var", "raw.var"])
def test_var_and_obs_column_name_uniqueness(adata, component_name):
    """checks that all column names in the 'var' and 'obs' DataFrames are unique"""
    adata_component: pd.DataFrame = getattr_anndata(adata, component_name)
    if adata_component is None:
        return
    duplicate_columns = adata_component.columns[adata_component.columns.duplicated()].tolist()
    assert not duplicate_columns, (
        f"Duplicate column names detected in 'adata.{component_name}' DataFrame. All "
        f"DataFrame column names must be unique. Duplicate columns: {duplicate_columns}"
    )


# convert Validator._check_deprecated_columns to a pytest similar to the above test
@pytest.mark.parametrize("component_name, component_def", pytest.cxg_schema_def["components"].items())
def test_check_deprecated_columns(adata, component_name, component_def):
    """check for columns or keys that have been deprecated"""
    component = getattr_anndata(adata, component_name)
    if component is None:
        return

    if "deprecated_columns" in component_def:
        deprecated_columns = []
        for column in component_def["deprecated_columns"]:
            if column in component:
                deprecated_columns.append(column)
        assert not deprecated_columns, (
            f"The following deprecated columns are present in '{component_name}': " f"{deprecated_columns}"
        )


@pytest.mark.parametrize("component_name", pytest.cxg_schema_def["components"].keys())
def test_invalid_columns(adata, component_name):
    """check for columns or keys that are not in the schema"""
    component = getattr_anndata(adata, component_name)
    if component is None:
        return

    invalid_columns = []
    for column in component:
        if column.startswith("__"):
            invalid_columns.append(column)
    assert not invalid_columns, (
        f"The following invalid columns are present in '{component_name}'. Fields that start with '__' are reserved. "
        f"Invalid columns: {invalid_columns}"
    )


def forbidden_column_iterator():
    for component_name, component_def in pytest.cxg_schema_def["components"].items():
        if "forbidden_columns" in component_def:
            for forbidden_column in component_def["forbidden_columns"]:
                yield component_name, forbidden_column


@pytest.mark.parametrize("component_name, forbidden_column", forbidden_column_iterator())
def test_for_forbidden_column(adata, component_name, forbidden_column):
    """check for columns that are reserved in components and validate that they are available as expected"""
    component = getattr_anndata(adata, component_name)
    if component is None:
        return
    assert forbidden_column not in component, f"The column '{forbidden_column}' is forbidden in '{component_name}'"


def reserved_column_iterator():
    for component_name, component_def in pytest.cxg_schema_def["components"].items():
        if "reserved_columns" in component_def:
            for reserved_column in component_def["reserved_columns"]:
                yield component_name, reserved_column


@pytest.mark.parametrize("component_name, reserved_column", reserved_column_iterator())
@ignoring_labels
def test_for_reserved_column(adata, component_name, reserved_column):
    """check for columns that are reserved in components and validate that they are available as expected"""
    component = getattr_anndata(adata, component_name)
    if component is None:
        return
    assert (
        reserved_column not in component
    ), f"The column '{reserved_column}' is a reserved name in '{component_name}'. Remove it from h5ad and try again."


def columns_of_column_iterator():
    for component_name, component_def in pytest.cxg_schema_def["components"].items():
        for attr_name in ["columns", "index"]:
            if attr_name in component_def:
                attr_def = component_def[attr_name]
                if "add_labels" in attr_def:
                    for label_def in attr_def["add_labels"]:
                        yield component_name, label_def["to_column"]


@pytest.mark.parametrize("component_name, reserved_name", columns_of_column_iterator())
def test_for_columns_mapped_to_other_columns(adata, component_name, reserved_name):
    """
    This method checks a single reserved column in adata.obs or adata.var and adds a message to error if
        it already exists

    check for columns that map to other columns, for post-upload annotation"""
    component = getattr_anndata(adata, component_name)
    if component is None:
        return
    assert reserved_name not in component, (
        f"Add labels error: Column '{reserved_name}' is a reserved column name "
        f"of '{component}'. Remove it from h5ad and try again."
    )


def _test_sparsity(x, x_name):
    """
    calculates sparsity of x and raw.x, if bigger than indicated in the schema and not a scipy sparse
    matrix, then adds to warnings"""
    max_sparsity = float(pytest.cxg_schema_def["sparsity"])

    # Check sparsity
    matrix_format = get_matrix_format(x)
    if matrix_format == "csr":
        return
    assert matrix_format != "unknown"

    nnz = count_matrix_nonzero(x)
    sparsity = 1 - nnz / np.prod(x.shape)
    if sparsity > max_sparsity:
        warnings.warn(
            f"Sparsity of '{x_name}' is {sparsity} which is greater than {max_sparsity}, "
            f"and it is not a 'scipy.sparse.csr_matrix'. It is STRONGLY RECOMMENDED "
            f"to use this type of matrix for the given sparsity.",
            UserWarning,
            stacklevel=2,
        )


def test_X_sparsity(adata):
    _test_sparsity(adata.X, "adata.X")


def test_raw_X_sparsity(adata):
    if adata.raw:
        _test_sparsity(adata.raw.X, "adata.raw.X")
    else:
        pytest.skip("No raw data found in 'adata.raw'.")


# TODO: test layers independently
def test_layers_sparsity(adata):
    if adata.layers:
        for layer_name, layer in adata.layers.items():
            _test_sparsity(layer, f"adata.layers[{layer_name}]")
    else:
        pytest.skip("No layers found in 'adata.layers'.")


def test_obsm_exists(adata):
    assert adata.obsm, "No embeddings found in 'adata.obsm'."


@pytest.mark.depends(on=["test_obsm_exists"])
def test_obms(adata, is_supported_spatial_assay):
    """
    Validates the embedding dictionary -- it checks that all values of adata.obsm are numpy arrays with the correct
    dimension. Adds errors to errors if any. Checks that the keys start with "X_", have no whitespace, and have
    a suffix at least 1 character long. For keys that don't start with "X_", we will run them through the same
    validation checks, but raise warnings instead of errors.
    """
    obsm_with_x_prefix = 0
    for key, value in adata.obsm.items():

        regex_pattern = r"^[a-zA-Z][a-zA-Z0-9_.-]*$"
        key_is_spatial = key.lower() == "spatial"

        unknown_key = False
        if key.startswith("X_"):
            obsm_with_x_prefix += 1
            if key.lower() == "x_spatial":
                check.fail(f"Embedding key in 'adata.obsm' {key} cannot be used. 'X_spatial' is a reserved key.")
            if not re.match(regex_pattern, key):
                check.fail(
                    f"Suffix for embedding key in 'adata.obsm' {key} does not match the regex pattern {regex_pattern}."
                )
        elif not key_is_spatial:
            if not re.match(regex_pattern, key):
                check.fail(f"Embedding key in 'adata.obsm' {key} does not match the regex pattern {regex_pattern}.")
            warnings.warn(
                f"Embedding key in 'adata.obsm' {key} is not 'spatial' nor does it start with 'X_'. Thus, it will "
                f"not be available in Explorer",
                UserWarning,
                stacklevel=2,
            )
            unknown_key = True

        if not isinstance(value, np.ndarray):
            check.fail(f"All embeddings have to be of 'numpy.ndarray' type, 'adata.obsm['{key}']' is {type(value)}.")
            continue

        if len(value.shape) < 2:
            check.fail(
                f"All embeddings must at least two dimensions. 'adata.obsm['{key}']' has a shape length of '"
                f"{len(value.shape)}'."
            )
        else:
            if value.shape[0] != adata.n_obs:
                check.fail(
                    f"All embeddings must have as many rows as cells. 'adata.obsm['{key}']' has rows='"
                    f"{value.shape[0]}'."
                )

            if unknown_key and value.shape[1] < 1:
                check.fail(
                    f"All unspecified embeddings must have at least one column. 'adata.obsm['{key}']' has columns='"
                    f"{value.shape[1]}'."
                )

            if not unknown_key and value.shape[1] < 2:
                check.fail(
                    f"All 'X_' and 'spatial' embeddings must have at least two columns. 'adata.obsm['{key}']' has "
                    f"columns='{value.shape[1]}'."
                )

        if not (np.issubdtype(value.dtype, np.integer) or np.issubdtype(value.dtype, np.floating)):
            check.fail(
                f"adata.obsm['{key}'] has an invalid data type. It should be "
                "float, integer, or unsigned integer of any precision (8, 16, 32, or 64 bits)."
            )
        else:
            if np.isinf(value).any():
                check.fail(f"adata.obsm['{key}'] contains positive infinity or negative infinity values.")

            if key_is_spatial and np.any(np.isnan(value)):
                check.fail("adata.obs['spatial] contains at least one NaN value.")
            elif np.all(np.isnan(value)):
                check.fail(f"adata.obsm['{key}'] contains all NaN values.")

    if not is_supported_spatial_assay and obsm_with_x_prefix == 0:
        check.fail("At least one embedding in 'obsm' has to have a key with an 'X_' prefix.")


def test_validate_uns(adata):
    uns = getattr_anndata(adata, "uns")
    obs = adata.obs
    # Mapping from obs column name to number of unique categorical values
    category_mapping = {}

    # Check for categorical dtypes in the dataframe directly
    for column_name in obs.columns:
        column = obs[column_name]
        if column.dtype.name == "category":
            category_mapping[column_name] = column.nunique()

    for key, value in uns.items():
        if any(
            isinstance(value, sparse_class)
            for sparse_class in (scipy.sparse.csr_matrix, scipy.sparse.csc_matrix, scipy.sparse.coo_matrix)
        ):
            if value.nnz == 0:  # number non-zero
                check.fail(f"uns['{key}'] cannot be an empty value.")
        elif value is not None and not isinstance(value, (np.bool_, bool, numbers.Number)) and len(value) == 0:
            check.fail(f"uns['{key}'] cannot be an empty value.")
        if key.endswith("_colors"):
            # 1. Verify that the corresponding categorical field exists in obs
            column_name = key.replace("_colors", "")
            obs_unique_values = category_mapping.get(column_name)
            if not obs_unique_values:
                error_message = f"Colors field uns[{key}] does not have a corresponding categorical field in obs"
                if column_name in obs.columns:
                    error_message += f". {column_name} is present but is dtype {obs[column_name].dtype.name}"
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
                check.fail(error_message)
                continue
            # 2. Verify that the value is a numpy array
            if value is None or not isinstance(value, np.ndarray):
                check.fail(f"Colors field uns['{key}'] must be of 'numpy.ndarray' type, it is {type(value)}")
                # Skip over all subsequent validations which expect a numpy array
                continue
            # 3. Verify that we have strings in the array
            all_strings = all(isinstance(color, str) for color in value)
            if not all_strings:
                check.fail(f"Colors in uns[{key}] must be strings. Found: {value} which are {value.dtype.name}")
                continue
            # 4. Verify that we have at least as many colors as unique values in the corresponding categorical field
            if len(value) < obs_unique_values:
                check.fail(
                    f"Annotated categorical field {key.replace('_colors', '')} must have at least {obs_unique_values} "
                    f"color options "
                    f"in uns[{key}]. Found: {value}"
                )
            # 5. Verify that either all colors are hex OR all colors are CSS4 named colors strings
            all_hex_colors = all(re.match(r"^#([0-9a-fA-F]{6})$", color) for color in value)
            all_css4_colors = all(color in mcolors.CSS4_COLORS for color in value)
            if not (all_hex_colors or all_css4_colors):
                check.fail(
                    f"Colors in uns[{key}] must be either all hex colors or all CSS4 named colors. Found: {value}"
                )


def components_iterator():
    for component_name, component_def in pytest.cxg_schema_def["components"].items():
        yield component_name, component_def


@pytest.mark.parametrize("component_name, component_def", components_iterator())
def test_components(adata, component_name, component_def):
    component = getattr_anndata(adata, component_name)

    # Skip if component does not exist: only useful for adata.raw.var
    if component is None:
        # Check for required components
        if component_def.get("required", False):
            check.fail(f"'{component_name}' is missing from adata and is required.")
            return
        pytest.skip(f"'{component_name}' is not present in adata.")
    elif component_def["type"] == "dataframe":
        _validate_dataframe(adata, component_name, component_def)
    elif component_def["type"] == "dict":
        _validate_dict(adata, component, component_name, component_def)
    elif component_def["type"] == "annotation_mapping":
        _validate_annotation_mapping(component_name, component)


def _validate_annotation_mapping(component_name: str, component: Mapping):
    for key, value in component.items():
        # Check for empty ndarrays
        if isinstance(value, np.ndarray) and not value.size:
            check.fail(f"The size of the ndarray stored for a 'adata.{component_name}['{key}']' MUST NOT be zero.")


def _validate_dataframe(adata, df_name, df_def):
    """
    Verifies the dataframe follows the schema. Adds errors to errors if any

    :param str df_name: Name of dataframe in the adata (e.g. "obs")

    :rtype None
    """

    def _format_error_message(msg):
        return f"Column '{column_name}' in dataframe '{df_name}' {msg}"

    df = getattr_anndata(adata, df_name)

    # Validate index if needed
    if "index" in df_def:
        _validate_column(pd.Series(df.index), "index", df_name, df_def["index"])

    low_rows_threshold = df_def.get("warn_if_less_than_rows")
    if low_rows_threshold is not None:
        num_rows = df.shape[0]
        if num_rows < low_rows_threshold:
            warnings.warn(
                f"Dataframe '{df_name}' only has {num_rows} rows. "
                f"Features SHOULD NOT be filtered from expression matrix.",
                UserWarning,
                stacklevel=2,
            )

    for column_name in df.columns:
        column = df[column_name]
        if column.dtype.name != "category":
            # Check for columns with mixed values, which is not supported by anndata
            value_types = {type(x) for x in column.values}
            if len(value_types) != 1:
                check.fail(_format_error_message(f"cannot contain mixed types. Found {value_types}."))
        else:
            # Check for columns that have a category defined 0 times (obs only)
            if df_name == "obs":
                for category in column.dtype.categories:
                    if category not in column.values:
                        warnings.warn(
                            _format_error_message(
                                f"contains a category '{category}' with "
                                f"zero observations. These categories will be removed when `--add-labels` flag is "
                                f"present."
                            ),
                            UserWarning,
                            stacklevel=2,
                        )
            categorical_types = {type(x) for x in column.dtype.categories.values}
            # Check for columns that have illegal categories, which are not supported by anndata
            blocked_categorical_types = {bool}
            illegal_categorical_types = categorical_types & blocked_categorical_types
            if illegal_categorical_types:
                check.fail(_format_error_message(f"contains {illegal_categorical_types=}."))
            # Check for categorical column has mixed types, which is not supported by anndata
            categorical_types = {type(x) for x in column.dtype.categories.values}
            if len(categorical_types) > 1:
                check.fail(
                    _format_error_message(
                        f"contains {len(categorical_types)} categorical types. Only one type is allowed."
                    )
                )

    # Validate columns
    if "columns" in df_def:
        for column_name in df_def["columns"]:
            if column_name not in df.columns:
                check.fail(f"Dataframe '{df_name}' is missing column '{column_name}'.")
                continue

            column_def = df_def["columns"][column_name]
            column = getattr(df, column_name)

            # First check if there are dependencies with other columns and work with a subset of the data if so
            if "dependencies" in column_def:
                column = _validate_column_dependencies(df, df_name, column_name, column_def["dependencies"])

            # If after validating dependencies there's still values in the column, validate them.
            if len(column) > 0:
                if "warning_message" in column_def:
                    warnings.warn(column_def["warning_message"], UserWarning, stacklevel=2)
                _validate_column(column, column_name, df_name, column_def)


def _validate_column(
    column: pd.Series, column_name: str, df_name: str, column_def: dict, default_error_message_suffix=None
):
    """
    Given a schema definition and the column of a dataframe, verify that the column satisfies the schema.
    If there are any errors, it adds them to errors

    :param pandas.Series column: Column of a dataframe to validate
    :param str column_name: Name of the column in the dataframe
    :param str df_name: Name of the dataframe
    :param dict column_def: schema definition for this specific column,
    e.g. schema_def["obs"]["columns"]["cell_type_ontology_term_id"]
    :param str default_error_message_suffix: default error message suffix to be added to errors found here

    :rtype None
    """

    def _format_error_message(msg):
        msg = f"Column '{column_name}' in dataframe '{df_name}' {msg}"
        return f"{msg} {default_error_message_suffix}" if default_error_message_suffix else msg

    if column_def.get("unique") and column.nunique() != len(column):
        check.fail(_format_error_message("is not unique."))

    if column_def.get("type") == "bool" and column.dtype != bool:
        check.fail(_format_error_message(f"must be boolean, not '{column.dtype.name}'."))

    if column_def.get("type") == "categorical":
        if column.dtype.name != "category":
            check.fail(_format_error_message(f"must be categorical, not {column.dtype.name}."))
        else:
            if column_def.get("subtype") == "str":
                if column.dtype.categories.dtype != "object" and column.dtype.categories.dtype != "string":
                    check.fail(
                        _format_error_message(f"must be object or string, not" f" {column.dtype.categories.dtype}.")
                    )
                else:
                    if any(len(cat.strip()) == 0 for cat in column.dtype.categories):
                        check.fail(_format_error_message("must not contain empty values."))

            # check for null values--skip on column defs with enums, since it will already be part of that check
            if not column_def.get("enum") and column.isnull().any():
                check.fail(_format_error_message("must not contain NaN values."))

    if column_def.get("type") == "feature_is_filtered":
        _validate_column_feature_is_filtered(column, column_name, df_name)

    if column_def.get("type") == "genetic_ancestry_value":
        _validate_individual_genetic_ancestry_value(column, column_name)

    if "enum" in column_def:
        bad_enums = [v for v in column.drop_duplicates() if v not in column_def["enum"]]
        if bad_enums:
            check.fail(
                _format_error_message(
                    f"contains invalid values " f"'{bad_enums}'. Values must be one of {column_def['enum']}"
                )
            )

    if column_def.get("type") == "feature_id":
        # Validates each id
        for feature_id in column:
            _validate_feature_id(feature_id, df_name)

    if column_def.get("type") == "curie":
        # Check for NaN values
        if column.isnull().any():
            check.fail(_format_error_message("must not contain NaN values."))
            return

        if "curie_constraints" in column_def:
            for term_str in column.drop_duplicates():
                _validate_curie_str(term_str, column_name, column_def["curie_constraints"])


def _validate_column_feature_is_filtered(column: pd.Series, column_name: str, df_name: str):
    """
    Validates the "is_feature_filtered" in adata.var. This column must be bool, and for genes that are set to
    True, their expression values in X must be 0.
    If there are any errors, it adds them to errors.

    :rtype none
    """

    if column.dtype != bool:
        check.fail(f"Column '{column_name}' in dataframe '{df_name}' must be boolean, not '{column.dtype.name}'.")
        return

    if sum(column) > 0:
        n_nonzero = count_matrix_nonzero(adata.X[:, column])

        if n_nonzero > 0:
            check.fail(
                f"Some features are 'True' in '{column_name}' of dataframe '{df_name}', but there are "
                f"{n_nonzero} non-zero values in the corresponding columns of the matrix 'X'. All values for "
                f"these features must be 0."
            )


def _validate_individual_genetic_ancestry_value(column: pd.Series, column_name: str):
    """
    The following fields are valid for genetic_ancestry_value columns:
    - float values between 0 and 1
    - float('nan')
    """
    if column.dtype != float:
        check.fail(f"Column '{column_name}' in obs must be float, not '{column.dtype.name}'.")
        return

    def is_individual_value_valid(value):
        if isinstance(value, (float, int)) and 0 <= value <= 1:
            return True
        # Ensures only float('nan') or numpy.nan is valid, None is invalid
        if isinstance(value, float) and pd.isna(value):
            return True
        return False

    # Identify invalid values
    invalid_values = column[~column.map(is_individual_value_valid)]

    if not invalid_values.empty:
        check.fail(
            f"Column '{column_name}' in obs contains invalid values: {invalid_values.to_list()}. "
            f"Valid values are floats between 0 and 1 or float('nan')."
        )


def _validate_feature_id(feature_id: str, df_name: str):
    """
    Validates a feature id, i.e. checks that it's present in the reference
    If there are any errors, it adds them to errors and adds it to the list of invalid features

    :param str feature_id: the feature id to be validated
    :param str df_name: name of dataframe the feauter id comes from (var or raw.var)

    :rtype none
    """

    organism = gencode.get_organism_from_feature_id(feature_id)

    if not organism:
        check.fail(
            f"Could not infer organism from feature ID '{feature_id}' in '{df_name}', " f"make sure it is a valid ID."
        )
        return
    gene_checkers = GENE_CHECKERS
    if organism not in gene_checkers:
        gene_checkers[organism] = gencode.GeneChecker(organism)

    if not gene_checkers[organism].is_valid_id(feature_id):
        check.fail(f"'{feature_id}' is not a valid feature ID in '{df_name}'.")


def _validate_curie_str(term_str: str, column_name, curie_constraints: dict) -> None:
    """
    Validate a curie str based on some constraints. If there are any errors, it adds them to errors

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
        check.fail(f"'{term_str}' in '{column_name}' is not a valid value of '{column_name}'.")
        return

    if not isinstance(term_str, str):
        check.fail(f"'{term_str}' in '{column_name}' is not a valid ontology term value, it must be a string.")
        return

    # if multi_term is defined, split str into individual ontology terms and validate each
    if "multi_term" in curie_constraints:
        delimiter = curie_constraints["multi_term"]["delimiter"]
        term_ids = term_str.split(delimiter)
        for term_id in term_ids:
            _validate_curie(term_id, column_name, curie_constraints)
        if (
            curie_constraints["multi_term"].get("sorted", False)
            and len(term_ids) > 1
            and not all(term_ids[i].strip() <= term_ids[i + 1].strip() for i in range(len(term_ids) - 1))
        ):
            check.fail(f"'{term_str}' in '{column_name}' is not in ascending lexical order.")
        if len(set(term_ids)) != len(term_ids):
            check.fail(f"'{term_str}' in '{column_name}' contains duplicates.")
    else:
        _validate_curie(term_str, column_name, curie_constraints)


def _validate_curie(term_id: str, column_name: str, curie_constraints: dict):
    """
    Validate a single curie term id based on some constraints.
    If there are any errors, it adds them to errors

    :param str term_id: the curie term id to validate
    :param str column_name: Name of the column in the dataframe
    :param dict curie_constraints: constraints for the curie term to be validated,
    this part of the schema definition

    :rtype None
    """
    # If the term id does not belong to an allowed ontology, the subsequent checks are redundant
    if not _validate_curie_ontology(term_id, column_name, curie_constraints["ontologies"]):
        return

    # Check if term_id is forbidden by schema definition despite being a valid ontology term
    if "forbidden" in curie_constraints:
        if "terms" in curie_constraints["forbidden"] and term_id in curie_constraints["forbidden"]["terms"]:
            check.fail(f"'{term_id}' in '{column_name}' is not allowed.")
            return
        if "ancestors" in curie_constraints["forbidden"] and _has_forbidden_curie_ancestor(
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
            and _validate_curie_ancestors(term_id, curie_constraints["allowed"]["ancestors"])
        ):
            is_allowed = True

        if not is_allowed:
            check.fail(f"'{term_id}' in '{column_name}' is not an allowed term id.")


def _validate_curie_ontology(term_id: str, column_name: str, allowed_ontologies: List[str]) -> bool:
    """
    Validate a single curie term id belongs to specified ontologies. If it does belong to an allowed ontology
    verifies that it is not deprecated (obsolete).
    If there are any errors, it adds them to errors

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
            check.fail(f"'{term_id}' in '{column_name}' is a deprecated term id of '{ontology_name}'.")
            return False

    if sum(checks) == 0:
        check.fail(
            f"'{term_id}' in '{column_name}' is not a valid ontology term id of '{', '.join(allowed_ontologies)}'."
        )
        return False
    return True


def _validate_curie_ancestors(
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
                is_descendant = ancestor in ONTOLOGY_PARSER.get_term_ancestors(term_id, inclusive)
                checks.append(is_descendant)

    if True not in checks:
        return False
    return True


def _validate_column_dependencies(
    df: pd.DataFrame, df_name: str, column_name: str, dependencies: List[dict]
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
        terms_to_match = set()
        column_to_match = dependency_def["rule"]["column"]
        if "match_ancestors_inclusive" in dependency_def["rule"]:
            ancestors = dependency_def["rule"]["match_ancestors_inclusive"]["ancestors"]
            for ancestor in ancestors:
                terms_to_match.update(ONTOLOGY_PARSER.get_term_descendants(ancestor, include_self=True))
        if "match_exact" in dependency_def["rule"]:
            terms_to_match.update(dependency_def["rule"]["match_exact"]["terms"])
        try:
            match_query = df[column_to_match].isin(terms_to_match)
            match_df = df[match_query]
            column = getattr(match_df, column_name)
            error_message_suffix = dependency_def.get("error_message_suffix", None)
            if not error_message_suffix:
                matched_values = list(getattr(match_df, column_to_match).unique())
                error_message_suffix = f"when '{column_to_match}' is in {matched_values}"
        except KeyError:
            check.fail(
                f"Checking values with dependencies failed for adata.{df_name}['{column_name}'], "
                f"this is likely due to missing dependent column in adata.{df_name}."
            )
            return pd.Series(dtype=np.float64)

        all_rules.append(match_query)
        _validate_column(column, column_name, df_name, dependency_def, error_message_suffix)

    # Return column of data that was not matched by any of the rules
    column = getattr(df[~np.logical_or.reduce(all_rules)], column_name)

    return column


def _validate_dict(adata, dictionary: dict, dict_name: str, dict_def: dict):
    """
    Verifies the dictionary follows the schema. Adds errors to errors if any

    :param Anndata adata: The adata object
    :param str dictionary: The dictionary to validate
    :param str dict_name: Name of dictionary in the adata (e.g. "uns")
    :param str dict_def: The schema definition for this specific dictionary

    :rtype None
    """

    for key, value_def in dict_def["keys"].items():
        if key not in dictionary:
            if value_def.get("required", False):
                check.fail(f"'{key}' in '{dict_name}' is not present.")
            continue

        value = dictionary[key]

        if value_def["type"] == "string":
            if not _validate_str_in_dict(value, dict_name, key):
                continue

            if "enum" in value_def:
                _validate_enum_in_dict(value, value_def["enum"], dict_name, key)

        if value_def["type"] == "match_obsm_keys":
            if not _validate_str_in_dict(value, dict_name, key):
                continue

            if value not in adata.obsm:
                check.fail(f"'{value}' in '{dict_name}['{key}']' is not valid, " f"it must be a key of 'adata.obsm'.")

        if value_def["type"] == "list":
            if not (isinstance(value, (list, np.ndarray))):
                check.fail(f"'{value}' in '{dict_name}['{key}']' is not valid, " f"it must be a list or numpy array.")
                continue

            _validate_list(key, value, value_def["element_type"])


def _validate_str_in_dict(value, dict_name: str, key: str):
    """
    Validates that a value from a dictionary is a string and it does not have leading, trailing or double spaces.

    :param str value: The dictionary to validate
    :param str dict_name: Name of dictionary in the adata (e.g. "uns")
    :param str key: The key in the dictionary

    :rtype bool
    :return True if passed, False otherwise
    """

    is_valid = True

    if not isinstance(value, str):
        check.fail(f"'{value}' in '{dict_name}['{key}']' is not valid, it must be a string.")

        is_valid = False
    else:
        if value != value.rstrip():
            check.fail(f"'{value}' in '{dict_name}['{key}']' is not valid, it contains trailing spaces.")

            is_valid = False

        if value != value.lstrip():
            check.fail(f"'{value}' in '{dict_name}['{key}']' is not valid, it contains leading spaces.")

            is_valid = False

        if value.strip() != " ".join(value.split()):
            check.fail(f"'{value}' in '{dict_name}['{key}']' is not valid, it contains double spaces.")

            is_valid = False

    return is_valid


def _validate_list(adata: AnnData, list_name: str, current_list: List[str], element_type: str):
    """
    Validates the elements of a list based on the type definition.

    :param AnnData adata: The adata object
    :param str list_name: name of list to use for error messages (if any)
    :param str current_list: the list to be validated
    :param str element_type: type to be validated

    :rtype None
    """

    for i in current_list:
        if element_type == "match_obs_columns" and i not in adata.obs.columns:
            check.fail(
                f"Value '{i}' of list '{list_name}' is not a column in 'adata.obs'."
            )  # TODO this should not assume
            # obs.


def _validate_enum_in_dict(value, enum: List[str], dict_name: str, key: str):
    """
    Validates that a value from a dictionary is part of a list.

    :param str value: The dictionary to validate
    :param List[str] enum: The allowed values
    :param str dict_name: Name of dictionary in the adata (e.g. "uns")
    :param str key: The key in the dictionary

    :rtype  None
    """

    if value not in enum:
        check.fail(f"'{value}' in '{dict_name}['{key}']' is not valid. " f"Allowed terms: {enum}.")


def _has_forbidden_curie_ancestor(term_id: str, column_name: str, forbidden_def: Dict[str, List[str]]) -> bool:
    """
    Validate if a single curie term id is a descendant term of any forbidden ancestors.
    If there is a forbidden ancestor detected, it adds it to errors.

    :param str term_id: the curie term id to validate
    :param str column_name: original column name in adata where the term_id comes from (used for error messages)
    :param Dict[str, List[str] forbidden_def: mapping of ontologies to list of ancestor terms to validate against

    :returns bool
    """
    for ontology_name in forbidden_def:
        for ancestor in forbidden_def[ontology_name]:
            if ancestor in ONTOLOGY_PARSER.get_term_ancestors(term_id):
                check.fail(
                    f"'{term_id}' in '{column_name}' is not allowed. Descendant terms of "
                    f"'{ancestor}' are not allowed."
                )
                return True
    return False


def _are_descendants_of(adata, component: str, column: str, ancestors: List[str]) -> bool:
    """
    Checks if elements in the specified column of the component (e.g. 'assay_ontology_term_id' of 'adata.obs') are
    descendants of the given ancestors.

    Ancestors checks are inclusive, meaning that a value is its own ancestor as well.

    :param str component: the name of the component that's been checked.
    :param str column: Column in the component to check
    :param List[str] ancestors: List of ancestors

    :rtype bool
    :return True if any value in column is a descendant of any ancestor.
    """

    curies = getattr(getattr(adata, component), column)
    curies = curies.drop_duplicates()

    for curie in curies:
        if ONTOLOGY_PARSER.is_valid_term_id(curie):
            curie_ancestors = ONTOLOGY_PARSER.get_term_ancestors(curie, include_self=True)
            if bool(set(curie_ancestors) & set(ancestors)):
                return True

    return False


def test_validate_x_raw_x_dimensions(adata, raw_x_loc):
    """
    Validates that X and raw.X have the same shape, if raw.X exists. Adds errors to errors if any.
    """

    if raw_x_loc == "raw.X":
        if adata.n_vars != adata.raw.n_vars:
            check.fail(f"Number of genes in X ({adata.n_vars}) is different " f"than raw.X ({adata.raw.n_vars}).")
        else:
            if not (adata.var.index == adata.raw.var.index).all():
                check.fail("Index of 'raw.var' is not identical to index of 'var'.")
        if adata.n_obs != adata.raw.n_obs:
            check.fail(f"Number of cells in X ({adata.n_obs}) is different " f"than raw.X ({adata.raw.n_obs}).")
        else:
            if not (adata.obs_names == adata.raw.obs_names).all():
                check.fail("Cells in X and raw.X are different.")


def test_raw_check(adata):
    # Asses if we actually need to perform validation of raw based on the rules in the schema
    # As of now, this means that we only do validation of raw if it's RNA data
    checks = []
    for component, component_rules in pytest.cxg_schema_def["raw"].items():
        for column, column_rules in component_rules.items():
            for rule, rule_def in column_rules.items():
                if rule == "not_descendants_of":
                    for _, ancestors in rule_def.items():
                        # TODO: should ontology_name be used or removed from the schema? It's not currently used.
                        checks.append(not _are_descendants_of(adata, component, column, ancestors))
                else:
                    raise ValueError(f"'{rule}' rule in raw definition of the schema is not implemented ")


@pytest.mark.depends(on=["test_validate_x_raw_x_dimensions", "test_raw_check"])
def test_validate_raw(has_valid_raw, raw_x_loc):
    """
    Validates raw only if the rules in the schema definition are fulfilled and that X and raw.X have the same shape
    The validation entails checking that:
     1. X and raw.X have the same column and row indices
     2. there's an expression matrix containing raw (integer) values, first in adata.raw.X and then adata.X if
     the former does not exist.
     3. For applicable assays, checks that each row has at least one non-zero value

    :rtype None
    """

    # If both "raw.X" and "X" exist but neither are raw
    # This is testing for when sometimes data contributors put a normalized matrix in both "X" and "raw.X".
    if not has_valid_raw and raw_x_loc == "raw.X":
        check.fail("Raw data may be missing: data in 'raw.X' does not meet schema requirements.")

    # Only "X" exists but it's not raw
    # This is testing for when there is only a normalized matrix in "X" and there is no "raw.X".
    if not has_valid_raw and raw_x_loc == "X":
        check.fail("Raw data is missing: there is only a normalized matrix in X and no raw.X")

    # If raw data is in X and there is nothing in raw.X (i.e. normalized values are not provided), then
    # add a warning because normalized data for RNA data is STRONGLY RECOMMENDED
    if has_valid_raw and raw_x_loc == "X":
        warnings.warn(
            "Only raw data was found, i.e. there is no 'raw.X'. "
            "It is STRONGLY RECOMMENDED that 'final' (normalized) data is provided.",
            UserWarning,
            stacklevel=2,
        )


@pytest.fixture(scope="session")
def raw_x_loc(adata):
    """
    gets raw x location (best guess, i.e. not guarantee it's actually raw)
    """
    if adata.raw:
        return "raw.X"
    else:
        return "X"


@pytest.fixture(scope="session")
def raw_x(adata):
    """
    gets raw x (best guess, i.e. not guarantee it's actually raw)
    """

    if adata.raw:
        return adata.raw.X
    else:
        return adata.X


@pytest.fixture(scope="session")
def visium_and_is_single_true_matrix_size_with_error_suffix(
    adata, is_visium_including_descendants
) -> Tuple[Optional[int], Optional[str]]:
    """
    Returns the required matrix size based on assay type, if applicable, else returns None.
    """
    visium_and_is_single_true_matrix_size = None
    visium_error_suffix = None
    # Visium 11M's raw matrix size is distinct from other visium assays
    if bool(
        adata.obs["assay_ontology_term_id"]
        .apply(lambda t: is_ontological_descendant_of(ONTOLOGY_PARSER, t, ASSAY_VISIUM_11M, True))
        .astype(bool)
        .any()
    ):
        visium_error_suffix = f"{ERROR_SUFFIX_VISIUM_11M} and {ERROR_SUFFIX_IS_SINGLE}"
        visium_and_is_single_true_matrix_size = VISIUM_11MM_AND_IS_SINGLE_TRUE_MATRIX_SIZE
    elif is_visium_including_descendants:
        visium_error_suffix = f"{ERROR_SUFFIX_VISIUM} and {ERROR_SUFFIX_IS_SINGLE}"
        visium_and_is_single_true_matrix_size = VISIUM_AND_IS_SINGLE_TRUE_MATRIX_SIZE
    return visium_and_is_single_true_matrix_size, visium_error_suffix


@pytest.fixture(scope="session")
def hires_max_dimension_size(adata, is_visium_including_descendants) -> Optional[int]:
    """
    Returns the restricted hires image dimension based on assay type, if applicable, else returns None.
    """
    _hires_max_dimension_size = None
    # Visium 11M's max dimension size is distinct from other visium assays
    if bool(
        adata.obs["assay_ontology_term_id"]
        .apply(lambda t: is_ontological_descendant_of(ONTOLOGY_PARSER, t, ASSAY_VISIUM_11M, True))
        .astype(bool)
        .any()
    ):
        _visium_error_suffix = ERROR_SUFFIX_VISIUM_11M
        _hires_max_dimension_size = SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE_VISIUM_11MM
    elif is_visium_including_descendants:
        _visium_error_suffix = ERROR_SUFFIX_VISIUM
        _hires_max_dimension_size = SPATIAL_HIRES_IMAGE_MAX_DIMENSION_SIZE
    return _hires_max_dimension_size


def _is_valid_visium_image_shape(image: np.ndarray) -> bool:
    """
    Determine if the image has shape (,,3 or 4); image is expected to be a 3D numpy array
    with the size of the last dimension being three or four.

    :param np.ndarray image: the image to check the shape of.

    :return True if image has shape (,,3 or 4), False otherwise.
    :rtype bool
    """
    return len(image.shape) == 3 and image.shape[2] in [3, 4]


def has_no_extra_keys(dictionary: dict, allowed_keys: List[str]) -> bool:
    """
    Determine if the dictionary has only the given allowed keys. Keys can be missing (required
    checks are executed separately) but no additional keys are allowed.

    :param dict dictionary: the dictionary to check.
    :param List[str] allowed_keys: the list of allowed keys.

    :rtype bool
    """
    return set(dictionary.keys()).issubset(allowed_keys)


def test_spatial_embeddings(adata, is_single):
    has_spatial_embedding = "spatial" in adata.obsm
    if is_single and not has_spatial_embedding:
        check.fail("'spatial' embedding is required in 'adata.obsm' if " "adata.uns['spatial']['is_single'] is True.")
    elif is_single is None and has_spatial_embedding:
        check.fail(
            "'spatial' embedding is forbidden in 'adata.obsm' if " "adata.uns['spatial']['is_single'] is not set."
        )


@pytest.mark.depends(on=[""])
def test_validate_spatial_assay_ontology_term_id(adata, is_supported_spatial_assay):
    """
    If assay is spatial, all assay ontology term ids should be identical.

    :rtype none
    """
    # Identical check requires assay_ontology_term_id.
    obs = getattr_anndata(adata, "obs")
    if obs is None or "assay_ontology_term_id" not in obs:
        return

    # Identical check is only applicable to spatial datasets.
    if not is_supported_spatial_assay:
        return

    # Validate assay ontology term ids are identical.
    term_count = obs["assay_ontology_term_id"].nunique()
    if term_count > 1:
        check.fail(f"When {ERROR_SUFFIX_SPATIAL}" ", all observations must contain the same value.")


def test_validate_spatial_cell_type_ontology_term_id(
    adata, is_supported_spatial_assay, is_single, is_visium_and_is_single_true
):
    """
    Validate cell type ontology term id is "unknown" if Visium, is_single is True and in_tissue is 0.

    :rtype none
    """
    # skip checks if not a valid spatial assay with a corresponding "in_tissue" column
    if not is_visium_and_is_single_true:
        pytest.skip("not a valid spatial assay")
    elif is_visium_and_is_single_true and "in_tissue" not in adata.obs.columns:
        pytest.skip("valid spatial assay, but missing 'in_tissue' column")

    # Validate all out of tissue (in_tissue==0) spatial spots have unknown cell ontology term
    is_spatial = (
        adata.obs["assay_ontology_term_id"]
        .apply(lambda assay: is_ontological_descendant_of(ONTOLOGY_PARSER, assay, ASSAY_VISIUM, True))
        .astype(bool)
    )
    is_not_tissue = adata.obs["in_tissue"] == 0
    is_not_unknown = adata.obs["cell_type_ontology_term_id"] != "unknown"
    if (is_spatial & is_not_tissue & is_not_unknown).any():
        check.fail(
            f"obs['cell_type_ontology_term_id'] must be 'unknown' when "
            f"{ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_IN_TISSUE_0}."
        )


def test_validate_spatial_is_primary_data(adata, is_single):
    """
    Validate is_primary_data for spatial datasets.
    """
    obs = getattr_anndata(adata, "obs")
    if obs is None or "is_primary_data" not in obs:
        return
    if is_single is False and obs["is_primary_data"].any():
        check.fail("When uns['spatial']['is_single'] is False, obs['is_primary_data'] must be False for all rows.")


@pytest.fixture(scope="session")
def tissue_position_maxes(adata, is_visium_and_is_single_true) -> Optional[Tuple[int, int]]:
    if is_visium_and_is_single_true:
        # visium 11 has different requirements than other visium
        if (
            adata.obs["assay_ontology_term_id"]
            .apply(lambda t: is_ontological_descendant_of(ONTOLOGY_PARSER, t, ASSAY_VISIUM_11M, True))
            .astype(bool)
            .any()
        ):
            return (
                VISIUM_11MM_TISSUE_POSITION_MAX_ROW,
                VISIUM_11MM_TISSUE_POSITION_MAX_COL,
            )
        else:
            return VISIUM_TISSUE_POSITION_MAX_ROW, VISIUM_TISSUE_POSITION_MAX_COL
    pytest.skip("not a valid spatial assay")


def test_validate_spatial_tissue_positions(adata, tissue_position_maxes, is_visium_and_is_single_true):
    """
    Validate tissue positions of spatial datasets.
    """

    def validate_spatial_tissue_position(tissue_position_name: str, min: int, max: int):
        """
        Validate tissue position is allowed and required, and are integers within the given range. Validation is not
        defined in
        schema definition yaml.

        :rtype none
        """
        # Tissue position is foribidden if assay is not Visium and is_single is True.
        if tissue_position_name in adata.obs and (
            not is_visium_and_is_single_true
            or (
                ~(
                    adata.obs["assay_ontology_term_id"]
                    .apply(lambda t: is_ontological_descendant_of(ONTOLOGY_PARSER, t, ASSAY_VISIUM, True))
                    .astype(bool)
                )
                & (adata.obs[tissue_position_name].notnull())
            ).any()
        ):
            check.fail(f"obs['{tissue_position_name}'] {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_FORBIDDEN}.")
            return

        # Exit if we're not dealing with Visium and _is_single True as no further checks are necessary.
        if not is_visium_and_is_single_true:
            return

        # At this point, is_single is True and:
        # - there's at least one row with Visum, tissue position column is required
        # - for any Visium row, tissue position is required.
        if (
            tissue_position_name not in adata.obs
            or (
                (
                    adata.obs["assay_ontology_term_id"]
                    .apply(lambda t: is_ontological_descendant_of(ONTOLOGY_PARSER, t, ASSAY_VISIUM, True))
                    .astype(bool)
                )
                & (adata.obs[tissue_position_name].isnull())
            ).any()
        ):
            check.fail(f"obs['{tissue_position_name}'] {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_REQUIRED}.")
            return

        # Tissue position must be an int.
        obs_tissue_position = adata.obs.get(tissue_position_name)
        if not np.issubdtype(obs_tissue_position.dtype, np.integer):
            check.fail(f"obs['{tissue_position_name}'] must be of int type, it is {obs_tissue_position.dtype}.")
            return

        # Tissue position must be within the given range.
        if not ((obs_tissue_position >= min) & (obs_tissue_position <= max)).all():
            if tissue_position_name == "in_tissue":
                error_message_token = f"{min} or {max}"
            else:
                error_message_token = f"between {min} and {max}"
            check.fail(
                f"obs['{tissue_position_name}'] must be {error_message_token}, the min and max are "
                f"{obs_tissue_position.min()} and {obs_tissue_position.max()}. "
                f"This must be the value of the column tissue_positions_in_tissue from the tissue_positions_list.csv "
                f"or tissue_positions.csv."
            )

    validate_spatial_tissue_position("array_col", 0, tissue_position_maxes[1])
    validate_spatial_tissue_position("array_row", 0, tissue_position_maxes[0])
    validate_spatial_tissue_position("in_tissue", 0, 1)


def test_spatial_uns(adata, is_supported_spatial_assay, is_visium_and_is_single_true):
    """
    Validate uns spatial-related values of the AnnData object. Validation is not defined in schema definition yaml.
    Errors are added to errors.

    :rtype none
    """

    def validate_spatial_image_shape(image_name: str, image: np.ndarray, max_dimension: int = None):
        """
        Validate the spatial image is of shape (,,3 or 4) and has a max dimension, if specified. A spatial image
        is either spatial[library_id]['images']['hires'] or spatial[library_id]['images']['fullres']. Errors
        are added to errors if any.

        :param str image_name: the name of the image, either "hires" or "fullres".
        :param np.ndarray image: the image to validate.
        :param int max_dimension: the largest allowed dimension of the image, optional.

        :rtype None
        """
        # Image must be an ndarray.
        if not isinstance(image, np.ndarray):
            check.fail(
                f"uns['spatial'][library_id]['images']['{image_name}'] must be of numpy.ndarray type, "
                f"it is {type(image)}."
            )
            return

        # Confirm type of ndarray is uint8.
        if image.dtype != np.uint8:
            check.fail(
                f"uns['spatial'][library_id]['images']['{image_name}'] must be of type numpy.uint8, "
                f"it is {image.dtype}."
            )

        # Confirm shape of image is valid: allowed shape is (,,3 or 4).
        if not _is_valid_visium_image_shape(image):
            check.fail(
                f"uns['spatial'][library_id]['images']['{image_name}'] must have a length of 3 and "
                "either 3 (RGB color model for example) or 4 (RGBA color model for example) for its last dimension, "
                f"it has shape {image.shape}."
            )

        # Confirm max dimension of image, if specified, is valid.
        if max_dimension is not None and max(image.shape) != max_dimension:
            check.fail(
                f"The largest dimension of uns['spatial'][library_id]['images']['{image_name}'] must be "
                f"{max_dimension} pixels, it has a largest dimension of {max(image.shape)} pixels."
            )

    # Exit if uns is not specified. Error is reported in core validate functionality.
    uns_component = getattr_anndata(adata, "uns")
    if uns_component is None:
        pytest.skip("uns is not specified")

    # uns spatial validation is dependent on obs.assay_ontology_term_id; exit if not specified. Error is
    # reported in core validate functionality.
    obs_component = getattr_anndata(adata, "obs")
    if obs_component is None or "assay_ontology_term_id" not in obs_component:
        pytest.skip("obs['assay_ontology_term_id'] is not specified")

    # spatial is forbidden if assay it not a supported spatial assay.
    uns_spatial = adata.uns.get("spatial")
    if uns_spatial is not None and not is_supported_spatial_assay:
        check.fail(f"uns['spatial'] is only allowed when {ERROR_SUFFIX_SPATIAL}")
        return

    # Exit if we aren't dealing with a supported spatial assay as no further checks are necessary.
    if not is_supported_spatial_assay:
        return

    # spatial is required for supported spatial assays.
    if not isinstance(uns_spatial, dict):
        check.fail("A dict in uns['spatial'] is required when " f"{ERROR_SUFFIX_SPATIAL}.")
        return

    # is_single is required.
    if "is_single" not in uns_spatial:
        check.fail("uns['spatial'] must contain the key 'is_single'.")
        # Exit if is_single is missing as all further checks are dependent on its value.
        return

    # is_single must be a boolean.
    uns_is_single = uns_spatial["is_single"]
    if not isinstance(uns_is_single, (np.bool_, bool)):
        check.fail(f"uns['spatial']['is_single'] must be of boolean type, it is {type(uns_is_single)}.")
        # Exit if is_single is not valid as all further checks are dependent on its value.
        return

    # Check there is at most one library_id.
    uns_spatial_keys = list(uns_spatial.keys())
    library_ids = list(filter(lambda x: x != "is_single", uns_spatial_keys))
    if len(library_ids) > 1:
        check.fail(
            "uns['spatial'] must contain only two top-level keys: 'is_single' and a library_id. "
            f"More than two top-level keys detected: {library_ids}."
        )
        # Exit if there is more than one library_id as we don't know which library_id to validate.
        return

    # library_id is forbidden if assay is not Visium or is_single is false.
    if len(library_ids) > 0 and not is_visium_and_is_single_true:
        check.fail(f"uns['spatial'][library_id] {ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE_FORBIDDEN}.")
        # Exit as library_id is not allowed.
        return

    # Exit if we're not dealing with Visium and _is_single True as no further checks are necessary.
    if not is_visium_and_is_single_true:
        return

    # library_id is required if assay is Visium and is_single is True.
    if len(library_ids) == 0:
        check.fail(
            f"uns['spatial'] must contain at least one key representing the library_id when "
            f"{ERROR_SUFFIX_VISIUM_AND_IS_SINGLE_TRUE}."
        )
        # Exit as library_id is missing.
        return

    # Confirm shape of library_id is valid: allowed keys are images and scalefactors.
    library_id_key = library_ids[0]
    uns_library_id = uns_spatial[library_id_key]
    if not isinstance(uns_library_id, dict):
        check.fail("uns['spatial'][library_id] must be a dictionary.")
        return
    elif not has_no_extra_keys(uns_library_id, ["images", "scalefactors"]):
        check.fail(
            "uns['spatial'][library_id] can only contain the keys 'images' and 'scalefactors'."
            f"Detected keys: {list(uns_library_id.keys())}."
        )

    # images is required.
    if "images" not in uns_library_id:
        check.fail("uns['spatial'][library_id] must contain the key 'images'.")
    # images is specified: proceed with validation of images.
    elif not isinstance(uns_library_id["images"], dict):
        check.fail("uns['spatial'][library_id]['images'] must be a dictionary.")
    else:
        # Confirm shape of images is valid: allowed keys are fullres and hires.
        uns_images = uns_library_id["images"]

        if not has_no_extra_keys(uns_images, ["fullres", "hires"]):
            check.fail(
                "uns['spatial'][library_id]['images'] can only contain the keys 'fullres' and 'hires'."
                f"Detected keys: {list(uns_images.keys())}."
            )

        # hires is required.
        if "hires" not in uns_images:
            check.fail("uns['spatial'][library_id]['images'] must contain the key 'hires'.")
        # hires is specified: proceed with validation of hires.
        else:
            _max_size = hires_max_dimension_size
            validate_spatial_image_shape("hires", uns_images["hires"], _max_size)

        # fullres is optional.
        uns_fullres = uns_images.get("fullres")
        if uns_fullres is None:
            # Warn if no fullres is specified as it is strongly recommended.
            warnings.append(
                "No uns['spatial'][library_id]['images']['fullres'] was found. "
                "It is STRONGLY RECOMMENDED that uns['spatial'][library_id]['images']['fullres'] is provided."
            )
        else:
            validate_spatial_image_shape("fullres", uns_fullres)

    # scalefactors is required.
    if "scalefactors" not in uns_library_id:
        check.fail("uns['spatial'][library_id] must contain the key 'scalefactors'.")
    # scalefactors is specified: proceed with validation of scalefactors.
    elif not isinstance(uns_library_id["scalefactors"], dict):
        check.fail("uns['spatial'][library_id]['scalefactors'] must be a dictionary.")
    else:
        # Confirm shape of scalefactors is valid: allowed keys are spot_diameter_fullres and tissue_hires_scalef.
        uns_scalefactors = uns_library_id["scalefactors"]
        if not has_no_extra_keys(uns_scalefactors, ["spot_diameter_fullres", "tissue_hires_scalef"]):
            check.fail(
                "uns['spatial'][library_id]['scalefactors'] can only contain the keys "
                "'spot_diameter_fullres' and 'tissue_hires_scalef'."
                f"Detected keys: {list(uns_scalefactors.keys())}."
            )

        # spot_diameter_fullres is required.
        if "spot_diameter_fullres" not in uns_scalefactors:
            check.fail("uns['spatial'][library_id]['scalefactors'] must contain the key 'spot_diameter_fullres'.")
        # spot_diameter_fullres is specified: proceed with validation.
        else:
            spot_diameter_fullres = uns_scalefactors["spot_diameter_fullres"]
            if not isinstance(spot_diameter_fullres, float):
                check.fail(
                    "uns['spatial'][library_id]['scalefactors']['spot_diameter_fullres'] must be of type float, it is "
                    f"{type(spot_diameter_fullres)}. This must be the value of the spot_diameter_fullres field from "
                    f"scalefactors_json.json"
                )

        # tissue_hires_scalef is required.
        if "tissue_hires_scalef" not in uns_scalefactors:
            check.fail("uns['spatial'][library_id]['scalefactors'] must contain the key 'tissue_hires_scalef'.")
        # tissue_hires_scalef is specified: proceed with validation.
        else:
            tissue_hires_scalef = uns_scalefactors["tissue_hires_scalef"]
            if not isinstance(tissue_hires_scalef, float):
                check.fail(
                    "uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef'] must be of type float, it is "
                    f"{type(tissue_hires_scalef)}. This must be the value of the tissue_hires_scalef field from "
                    f"scalefactors_json.json"
                )


@pytest.fixture(scope="session")
def h5ad_path(request):
    return request.config.getoption("--dataset")


@pytest.fixture(scope="session")
def adata(h5ad_path: str) -> AnnData:
    return read_h5ad(h5ad_path)


@pytest.fixture(scope="session")
def is_single(adata):
    """
    Determine value of uns.spatial.is_single. None if non-spatial.

    :return Value of uns.spatial.is_single if specified, None otherwise.
    :rtype bool | None
    """
    return (
        adata.uns["spatial"]["is_single"]
        if hasattr(adata, "uns")
        and "spatial" in adata.uns
        and isinstance(adata.uns["spatial"], dict)
        and "is_single" in adata.uns["spatial"]
        else None
    )


@pytest.fixture(scope="session")
def is_supported_spatial_assay(adata, is_visium_including_descendants) -> bool:
    """
    Determine if the assay_ontology_term_id is either Visium (EFO:0010961) or Slide-seqV2 (EFO:0030062).

    :return True if assay_ontology_term_id is Visium or Slide-seqV2, False otherwise.
    :rtype bool
    """
    try:
        _spatial = (
            is_visium_including_descendants
            or adata.obs.assay_ontology_term_id.isin([ASSAY_SLIDE_SEQV2]).astype(bool).any()
        )
        is_spatial = bool(_spatial)
    except AttributeError:
        # specific error reporting will occur downstream in the validation
        is_spatial = False
    return is_spatial


@pytest.fixture(scope="session")
def is_visium_and_is_single_true(is_single, is_visium) -> bool:
    """
    Determine if the assay_ontology_term_id is Visium (EFO:0010961) and uns.spatial.is_single is True.

    :return True if assay_ontology_term_id is Visium and is_single_cell is True, False otherwise.
    :rtype bool
    """
    return bool(is_visium and is_single)


@pytest.fixture(scope="session")
def has_valid_raw(
    adata, raw_x, is_visium_and_is_single_true, visium_and_is_single_true_matrix_size_with_error_suffix
) -> bool:
    """
    Checks if the non-zero values for the raw matrix (adata.X or adata.raw.X)
    are positive integers stored as numpy.float32. Also validates that every row contains at least one non-zero
    value.

    Returns False if at least one value / row does not meet requirements. True otherwise.

    Since this process is memory intensive, it will return a cache value if this function has been called before.
    If calculation needs to be repeated use `force = True`

    :rtype bool
    """
    # TODO this could be split up into the visium portion and the general portion
    # NOTE: "if force" removed in pytest version.
    raw_layer_exists = None
    visium_and_is_single_true_matrix_size, visium_error_suffix = visium_and_is_single_true_matrix_size_with_error_suffix
    # Get potential raw_X
    if raw_x.dtype != np.float32:
        check.fail("Raw matrix values must have type numpy.float32.")
        return False

    matrix_format = get_matrix_format(raw_x)
    assert matrix_format != "unknown"
    raw_layer_exists = True
    is_sparse_matrix = matrix_format in SPARSE_MATRIX_TYPES

    if is_visium_and_is_single_true and raw_x.shape[0] != visium_and_is_single_true_matrix_size:
        raw_layer_exists = False
        check.fail(
            f"When {visium_error_suffix}, the raw matrix must be the "
            f"unfiltered feature-barcode matrix 'raw_feature_bc_matrix'. It must have exactly "
            f"{visium_and_is_single_true_matrix_size} rows. Raw matrix row count is "
            f"{raw_x.shape[0]}."
        )

    if is_visium_and_is_single_true and "in_tissue" in adata.obs and 0 in adata.obs["in_tissue"].values:
        raw_layer_exists = validate_raw_data_with_in_tissue_0(adata, raw_x, is_sparse_matrix, raw_layer_exists)
    else:
        raw_layer_exists = validate_raw_data(raw_x, is_sparse_matrix, raw_layer_exists)

    return raw_layer_exists


def validate_raw_data(x: DaskArray, is_sparse_matrix: bool, raw_layer_exists: bool) -> bool:
    """
    Validates the data values in the raw matrix. Matrix size is chunked for large matrices.

    :param x: raw matrix
    :param is_sparse_matrix: bool indicating if the matrix is sparse {csc, csr, coo}
    """

    def validate_chunk(matrix_chunk: Union[np.ndarray, sparse.spmatrix], is_sparse_matrix: bool) -> np.array:
        chunk_has_row_of_zeros = False
        chunk_has_invalid_nonzero_value = False
        if not chunk_has_row_of_zeros:
            if is_sparse_matrix:
                row_indices, _ = matrix_chunk.nonzero()
                if len(set(row_indices)) != matrix_chunk.shape[0]:
                    chunk_has_row_of_zeros = True
            # else, must be dense matrix, confirm that all rows have at least 1 nonzero value
            elif not all(np.apply_along_axis(np.any, axis=1, arr=matrix_chunk)):
                chunk_has_row_of_zeros = True

        if not chunk_has_invalid_nonzero_value and matrix_has_invalid_nonzero_values(matrix_chunk):
            chunk_has_invalid_nonzero_value = True

        return np.array([np.array([chunk_has_row_of_zeros, chunk_has_invalid_nonzero_value], dtype=object)])

    if len(x.chunks[0]) > 1:
        results = map_blocks(validate_chunk, x, is_sparse_matrix, dtype=object).compute()
        # Combine the results from all chunks
        has_row_of_zeros = any(chunk_result[0] for chunk_result in results)
        has_invalid_nonzero_value = any(chunk_result[1] for chunk_result in results)
    else:
        has_row_of_zeros, has_invalid_nonzero_value = validate_chunk(x.compute(), is_sparse_matrix)[0]

    if has_row_of_zeros:
        raw_layer_exists = False
        check.fail("Each cell must have at least one non-zero value in its row in the raw matrix.")
    if has_invalid_nonzero_value:
        raw_layer_exists = False
        check.fail("All non-zero values in raw matrix must be positive integers of type numpy.float32.")
    return raw_layer_exists


def validate_raw_data_with_in_tissue_0(
    adata: AnnData, x: DaskArray, is_sparse_matrix: bool, raw_layer_exists: bool
) -> bool:
    """
    Special case validation checks for Visium data with is_single = True and in_tissue column in obs where in_tissue
    has at least one value 0.

    :param x: raw matrix
    :param is_sparse_matrix: bool indicating if the matrix is sparse
    """

    def validate_chunk(
        matrix_chunk: Union[np.ndarray, sparse.spmatrix], is_sparse_matrix: bool, block_info: dict = None
    ) -> np.array:
        chunk_has_tissue_0_non_zero_row = False
        chunk_has_tissue_1_zero_row = False
        chunk_has_invalid_nonzero_values = False
        chunk_start_row = block_info[0]["array-location"][0][0] if (block_info and block_info.get(0)) else 0
        if matrix_has_invalid_nonzero_values(matrix_chunk):
            chunk_has_invalid_nonzero_values = True
        if is_sparse_matrix:
            nonzero_row_indices, _ = matrix_chunk.nonzero()
        else:  # must be dense matrix
            nonzero_row_indices = np.where(np.any(matrix_chunk != 0, axis=1))[0]
        for i in range(matrix_chunk.shape[0]):
            if chunk_has_tissue_0_non_zero_row and chunk_has_tissue_1_zero_row:
                # exit inner loop early
                break
            unchunked_i = i + chunk_start_row
            if (
                not chunk_has_tissue_0_non_zero_row
                and i in nonzero_row_indices
                and adata.obs["in_tissue"].iloc[unchunked_i] == 0
            ):
                chunk_has_tissue_0_non_zero_row = True
            elif (
                not chunk_has_tissue_1_zero_row
                and i not in nonzero_row_indices
                and adata.obs["in_tissue"].iloc[unchunked_i] == 1
            ):
                chunk_has_tissue_1_zero_row = True
        return np.array(
            [
                np.array(
                    [
                        chunk_has_tissue_0_non_zero_row,
                        chunk_has_tissue_1_zero_row,
                        chunk_has_invalid_nonzero_values,
                    ],
                    dtype=object,
                )
            ]
        )

    if len(x.chunks[0]) > 1:
        results = map_blocks(validate_chunk, x, is_sparse_matrix, dtype=object).compute()
        # Combine the results from all chunks
        has_tissue_0_non_zero_row = any(chunk_result[0] for chunk_result in results)
        has_tissue_1_zero_row = any(chunk_result[1] for chunk_result in results)
        has_invalid_nonzero_values = any(chunk_result[2] for chunk_result in results)
    else:
        has_tissue_0_non_zero_row, has_tissue_1_zero_row, has_invalid_nonzero_values = validate_chunk(
            x.compute(), is_sparse_matrix
        )[0]

    if not has_tissue_0_non_zero_row:
        raw_layer_exists = False
        check.fail(
            "If obs['in_tissue'] contains at least one value 0, then there must be at least "
            "one row with obs['in_tissue'] == 0 that has a non-zero value in the raw matrix."
        )
    if has_tissue_1_zero_row:
        raw_layer_exists = False
        check.fail(
            "Each observation with obs['in_tissue'] == 1 must have at least one "
            "non-zero value in its row in the raw matrix."
        )
    if has_invalid_nonzero_values:
        raw_layer_exists = False
        check.fail("All non-zero values in raw matrix must be positive integers of type numpy.float32.")
    return raw_layer_exists


def matrix_has_invalid_nonzero_values(x: Union[np.ndarray]) -> bool:
    """
    Checks whether the matrix has invalid non-zero values. The matrix must have all non-zero values as positive
    integers (type is numpy.float32).

    :param x: The matrix to validate

    :rtype bool
    :return True if any non-zero values are invalid given validation rules, False otherwise
    """
    data = x if isinstance(x, np.ndarray) else x.data
    return np.any((data % 1 > 0) | (data < 0))


@pytest.fixture(scope="session")
def is_visium(
    adata,
) -> bool:
    """
    Determine if the assay_ontology_term_id is Visium (EFO:0010961).

    :return True if assay_ontology_term_id is Visium, False otherwise.
    :rtype bool
    """
    assay_ontology_term_id = adata.obs.get("assay_ontology_term_id")
    is_visium = assay_ontology_term_id is not None and (assay_ontology_term_id == ASSAY_VISIUM).any()
    return is_visium


@pytest.fixture(scope="session")
def is_visium_including_descendants(adata, is_visium) -> bool:
    """
    Determine if the assay_ontology_term_id is Visium (inclusive descendant of EFO:0010961).
    Returns True if ANY assay_ontology_term_id is a Visium descendant

    :return True if assay_ontology_term_id is Visium, False otherwise.
    :rtype bool
    """
    _assay_key = "assay_ontology_term_id"

    # only compute if not already stored
    if is_visium is None and _assay_key in adata.obs.columns:
        # check if any assay_ontology_term_ids are descendants of VISIUM
        is_visium = bool(
            adata.obs[_assay_key]
            .astype("string")
            .apply(lambda assay: is_ontological_descendant_of(ONTOLOGY_PARSER, assay, ASSAY_VISIUM, True))
            .astype(bool)
            .any()
        )

    return is_visium
