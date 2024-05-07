import collections
import json
import re
from copy import deepcopy
from typing import Mapping, Tuple

import anndata as ad
import cerberus
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from cellxgene_schema import gencode
from tests.fixtures.examples_validate import adata, h5ad_valid
from utils import SPARSE_MATRIX_TYPES, get_matrix_format


class CustomErrorHandler(cerberus.errors.BasicErrorHandler):
    def add(self, error):
        """
        overriding because the deepcopy need to be skipped to avoid h5 Errors
        """
        self._rewrite_error_path(error)

        if error.is_logic_error:
            self._insert_logic_error(error)
        elif error.is_group_error:
            self._insert_group_error(error)
        elif error.code in self.messages:
            self._insert_error(error.document_path, self._format_message(error.field, error))

    def _format_message(self, field, error):
        prefix = "Warning" if "warn" in error.schema_path else "Error"
        return f"{prefix}: " + self.messages[error.code].format(
            *error.info, constraint=error.constraint, field=field, value=error.value
        )


class MyValidator(cerberus.Validator):
    types_mapping = cerberus.Validator.types_mapping.copy()
    types_mapping["anndata"] = cerberus.TypeDefinition("anndata", (ad.AnnData,), ())
    types_mapping["dataframe"] = cerberus.TypeDefinition("dataframe", (pd.DataFrame,), ())
    types_mapping["ndarray"] = cerberus.TypeDefinition("ndarray", (np.ndarray,), ())

    def __init__(self, *args, **kwargs):
        self.gene_checkers = {}
        super(MyValidator, self).__init__(*args, **kwargs)

    def _validate_forbidden(self, forbidden_values, field, value):
        """{'type': ['list', 'dict']}"""
        if isinstance(forbidden_values, dict) and isinstance(value, str):
            case_sensative = forbidden_values.get("case_sensative", False)
            suffix = forbidden_values.get("suffix", "")
            prefix = forbidden_values.get("prefix", "")
            regex = forbidden_values.get("regex", "")
            values = forbidden_values.get("values", [])
            regex_flags = 0
            if case_sensative:
                suffix = suffix.lower()
                prefix = prefix.lower()
                value = value.lower()
                regex_flags = re.IGNORECASE
                values = [v.lower() for v in values]
            if suffix and value.endswith(suffix):
                self._error(field, f"Value '{value}' must not end with '{suffix}'.")
            if prefix and value.startswith(prefix):
                self._error(field, f"Value '{value}' must not start with '{prefix}'.")
            if regex and re.match(regex, value, flags=regex_flags):
                self._error(field, f"Value '{value}' must not match the regex '{regex}'.")
            if values and value in values:
                self._error(field, f"Value '{value}' must not be one of '{values}'.")
        super(MyValidator, self)._validate_forbidden(forbidden_values, field, value)

    def _validate_allowed(self, allowed_values, field, value):
        """{'type': ['list', 'dict']}"""
        if isinstance(allowed_values, dict) and isinstance(value, str):
            case_sensative = allowed_values.get("case_sensative", False)
            suffix = allowed_values.get("suffix", "")
            prefix = allowed_values.get("prefix", "")
            regex = allowed_values.get("regex", "")
            values = allowed_values.get("values", [])
            regex_flags = 0
            if case_sensative:
                suffix = suffix.lower()
                prefix = prefix.lower()
                value = value.lower()
                regex_flags = re.IGNORECASE
                values = [v.lower() for v in values]
            if suffix and not value.endswith(suffix):
                self._error(field, f"Value '{value}' must end with '{suffix}'.")
            if prefix and not value.startswith(prefix):
                self._error(field, f"Value '{value}' must start with '{prefix}'.")
            if regex and not re.match(regex, value, flags=regex_flags):
                self._error(field, f"Value '{value}' must match the regex '{regex}'.")
            if values and value not in values:
                self._error(field, f"Value '{value}' must be one of '{values}'.")
        super(MyValidator, self)._validate_allowed(allowed_values, field, value)

    def _validate_attributes_schema(self, schemas: dict, field: str, _object: object) -> None:
        """
        The rule's arguments are validated against this schema:
        {'type': 'dict'}
        """
        if isinstance(_object, object):
            validator = self._get_child_validator(
                document_crumb=field,
                schema_crumb=(field, "attributes_schema"),
                schema=schemas,
            )
            document = dict()
            for k in schemas:
                try:
                    value = getattr(_object, k)
                except AttributeError:
                    if schemas[k].get("required", False):
                        self._error(field, f"Attribute '{k}' not found in {field}.")
                else:
                    document[k] = value
            if not validator(document, normalize=False):
                self._error(validator._errors)

    def _validate_index_schemas(self, schemas, field, _list):
        """
        The rule's arguments are validated against this schema:
        {'type': 'dict'}
        """
        if hasattr(_list, "__getitem__"):
            validator = self._get_child_validator(
                document_crumb=field,
                schema_crumb=(field, "index_schemas"),
                schema=schemas,
            )
            document = dict()
            for k in schemas:
                try:
                    value = _list[k]
                except IndexError as ex:
                    self._error(field, ex.args[0])
                    continue
                else:
                    document[k] = value
            if not validator(document, normalize=False):
                self._error(validator._errors)

    def _validate_warn(self, schema, field, value):
        """
        The rule's arguments are validated against this schema:
        {'type': 'dict'}
        """
        # all errrors will be marked as warnings.
        validator = self._get_child_validator(
            schema_crumb=("warn"),
            schema={field: schema},
        )
        if not validator({field: value}, normalize=False):
            self._error(validator._errors)

    def _validate_columns(self, column_schemas, field: str, df: pd.DataFrame):
        """
        The rule's arguments are validated against this schema:
        {'type': 'list'}
        """

        validator = self._get_child_validator(
            document_crumb=field,
            schema_crumb=(field, "columns"),
            schema=column_schemas,
        )
        document = dict()
        for k in column_schemas:
            try:
                value = df[k]
            except KeyError:
                if column_schemas[k].get("required", False):
                    self._error(field, f"Required column '{k}' not found in {field}.")
            else:
                document[k] = value
        if not validator(document, normalize=False):
            self._error(validator._errors)

    def _validate_index(self, index_schema: dict, field: str, df: pd.DataFrame):
        """
        The rule's arguments are validated against this schema:
        {'type': 'dict'}
        """
        validator = self._get_child_validator(
            document_crumb=field,
            schema_crumb=(field, "index"),
            schema={"index": index_schema},
        )
        document = dict(index=pd.Series(df.index))
        if not validator(document, normalize=False):
            self._error(validator._errors)

    def _validate_encoding_version(self, contraint: str, field: str, value: ad.AnnData) -> None:
        """
        The rule's arguments are validated against this schema:
        {'type': 'string'}
        """
        import h5py

        try:
            with h5py.File(value.filename, "r") as f:
                encoding_dict = dict(f.attrs)
                encoding_version = encoding_dict.get("encoding-version")
                if encoding_version != contraint:
                    self._error(field, "The h5ad artifact was generated with an AnnData version different from 0.8.0.")
        except Exception as e:
            self._error(field, f"Error validating encoding of h5ad file: {e}")

    def _validate_forbidden_attributes(self, forbidden, field, value):
        """
        The rule's arguments are validated against this schema:
        {'type': 'list'}
        """
        for key in forbidden:
            try:
                getattr(value, key)
                self._error(field, f"Attribute '{key}' is forbidden.")
            except AttributeError:
                pass

    def _validate_curie(self, constraint: dict, field: str, value) -> None:
        """
        The rule's arguments are validated against this schema:
        {'type': 'dict'}
        """
        # Check for NaN values
        if value.isnull().any():
            self._errors(field, "must not contain NaN values.")
            return

        if constraint:
            for term_str in value.drop_duplicates():
                self._validate_curie_str(term_str, field, constraint)

    def _validate_dtype(self, criteria: dict, column_name: str, column: pd.Series):
        """
        The rule's arguments are validated against this schema:
        {'type': 'dict'}
        """
        dtype = criteria.get("type")
        if dtype == "boolean" and column.dtype != bool:
            self.errors.append(f"Column '{column_name}' must be boolean, not '{column.dtype.name}'.")
        elif dtype == "categorical":
            if column.dtype.name != "category":
                self.errors.append(f"Column '{column_name}' must be categorical, not {column.dtype.name}.")
            else:
                if criteria.get("subtype") == "string":
                    if column.dtype.categories.dtype not in ["object", "string"]:
                        self.errors.append(
                            f"Column '{column_name}' must be object or string, not" f" {column.dtype.categories.dtype}."
                        )
                    else:
                        if any(len(cat.strip()) == 0 for cat in column.dtype.categories):
                            self.errors.append(f"Column '{column_name}' must not contain empty values.")

                # check for null values--skip on column defs with enums, since it will already be part of that check
                if not criteria.get("enum") and column.isnull().any():
                    self.errors.append(f"Column '{column_name}' must not contain NaN values.")
        elif criteria.get("kind") and column.dtype.kind != criteria["kind"]:
            self.errors.append(
                f"Column '{column_name}' must be of kind '{criteria['kind']}', not '{column.dtype.kind}'."
            )

    def _check_with_match_obs_columns(self, field, value):
        for i in value:
            if i not in self.root_document["adata"].obs.columns:
                self._error(value, f"Value '{i}' of list '{field}' is not a column in 'adata.obs'.")

    def _check_with_match_obsm_key(self, field, value) -> None:
        if value not in self.root_document["adata"].obsm:
            self.errors.append(f"'{value}' in '{field}' is not valid, it must be a key of 'adata.obsm'.")

    @property
    def obs_category_mapping(self):
        if not hasattr(self, "_obs_category_mapping"):
            # Check for categorical dtypes in the dataframe directly
            category_mapping = {}
            df = self.root_document["adata"].obs
            for column_name in df.columns:
                column = df[column_name]
                if column.dtype.name == "category":
                    category_mapping[column_name] = column.nunique()
            self._obs_category_mapping = category_mapping
        return self._obs_category_mapping

    def _check_with_obs_category_mapping(self, field: str, value: np.ndarray) -> None:
        category_mapping = self.obs_category_mapping
        if not isinstance(field, str):
            self._error(field, "Field must be a string.")
            return
        column_name = field.replace("_colors", "").replace("_ontology_term_id_colors", "")
        obs_unique_values = category_mapping.get(column_name)
        if not obs_unique_values:
            self._error(field, f"Colors field uns['{field}'] does not have a corresponding categorical field in obs.")
            return
        if len(value) < obs_unique_values:
            self._error(
                field,
                f"Annotated categorical field {field.replace('_colors', '')} must have at least "
                f"{obs_unique_values} color options "
                f"in uns[{field}]. Found: {value}",
            )

    def _check_with_annotation_mapping(self, field: str, value: Mapping) -> None:
        """The size of the ndarray stored for a key MUST NOT be zero"""
        for key, _value in value.items():
            # Check for empty ndarrays
            if isinstance(_value, np.ndarray) and not _value.size:
                self._errors(field, f"The size of the ndarray stored for a 'adata.{field}['{key}']' MUST NOT be zero.")

    def _check_with_equal_to_X_rows(self, field, value):
        if value != self.root_document["adata"].X.shape[0]:
            self._error(field, "must be equal to the number of rows in 'adata.X'.")

    def _check_with_ndarray_not_any_ninf(self, field, value):
        if isinstance(value, np.ndarray) and np.any(np.isneginf(value)):
            self._error(field, "Array contains negative infinity values.")

    def _check_with_ndarray_not_any_inf(self, field, value):
        if isinstance(value, np.ndarray) and np.any(np.isinf(value)):
            self._error(field, "Array contains infinity values.")

    def _check_with_ndarray_not_all_nan(self, field, value):
        if isinstance(value, np.ndarray) and np.any(np.isnan(value)):
            self._error(field, "Array contains NaN values.")

    def _check_with_unique(self, field, value):
        if isinstance(value, (pd.Index, pd.Series)) and value.nunique() != len(value):
            duplicates = [item for item, count in collections.Counter(value).items() if count > 1]
            self._error(field, f"Values must be unique. Found duplicates: {duplicates}.")

    def _check_with_feature_id(self, field: str, feature_ids):
        """
        Validates feature ids, i.e. checks that they are present in the reference
        If there are any errors, it adds them to self._error and adds it to the list of invalid features

        :param str feature_ids: the feature ids to be validated
        :param str field: the name of the field or attribute

        """
        for feature_id in feature_ids:
            organism = gencode.get_organism_from_feature_id(feature_id)

            if not organism:
                self._error(
                    f"Could not infer organism from feature ID '{feature_id}' in '{field}', "
                    f"make sure it is a valid ID."
                )
                continue

            if organism not in self.gene_checkers:
                self.gene_checkers[organism] = gencode.GeneChecker(organism)

            if not self.gene_checkers[organism].is_valid_id(feature_id):
                self._error(f"'{feature_id}' is not a valid feature ID in '{field}'.")

    def _check_with_feature_is_filtered(self, field: str, column: pd.Series):
        """
        Validates the "is_feature_filtered" in adata.var. This column must be bool, and for genes that are set to
        True, their expression values in X must be 0.
        If there are any errors, it adds them to self._error.

        :rtype none
        """
        adata = self.root_document["adata"]
        if sum(column) > 0:
            n_nonzero = 0

            x_format = get_matrix_format(adata, adata.X)
            if x_format in SPARSE_MATRIX_TYPES:
                n_nonzero = adata.X[:, column].count_nonzero()
            elif x_format == "dense":
                n_nonzero = np.count_nonzero(adata.X[:, column])
            else:
                self._error(
                    field,
                    f"X matrix is of type {type(adata.X)}, validation of 'feature_is_filtered' "
                    f"cannot be completed.",
                )

            if n_nonzero > 0:
                self._error(
                    field,
                    f"Some features are 'True' in '{field}', but there are "
                    f"{n_nonzero} non-zero values in the corresponding columns of the matrix 'X'. All values for "
                    f"these features must be 0.",
                )

    def _check_with_is_supported_spatial_assay(self) -> bool:
        """
        Determine if the assay_ontology_term_id is either Visium (EFO:0010961) or Slide-seqV2 (EFO:0030062).

        :return True if assay_ontology_term_id is Visium or Slide-seqV2, False otherwise.
        :rtype bool
        """
        ASSAY_VISIUM = "EFO:0010961"
        ASSAY_SLIDE_SEQV2 = "EFO:0030062"

        if self.is_spatial is None:
            try:
                self.is_spatial = False
                if self.adata.obs.assay_ontology_term_id.isin([ASSAY_VISIUM, ASSAY_SLIDE_SEQV2]).any():
                    self.is_spatial = True
            except AttributeError:
                # specific error reporting will occur downstream in the validation
                self.is_spatial = False
        return self.is_spatial


def get_validator():
    obsm_schema = {
        "required": True,
        "anyof": [
            {
                "keysrules": {
                    "type": "string",
                    "regex": "^X_[a-zA-Z][a-zA-Z0-9_.-]*$",
                    "forbidden": {"values": ["x_spatial"], "case_sensative": True},
                },
                "valuesrules": {
                    "type": "ndarray",
                    "attributes_schema": {
                        "shape": {
                            "type": "list",
                            "minlength": 2,
                            "index_schemas": {
                                0: {"type": "integer", "check_with": "equal_to_X_rows"},
                                1: {"type": "integer", "min": 2},
                            },
                            "items": [
                                {"type": "integer", "check_with": "equal_to_X_rows"},
                                {"type": "integer", "min": 2},
                            ],
                        },
                        "dtype": {
                            "attributes_schema": {
                                "kind": {"type": "string", "allowed": ["f", "i", "u"]},
                            }
                        },
                        "size": {"type": "integer", "min": 1},
                    },
                    "check_with": ["ndarray_not_any_ninf", "ndarray_not_any_inf", "ndarray_not_all_nan"],
                },
            },
            {
                "keysrules": {
                    "type": "string",
                    "regex": "^(?!X_)[a-zA-Z][a-zA-Z0-9_.-]*$",
                },
                "valuesrules": {"type": "ndarray", "attributes_schema": {"size": {"type": "integer", "min": 1}}},
            },
        ],
    }
    var_schema_common = {
        "type": "dataframe",
        "index": {
            "check_with": ["unique", "feature_id"],
            "dtype": {"type": "string"},
        },
        "attributes_schema": {
            "shape": {"index_schemas": {0: {"warn": {"min": 2000}}}},  # row min length
        },
    }
    var_schema = deepcopy(var_schema_common)
    var_schema["attributes_schema"]["feature_is_filtered"] = {
        "required": True,
        "dtype": {"type": "boolean"},
        "check_with": "feature_is_filtered",
    }
    var_schema["required"] = True
    var_raw_schema = deepcopy(var_schema_common)
    # var_raw_schema["forbidden"] = ["feature_is_filtered"]

    uns_schema = {
        "type": "dict",
        "allow_unknown": True,
        "anyof": [
            {
                "keysrules": {
                    "type": "string",
                    "forbidden": [
                        "schema_version",
                        "citation",
                        "schema_reference",
                        "X_normalization",
                        "default_field",
                        "layer_descriptions",
                        "tags",
                        "version",
                        "contributors",
                        "preprint_doi",
                        "project_description",
                        "project_links",
                        "project_name",
                        "publication_doi",
                    ],
                },
                "valuesrules": {"anyof": [{"type": "boolean"}, {"empty": False}, {"nullable": True}]},
                "schema": {
                    "title": {"type": "string", "required": True},
                    "batch_condition": {
                        "type": ["list", "ndarray"],
                        "schema": {"type": "string"},
                        "check_with": "match_obs_columns",
                    },
                    "default_embedding": {"type": "string", "check_with": "match_obsm_key"},
                    "X_approximate_distribution": {"type": "string", "allowed": ["count", "normal"]},
                },
            },
            {
                "keysrules": {
                    "type": "string",
                    "regex": "^.*(_colors)|(_ontology_term_id_colors)$",
                    "forbidden": [
                        "assay_colors",
                        "cell_type_colors",
                        "development_stage_colors",
                        "disease_colors",
                        "organism_colors",
                        "self_reported_ethnicity_colors",
                        "sex_colors",
                        "tissue_colors",
                    ],
                },
                "valuesrules": {
                    "type": "ndarray",
                    "anyof": [
                        # Verify that either all colors are hex OR all colors are CSS4 named colors strings
                        {"schema": {"type": "string", "regex": "^#([0-9a-fA-F]{6})$"}},
                        {"schema": {"type": "string", "allowed": list(mcolors.CSS4_COLORS)}},
                    ],
                    "check_with": "obs_category_mapping",
                },
            },
        ],
    }
    schema = {
        "adata": {"type": "anndata", "required": True, "encoding_version": "0.1.0"},
        "obs": {
            "type": "dataframe",
            "required": True,
        },
        "var": var_schema,
        "obsm": obsm_schema,
        "obsp": {"check_with": "annotation_mapping"},
        "varm": {"check_with": "annotation_mapping"},
        "varp": {"check_with": "annotation_mapping"},
        "uns": uns_schema,
        "raw": {
            "attributes_schema": {
                # "X": {"type": "ndarray"},
                "var": var_raw_schema
            },
        },
    }
    return MyValidator(schema, error_handler=CustomErrorHandler())


def validate_anndata(adata: ad.AnnData, validator):
    document = dict(
        adata=adata,
        obs=adata.obs,
        var=adata.var,
        obsm=adata.obsm,
        obsp=adata.obsp,
        varm=adata.varm,
        varp=adata.varp,
        uns=adata.uns,
        raw=adata.raw,
    )
    if not validator.validate(document, normalize=False):
        errors = validator.errors
        print(json.dumps(separate_messages(errors), indent=4))


def separate_messages(data: dict) -> Tuple[dict, dict]:
    errors = collections.defaultdict(list)
    warnings = collections.defaultdict(list)

    def traverse_dict(d: dict, path: list):
        for key, value in d.items():
            path.append(key if isinstance(key, str) else str(key))
            if isinstance(value, list):
                for item in value:
                    if isinstance(item, str):
                        message_path = ".".join(path)
                        if item.startswith("Error:"):
                            errors[message_path].append(item)
                        elif item.startswith("Warning:"):
                            warnings[message_path].append(item)
                    elif isinstance(item, dict):
                        traverse_dict(item, path)
            elif isinstance(value, dict):
                traverse_dict(value, path)
            path.pop()

    traverse_dict(data, [])
    return errors, warnings


if __name__ == "__main__":
    VALIDATOR = get_validator()
    validate_anndata(adata, VALIDATOR)

    _adata = ad.read_h5ad(h5ad_valid, backed="r")
    validate_anndata(_adata, VALIDATOR)

    _adata.uns["title"] = [1, 2, 3]
    _adata.uns[1] = None
    _adata.uns["asdfads"] = "asfasdf"
    _adata.uns["project_name"] = "project_name"
    _adata.uns["X_approximate_distribution"] = "asdf"
    _adata.obsm["X_spatial"] = _adata.obsm["X_pca"]
    _adata.obsm["abcd"] = _adata.obsm["X_pca"]
    validate_anndata(_adata, VALIDATOR)
