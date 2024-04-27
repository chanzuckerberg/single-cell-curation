import json
from typing import Mapping

import anndata as ad
import cerberus
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from tests.fixtures.examples_validate import h5ad_valid


class CustomerErrorHandler(cerberus.errors.BasicErrorHandler):
    def add(self, error):
        self._rewrite_error_path(error)

        if error.is_logic_error:
            self._insert_logic_error(error)
        elif error.is_group_error:
            self._insert_group_error(error)
        elif error.code in self.messages:
            self._insert_error(error.document_path, self._format_message(error.field, error))


anndata_type = cerberus.TypeDefinition("anndata", (ad.AnnData,), ())
dataframe_type = cerberus.TypeDefinition("dataframe", (pd.DataFrame,), ())
ndarray_type = cerberus.TypeDefinition("ndarray", (np.ndarray,), ())


class MyValidator(cerberus.Validator):
    types_mapping = cerberus.Validator.types_mapping.copy()
    types_mapping["anndata"] = anndata_type
    types_mapping["dataframe"] = dataframe_type
    types_mapping["ndarray"] = ndarray_type

    def _validate_attributes(self, schemas: dict, field: str, _object: object) -> None:
        """
        The rule's arguments are validated against this schema:
        {'type': 'dict'}
        """
        if isinstance(_object, object):
            validator = self._get_child_validator(
                document_crumb=field,
                schema_crumb=(field, "attributes"),
                schema=schemas,
            )
            document = dict()
            for k in schemas:
                try:
                    value = getattr(_object, k)
                except AttributeError:
                    self._error(field, f"Attribute '{k}' not found in {field}.")
                    continue
                else:
                    document[k] = value
            if not validator(document, normalize=False):
                self._error(validator._errors)

    def _validate_match_obs_columns(self, match_obs_columns, field, value):
        """
        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if match_obs_columns:
            for i in value:
                if i not in self.root_document["adata"].obs.columns:
                    self._error(value, f"Value '{i}' of list '{field}' is not a column in 'adata.obs'.")

    def _validate_match_obsm_key(self, match_obsm_key, field, value) -> None:
        """
        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if match_obsm_key and value not in self.root_document["adata"].obsm:
            self.errors.append(f"'{value}' in '{field}' is not valid, it must be a key of 'adata.obsm'.")

    def _validate_encoding_version(self, contraint: str, field: str, value: ad.AnnData) -> None:
        """
        The rule's arguments are validated against this schema:
        {'type': 'string'}
        """
        import h5py

        with h5py.File(value.filename, "r") as f:
            encoding_dict = dict(f.attrs)
            encoding_version = encoding_dict.get("encoding-version")
            if encoding_version != contraint:
                self._error(field, "The h5ad artifact was generated with an AnnData version different from 0.8.0.")

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
        category_mapping = self.get_category_mapping()
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
        for key, _value in value.items():
            # Check for empty ndarrays
            if isinstance(_value, np.ndarray) and not _value.size:
                self._errors(field, f"The size of the ndarray stored for a 'adata.{field}['{key}']' MUST NOT be zero.")


def validate_anndata():
    schema = {
        "adata": {
            "type": "anndata",
            "required": True,
            "encoding_version": "0.1.0",
            "attributes": {
                "obs": {"type": "dataframe", "required": True},
                "var": {"type": "dataframe", "required": True},
                "obsm": {"check_with": "annotation_mapping"},
                "obsp": {"check_with": "annotation_mapping"},
                "varm": {"check_with": "annotation_mapping"},
                "varp": {"check_with": "annotation_mapping"},
                "uns": {"type": "dict", "required": True},
            },
        },
        "uns": {
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
                            "type": "list",
                            "schema": {"type": "string"},
                            "match_obs_columns": True,
                        },
                        "default_embedding": {"type": "string", "match_obsm_key": True},
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
        },
    }
    adata = ad.read_h5ad(h5ad_valid, backed="r")
    v = MyValidator(schema, error_handler=CustomerErrorHandler())
    if not v.validate(dict(adata=adata), normalize=False):
        print(json.dumps(v.errors, sort_keys=True, indent=4))

    adata.uns["title"] = [1, 2, 3]
    adata.uns[1] = None
    adata.uns["asdfads"] = "asfasdf"
    adata.uns["project_name"] = "project_name"
    adata.uns["X_approximate_distribution"] = "asdf"
    if not v.validate(dict(adata=adata, uns=adata.uns), normalize=False):
        print(json.dumps(v.errors, indent=4))


validate_anndata()
