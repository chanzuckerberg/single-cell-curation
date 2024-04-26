import anndata as ad
import cerberus
from tests.fixtures.examples_validate import h5ad_valid

anndata_type = cerberus.TypeDefinition("anndata", (ad.AnnData,), ())
adata = ad.read_h5ad(h5ad_valid, backed="r")


class MyValidator(cerberus.Validator):
    types_mapping = cerberus.Validator.types_mapping.copy()
    types_mapping["anndata"] = anndata_type

    def _validate_match_obs_columns(self, match_obs_columns, field, value):
        """
        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if match_obs_columns:
            for i in value:
                if i not in self.adata.obs.columns:
                    self._error(value, f"Value '{i}' of list '{field}' is not a column in 'adata.obs'.")

    def _validate_match_obs_keys(self, match_obs_keys, field, value):
        """
        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if match_obs_keys:
            for i in value:
                if i not in self.adata.obs.keys():
                    self._error(value, f"Value '{i}' of list '{field}' is not a key in 'adata.obs'.")

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


def validate_anndata():
    schema = {
        "adata": {"type": "anndata", "required": True, "encoding_version": "0.1.0"},
        "uns": {
            "type": "dict",
            "keysrules": {"type": "string"},
            "valuesrules": {"anyof": [{"type": "boolean"}, {"empty": False}, {"nullable": True}]},
            "contains": {
                "title": {"type": "string", "required": True},
                "batch_condition": {"type": "list", "match_obs_columns": True},
                "default_embedding": {"match_obs_keys": True},
                "X_approximate_distribution": {"type": "string", "allowed": ["count", "normal"]},
            },
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
    }
    v = MyValidator(schema)
    if not v.validate(
        dict(
            adata=adata,
            uns=adata.uns,
        )
    ):
        print(v.errors)
    return v


v = validate_anndata()
