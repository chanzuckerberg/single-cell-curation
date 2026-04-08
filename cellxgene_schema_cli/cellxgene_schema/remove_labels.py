import logging
import re

import anndata as ad
from pandas import DataFrame

from . import schema
from .utils import getattr_anndata

logger = logging.getLogger(__name__)


class AnnDataLabelRemover:
    """
    From valid h5ad, handles writing a new h5ad file with appended ontology/gene labels removed
    from adata.obs, adata.var, and adata.raw.var respectively as indicated in the schema definition
    """

    def __init__(self, adata: ad.AnnData = None):
        self.adata = adata
        self.schema_def = schema.get_schema_definition()

    def remove_labels(self):
        """
        Removes specified columns and keys from self.adata based on the schema definition,
        including handling genetic_perturbations, nested reserved_keys, pre-analysis datasets,
        perturbation datasets, and experimental condition fields.
        """
        for component_name in ["obs", "var", "raw.var", "uns", "genetic_perturbations"]:
            component = getattr_anndata(self.adata, component_name)
            if component is None:
                continue

            component_def = self.schema_def["components"].get(component_name, {})

            # Remove columns defined by add_labels
            if "columns" in component_def:
                for column_def in component_def["columns"].values():
                    if "add_labels" in column_def:
                        self._remove_columns(component, column_def)

            # Remove automatically annotated columns and keys
            for key in ("reserved_columns", "reserved_keys"):
                self._remove_reserved(component, component_def.get(key, []))

            # Remove columns defined by index
            if "index" in component_def:
                index_def = component_def["index"]
                if "add_labels" in index_def:
                    self._remove_columns(component, index_def)

            # Remove any labels added as dict keys
            if "keys" in component_def:
                for key, key_def in component_def["keys"].items():
                    if "add_labels" in key_def:
                        for label_def in key_def["add_labels"]:
                            key_to_remove = label_def["to_key"]
                            if key_to_remove in component:
                                del component[key_to_remove]
                    # Remove nested reserved_keys if present
                    if key in component and isinstance(component[key], dict):
                        # Handle nested keys structure (e.g., genetic_perturbations with key_pattern)
                        if "keys" in key_def:
                            # This key has nested keys (e.g., genetic_perturbation_ids with key_pattern)
                            nested_key_def = key_def["keys"]
                            for nested_key_def_item in nested_key_def.values():
                                if "reserved_keys" in nested_key_def_item:
                                    # Process each item in the nested dict (e.g., each perturbation ID)
                                    # The nested_key might have a key_pattern, so we process all keys in component[key]
                                    for item_key in component[key]:
                                        if isinstance(component[key][item_key], dict):
                                            self._remove_reserved(
                                                component[key][item_key], nested_key_def_item["reserved_keys"]
                                            )
                        # Also check for reserved_keys directly in key_def
                        if "reserved_keys" in key_def:
                            self._remove_reserved(component[key], key_def["reserved_keys"])

        # Remove intended_features from each genetic_perturbation entry if present.
        # This field is not part of the schema but may exist in legacy datasets.
        genetic_perturbations = self.adata.uns.get("genetic_perturbations")
        if genetic_perturbations:
            for pert_data in genetic_perturbations.values():
                if isinstance(pert_data, dict) and "intended_features" in pert_data:
                    del pert_data["intended_features"]

    def _remove_reserved(self, component, reserved):
        """
        Recursively removes reserved fields from a component, including key_pattern support.
        """
        if isinstance(reserved, dict):
            # Handle key_pattern
            if "key_pattern" in reserved and isinstance(component, dict):
                pattern = re.compile(reserved["key_pattern"])
                keys_to_remove = [k for k in component if pattern.match(k)]
                for k in keys_to_remove:
                    del component[k]
            # Recurse for other keys
            for k, v in reserved.items():
                if k != "key_pattern" and k in component:
                    self._remove_reserved(component[k], v)
        elif isinstance(reserved, list):
            for field in reserved:
                if field in component:
                    del component[field]

    def _remove_columns(self, component: DataFrame, subcomponent_definition: dict):
        """
        Given an adata component and subcomponent definition, this function deletes all existing columns in the
        self.adata component that are defined as added labels ('add_labels.to_column') in the subcomponent definition.

        :param pd.Dataframe component: dataframe within adata dataset (i.e. 'obs', 'var', 'raw.var')
        :param dict subcomponent_definition: yaml-defined schema for subcomponent of component (i.e. index, or a
                    particular column)

        :rtype None
        """
        for label_def in subcomponent_definition["add_labels"]:
            column_name = label_def["to_column"]
            if column_name in component:
                del component[column_name]
