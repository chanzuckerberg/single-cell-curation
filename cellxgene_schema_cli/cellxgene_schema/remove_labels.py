import logging

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
        From a valid (per cellxgene's schema) adata, this function parses the associated schema definition for
        'add_labels.to_column' annotations, and removes specified columns from self.adata
        """
        for component_name in ["obs", "var", "raw.var", "uns"]:
            # If the component does not exist, skip (this is for raw.var)
            component = getattr_anndata(self.adata, component_name)
            if component is None:
                continue

            # Doing it for columns
            if "columns" in self.schema_def["components"][component_name]:
                for column_def in self.schema_def["components"][component_name]["columns"].values():
                    if "add_labels" in column_def:
                        self._remove_columns(component, column_def)

            # Remove automatically annotated columns
            if "reserved_columns" in self.schema_def["components"][component_name]:
                for field in self.schema_def["components"][component_name]["reserved_columns"]:
                    del component[field]

            # Doing it for index
            if "index" in self.schema_def["components"][component_name]:
                index_def = self.schema_def["components"][component_name]["index"]
                if "add_labels" in index_def:
                    self._remove_columns(component, index_def)

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
