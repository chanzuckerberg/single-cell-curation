import logging
import os
import re
import sys
from base64 import b85encode
from functools import lru_cache
from typing import Any, Dict, List, Optional, Tuple, Union

import anndata as ad
import h5py
import numpy as np
from anndata.compat import DaskArray
from anndata.experimental import read_dispatched, read_elem_as_dask
from cellxgene_ontology_guide.ontology_parser import OntologyParser
from cellxgene_schema.matrix_utils import calculate_matrix_nonzero, check_non_csr_matrixes, determine_matrix_format
from xxhash import xxh3_64_intdigest

logger = logging.getLogger(__name__)

KB = 1024
MB = 1024 * KB
GB = 1024 * MB

"""
Ideally, these methods should all only live within matrix_utils. However, we
currently import these into single-cell-data-portal, so we need to keep these
in here for backwards compatibility until we can refactor that.
"""
SPARSE_MATRIX_TYPES = {"csr", "csc", "coo"}
SUPPORTED_SPARSE_MATRIX_TYPES = {"csr"}


def get_matrix_format(matrix: DaskArray) -> str:
    return determine_matrix_format(matrix)


def count_matrix_nonzero(matrix: DaskArray) -> int:
    return calculate_matrix_nonzero(matrix)


def check_non_csr_matrices(adata: ad.AnnData):
    return check_non_csr_matrixes(adata)


def replace_ontology_term(dataframe, ontology_name, update_map):
    column_name = f"{ontology_name}_ontology_term_id"
    if dataframe[column_name].dtype != "category":
        dataframe[column_name] = dataframe[column_name].astype("category")
    for old_term, new_term in update_map.items():
        if old_term in dataframe[column_name].cat.categories:
            # add new one if not already in category, else continue
            if new_term not in dataframe[column_name].cat.categories:
                dataframe[column_name] = dataframe[column_name].cat.add_categories(new_term)
            # replace in dataset
            dataframe.loc[dataframe[column_name] == old_term, column_name] = new_term
            # remove deprecated_term from category
            dataframe[column_name] = dataframe[column_name].cat.remove_categories(old_term)


def replace_ontology_term_uns(adata: ad.AnnData, ontology_name, update_map) -> ad.AnnData:
    key_name = f"{ontology_name}_ontology_term_id"

    for old_term, new_term in update_map.items():
        if adata.uns[key_name] == old_term:
            adata.uns[key_name] = new_term
    return adata


def replace_delimiter(dataframe, old_delimiter: str, new_delimiter: str, column_name: str):
    """
    Replace the delimiter in a column of a dataframe.

    :param dataframe: The dataframe containing the column to modify.
    :param old_delimiter: The delimiter to replace.
    :param new_delimiter: The new delimiter to use.
    :param column_name: The name of the column to modify.
    """
    if column_name not in dataframe.columns:
        raise KeyError(f"Column '{column_name}' not found")

    dataframe[column_name] = dataframe[column_name].str.replace(old_delimiter, new_delimiter, regex=False)


def move_column_from_obs_to_uns(adata: ad.AnnData, column_name: str) -> ad.AnnData:
    if column_name not in adata.obs:
        logger.warning(f"Column '{column_name}' not found in adata.obs, cannot migrate from obs to uns")
        return adata

    values = adata.obs[column_name].unique()

    if len(values) != 1:
        raise ValueError(f"Cannot migrate from obs to uns because '{column_name}' has multiple values: {values}")

    adata.uns[column_name] = values[0]
    del adata.obs[column_name]

    return adata


def map_ontology_term(dataframe, ontology_name, map_from_column, update_map):
    column_name = f"{ontology_name}_ontology_term_id"
    if dataframe[column_name].dtype != "category":
        dataframe[column_name] = dataframe[column_name].astype("category")
    for map_value, new_term in update_map.items():
        if new_term not in dataframe[column_name].cat.categories:
            dataframe[column_name] = dataframe[column_name].cat.add_categories(new_term)
        dataframe.loc[dataframe[map_from_column] == map_value, column_name] = new_term
    dataframe[column_name] = dataframe[column_name].cat.remove_unused_categories()


def remove_deprecated_features(*, adata: ad.AnnData, deprecated: List[str]) -> ad.AnnData:
    # Filter out genes that don't appear in the approved annotation
    var_to_keep = adata.var.index[~adata.var.index.isin(deprecated)].tolist()
    adata = adata[:, var_to_keep].copy()

    # Repeat much of the same steps for the raw.var, if it exists
    if adata.raw:
        raw_adata = ad.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs)
        var_to_keep = raw_adata.var.index[~raw_adata.var.index.isin(deprecated)].tolist()
        raw_adata = raw_adata[:, var_to_keep].copy()
        adata.raw = raw_adata
    return adata


def remap_deprecated_features(*, adata: ad.AnnData, remapped_features: Dict[str, str]) -> ad.AnnData:
    # Use remapped_terms to map to the new term ids
    adata.var.index = [remapped_features.get(val, val) for val in adata.var.index]

    # Repeat much of the same steps for the raw.var, if it exists
    if adata.raw:
        raw_adata = ad.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs)
        raw_adata.var.index = [remapped_features.get(val, val) for val in raw_adata.var.index]
        adata.raw = raw_adata
    return adata


def update_cell_line(dataframe, cellosaurus_term):
    cell_line_remapping = {
        "development_stage_ontology_term_id": "na",
        "sex_ontology_term_id": "na",
        "self_reported_ethnicity_ontology_term_id": "na",
        "donor_id": "na",
        "tissue_ontology_term_id": cellosaurus_term,
    }

    # first update all "cell culture" values to "cell line"
    if "cell line" not in dataframe["tissue_type"].cat.categories:
        dataframe["tissue_type"] = dataframe["tissue_type"].cat.add_categories("cell line")

    # replace in dataset
    dataframe.loc[dataframe["tissue_type"] == "cell culture", "tissue_type"] = "cell line"

    # remove deprecated_term from category
    dataframe["tissue_type"] = dataframe["tissue_type"].cat.remove_unused_categories()

    # per columns that now require na or cell_line info, set values appropriately
    for column, new_value in cell_line_remapping.items():
        if dataframe[column].dtype != "category":
            dataframe[column] = dataframe[column].astype("category")

        if new_value not in dataframe[column].cat.categories:
            dataframe[column] = dataframe[column].cat.add_categories(new_value)

        dataframe.loc[dataframe["tissue_type"] == "cell line", column] = new_value
        dataframe[column] = dataframe[column].cat.remove_unused_categories()


# Compiled regex pattern for parsing list indexing (e.g., "key[0]")
_LIST_INDEX_PATTERN = re.compile(r"^(.+?)\[(\d+)\]$")


def getattr_anndata(adata: ad.AnnData, attr: Optional[str] = None) -> Optional[Any]:
    """
    Extended getattr for AnnData objects that supports nested access using dot notation and list indexing.

    Supports:
    - Simple attributes: "obs", "var", "uns", "X", etc.
    - Nested dictionary access: "uns.key_name" -> adata.uns['key_name']
    - Column access: "obs.column_name" -> adata.obs['column_name']
    - List indexing: "uns.some_list[0]" -> adata.uns['some_list'][0]
    - Raw access: "raw.var", "raw.obs.column_name", etc.
    - Combined: "uns.some_list[0].nested_key" -> adata.uns['some_list'][0]['nested_key']

    Examples:
        getattr_anndata(adata, "obs") -> adata.obs
        getattr_anndata(adata, "obs.cell_type") -> adata.obs['cell_type']
        getattr_anndata(adata, "uns.species_ontology_term_id") -> adata.uns['species_ontology_term_id']
        getattr_anndata(adata, "uns.some_list[0]") -> adata.uns['some_list'][0]
        getattr_anndata(adata, "uns.some_list[0].nested_key") -> adata.uns['some_list'][0]['nested_key']
        getattr_anndata(adata, "raw.var") -> adata.raw.var
        getattr_anndata(adata, "raw.obs.cell_type") -> adata.raw.obs['cell_type']

    :param ad.AnnData adata: the anndata.AnnData object from which to extract an attribute
    :param Optional[str] attr: name of the attribute to extract, supports dot notation and bracket notation for nested access
    :return: the attribute or None if it does not exist
    :rtype: Optional[Any]
    """
    if attr is None or not attr:
        return None

    def parse_part(part: str) -> Tuple[str, Optional[int]]:
        """Parse a part that might contain bracket notation for list indexing."""
        match = _LIST_INDEX_PATTERN.match(part)
        if match:
            return match.group(1), int(match.group(2))
        return part, None

    def access_value(obj: Any, key: str, index: Optional[int] = None) -> Optional[Any]:
        """Access a value from an object, handling dict, list, DataFrame, or attribute access."""
        # First, get the value by key/attribute
        if isinstance(obj, dict):
            if key not in obj:
                return None
            value = obj[key]
        elif hasattr(obj, "columns"):
            # DataFrame - access column
            if key not in obj.columns:
                return None
            value = obj[key]
        else:
            # Try attribute access
            if not hasattr(obj, key):
                return None
            value = getattr(obj, key, None)
            if value is None:
                return None

        # If we have an index, apply it to the value
        if index is not None:
            if isinstance(value, (list, tuple)) or hasattr(value, "__getitem__"):
                try:
                    return value[index]
                except (IndexError, KeyError, TypeError):
                    return None
            return None

        return value

    # Split by dots to handle nested access
    parts = attr.split(".")
    if not parts or not any(parts):  # Handle empty string or only dots
        return None

    # Start with the AnnData object
    # Note: current changes type as we traverse (AnnData -> dict/DataFrame/list/etc.)
    current: Any = adata

    # Traverse the path
    for part in parts:
        if not part:  # Empty parts indicate invalid path (e.g., from "a..b")
            raise ValueError("Invalid attribute path: empty part found (consecutive dots not allowed)")

        # Check if we're accessing raw
        if part == "raw":
            if not hasattr(current, "raw") or current.raw is None:
                return None
            current = current.raw
            continue

        # Parse the part for potential list indexing
        key, index = parse_part(part)

        # Access the value
        current = access_value(current, key, index)
        if current is None:
            return None

    return current


def read_backed(f: h5py.File, chunk_size: int) -> ad.AnnData:
    """
    Read an AnnData object from a h5py.File object, reading in matrices (dense or sparse) as dask arrays. Does not
    read full matrices into memory.

    :param f: h5py.File object
    :param chunk_size: size of chunks to read matrices in
    :return: ad.AnnData object
    """

    def callback(func, elem_name: str, elem, iospec):
        if "/layers" in elem_name or "/uns" in elem_name or elem_name == "/X" or elem_name == "/raw/X":
            if iospec.encoding_type == "csr_matrix":
                n_vars = elem.attrs.get("shape")[1]
                return read_elem_as_dask(elem, chunks=(chunk_size, n_vars))
            elif iospec.encoding_type == "csc_matrix":
                n_obs = elem.attrs.get("shape")[0]
                return read_elem_as_dask(elem, chunks=(n_obs, chunk_size))
            elif iospec.encoding_type == "array" and len(elem.shape) == 2:
                n_vars = elem.shape[1]
                return read_elem_as_dask(elem, chunks=(chunk_size, n_vars))
            else:
                return func(elem)
        else:
            return func(elem)

    adata = read_dispatched(f, callback=callback)

    return adata


def read_h5ad(h5ad_path: Union[str, bytes, os.PathLike], chunk_size: int = 5000) -> ad.AnnData:
    """
    Reads h5ad into adata
    :params Union[str, bytes, os.PathLike] h5ad_path: path to h5ad to read

    :rtype None
    """
    try:
        f = h5py.File(h5ad_path)
        adata = read_backed(f, chunk_size)

    except (OSError, TypeError):
        logger.info(f"Unable to open '{h5ad_path}' with AnnData")
        sys.exit(1)

    return adata


def get_hash_digest_column(dataframe):
    """
    Get column with hash digest for each row in dataframe.
    """
    df_index = dataframe.index.to_series()
    assert df_index.is_unique
    return (
        df_index.map(xxh3_64_intdigest)
        .astype(np.uint64)
        .apply(lambda v: b85encode(v.to_bytes(8, "big")).decode("ascii"))
    )


@lru_cache()
def is_ontological_descendant_of(onto: OntologyParser, term: str, target: str, include_self: bool = True) -> bool:
    """
    Determines if :term is an ontological descendant of :target and whether to include :term==:target.

    This function is cached and is safe to call many times.

    #TODO:[EM] needs testing
    """
    return term in set(onto.get_term_descendants(target, include_self))


@lru_cache()
def get_descendants(onto: OntologyParser, term: str, include_self: bool = True) -> List[str]:
    return onto.get_term_descendants(term, include_self=True)


def get_chunks(*, step_size: int, total_size: int) -> List[tuple]:
    """
    Get chunks of size step_size from 0 to total_size

    :param step_size: size of each chunk
    :param total_size: total size of the object we're reading
    :return: list of tuples with start and end indices for each chunk
    """
    return [(i, min(i + step_size, total_size)) for i in range(0, total_size, step_size)]
