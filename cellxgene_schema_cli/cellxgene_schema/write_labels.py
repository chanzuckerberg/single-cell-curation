import logging
import traceback
from typing import Dict, List, Optional

import anndata
import pandas as pd

from . import gencode, schema
from .annotate_guides import annotate_perturbations_in_h5ad
from .env import SCHEMA_REFERENCE_BASE_URL, SCHEMA_REFERENCE_FILE_NAME
from .gencode import get_gene_checker
from .ontology_parser import ONTOLOGY_PARSER
from .utils import get_hash_digest_column, getattr_anndata

logger = logging.getLogger(__name__)


class AnnDataLabelAppender:
    """
    From valid h5ad, handles writing a new h5ad file with ontology/gene labels added
    to adata.obs and adata.var respectively as indicated in the schema definition
    """

    def __init__(self, adata: anndata.AnnData, pre_analysis: bool = False):
        """
        From a list of ids and defined constraints, creates a mapping dictionary {id: label, ...}

        :param str file_name: Path to h5ad file
        :param bool pre_analysis: Whether the dataset is pre-analysis
        """
        self.adata = adata
        self.pre_analysis = pre_analysis
        self.schema_version = schema.get_current_schema_version()
        self.schema_def = schema.get_schema_definition()
        self.errors = []

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

                if type(value_1) is not type(value_2):
                    raise ValueError("Inconsistent types, impossible to merge")

                if isinstance(value_2, str):
                    if key == "error_message_suffix":
                        merged_dict[key] = value_1 + " " + value_2
                        continue
                    if not value_2 == value_1:
                        raise ValueError(f"Strings types in dependencies cannot be different, {value_1} and {value_2}")

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

    def _get_ontology_term_label(
        self, term_id: str, allowed_ontologies: List[str], stripped_prefixes: Optional[dict] = None
    ) -> str:
        """
        Fetches human-readable label corresponding to a single ontology term_id, if it is a term_id defined in one of the
         allowed_ontologies. Raises ValueError if no ontology term label is found (should've triggered validation
         error if term_id is not valid).

        When stripped_prefixes is provided, terms matching a configured prefix (e.g. "anti-" for UniProt) are
        stripped before lookup and the prefix is prepended to the returned label
        (e.g. "anti-uniprot:Q99467" → "anti-" + label("uniprot:Q99467") = "anti-CD180_HUMAN").

        :param term_id: str single ontology term ID
        :param allowed_ontologies: List[str] list of ontologies to check for term_id label in
        :param stripped_prefixes: Optional[dict] prefix stripping rules from curie_constraints,
            e.g. {"UniProt": ["anti-"]}
        :return: str term label
        """
        label_prefix = ""
        term_to_look_up = term_id
        if stripped_prefixes:
            for _, prefixes in stripped_prefixes.items():
                for prefix in prefixes:
                    if term_id.startswith(prefix):
                        label_prefix = prefix
                        term_to_look_up = term_id[len(prefix) :]
                        break
                if label_prefix:
                    break

        for ontology_name in allowed_ontologies:
            if ontology_name == "NA":
                continue
            elif ONTOLOGY_PARSER.is_valid_term_id(term_to_look_up, ontology_name):
                return label_prefix + ONTOLOGY_PARSER.get_term_label(term_to_look_up)
        raise ValueError(f"Add labels error: Unable to get label for '{term_id}'")

    def _get_mapping_dict_curie(self, ids: List[str], curie_constraints: dict) -> Dict[str, str]:
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
        multi_term_def = curie_constraints.get("multi_term")
        delimiter = None if multi_term_def is None else multi_term_def["delimiter"]
        stripped_prefixes = curie_constraints.get("stripped_prefixes")

        # Map term_ids to their human-readable ontology labels
        for term_id in ids:
            # If there are exceptions the label should be the same as the id
            if "exceptions" in curie_constraints and term_id in curie_constraints["exceptions"]:
                mapping_dict[term_id] = term_id
                continue

            if delimiter is not None:
                labels = [
                    self._get_ontology_term_label(term, allowed_ontologies, stripped_prefixes)
                    for term in term_id.split(delimiter)
                ]
                mapping_dict[term_id] = delimiter.join(labels)
            else:
                mapping_dict[term_id] = self._get_ontology_term_label(term_id, allowed_ontologies, stripped_prefixes)

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
            organism = gencode.get_organism_from_feature_id(i)
            mapping_dict[i] = get_gene_checker(organism).get_symbol(i)

        return mapping_dict

    def _get_mapping_dict_feature_reference(self, ids: List[str]) -> Dict[str, Optional[gencode.SupportedOrganisms]]:
        """
        Creates a mapping dictionary of gene/feature IDs and NCBITaxon curies

        :param list[str] ids: Gene/feature IDs use for mapping

        :return a mapping dictionary: {id: label, ...}
        :rtype dict
        """

        mapping_dict = {}

        for i in ids:
            organism = gencode.get_organism_from_feature_id(i)
            mapping_dict[i] = organism.value

        return mapping_dict

    def _get_mapping_dict_feature_type(self, ids: List[str]) -> Dict[str, str]:
        """
        Creates a mapping dictionary of gene/feature IDs and its feature type.

        :param list[str] ids: Gene/feature IDs use for mapping

        :return a mapping dictionary: {id: feature_type, ...}
        :rtype dict
        """

        mapping_dict = {}

        for i in ids:
            organism = gencode.get_organism_from_feature_id(i)
            mapping_dict[i] = get_gene_checker(organism).get_type(i)

        return mapping_dict

    def _get_mapping_dict_feature_biotype(self, ids: List[str]) -> Dict[str, str]:
        """
        Creates a mapping dictionary of feature IDs and biotype ("gene" or "spike-in")

        :param list[str] ids: feature IDs use for mapping

        :return a mapping dictionary: {id: "gene", id: "spike-in", ...}
        :rtype dict
        """
        mapping_dict = {}

        for i in ids:
            if i.startswith("ERCC"):
                mapping_dict[i] = "spike-in"
            else:
                mapping_dict[i] = "gene"

        return mapping_dict

    def _get_mapping_dict_feature_length(self, ids: List[str]) -> Dict[str, int]:
        """
        Creates a mapping dictionary of feature IDs and feature length, fetching from pre-calculated gene info CSVs
        derived from GENCODE mappings for supported organisms. Set to 0 for non-gene features.

        :param list[str] ids: feature IDs use for mapping

        :return a mapping dictionary: {id: <int>, id: 0, ...}
        :rtype dict
        """
        mapping_dict = {}

        for i in ids:
            organism = gencode.get_organism_from_feature_id(i)
            mapping_dict[i] = get_gene_checker(organism).get_length(i)

        return mapping_dict

    def _get_labels(
        self,
        component: str,
        column: str,
        column_definition: dict,
        label_type: dict,
    ) -> pd.Series:
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
        current_df = getattr_anndata(self.adata, component)

        if column == "index":
            original_column = pd.Series(current_df.index)
            original_column.index = current_df.index
        else:
            original_column = getattr(current_df, column)

        ids = getattr(current_df, column).drop_duplicates().tolist()

        # Flatten column definition (will do so if there are dependencies in the definition
        column_definition = self._flatten_column_def_with_dependencies(column_definition)

        if label_type == "curie":
            if "curie_constraints" not in column_definition:
                raise ValueError(
                    f"Schema definition error: 'add_labels' with type 'curie' was found for '{column}' "
                    "but no curie constraints were found for the labels"
                )

            mapping_dict = self._get_mapping_dict_curie(ids, column_definition["curie_constraints"])

        elif label_type == "feature_id":
            mapping_dict = self._get_mapping_dict_feature_id(ids=ids)

        elif label_type == "feature_reference":
            mapping_dict = self._get_mapping_dict_feature_reference(ids=ids)

        elif label_type == "feature_biotype":
            mapping_dict = self._get_mapping_dict_feature_biotype(ids=ids)

        elif label_type == "feature_length":
            mapping_dict = self._get_mapping_dict_feature_length(ids=ids)

        elif label_type == "feature_type":
            mapping_dict = self._get_mapping_dict_feature_type(ids=ids)

        else:
            raise TypeError(f"'{label_type}' is not supported in 'add-labels' functionality")

        new_column = original_column.copy().map(mapping_dict).astype("category")

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
            new_column = self._get_labels(component, column, column_definition, label_def["type"])
            new_column_name = label_def["to_column"]
            getattr_anndata(self.adata, component)[new_column_name] = new_column

    def _add_labels(self):
        """
        Add columns to dataset dataframes based on values in other columns, as defined in schema definition yaml.
        """
        pre_analysis_constraints = self.schema_def.get("pre_analysis", {}) if self.pre_analysis else {}

        for component in ["obs", "var", "raw.var"]:
            # If the component does not exist, skip (this is for raw.var)
            if getattr_anndata(self.adata, component) is None:
                continue

            # Skip components not allowed in pre-analysis mode (e.g. obsm)
            comp_pa = pre_analysis_constraints.get(component, {})
            if not comp_pa.get("allowed", True):
                continue

            # Doing it for columns
            if "columns" in self.schema_def["components"][component]:
                component_df = getattr_anndata(self.adata, component)
                pa_forbidden_keys = {k for k, v in comp_pa.get("keys", {}).items() if not v.get("allowed", True)}
                for column, column_def in self.schema_def["components"][component]["columns"].items():
                    if "add_labels" in column_def:
                        # Skip optional columns that are not present in the dataframe
                        if not column_def.get("required", True) and column not in component_df.columns:
                            continue
                        # Skip columns explicitly forbidden in pre-analysis mode
                        if column in pa_forbidden_keys:
                            continue
                        self._add_column(component, column, column_def)

            # Doing it for index
            index_def = self.schema_def["components"][component]["index"]
            if "add_labels" in index_def:
                self._add_column(component, "index", index_def)

        uns_def = self.schema_def["components"]["uns"]
        uns_pa_forbidden_keys = {
            k for k, v in pre_analysis_constraints.get("uns", {}).get("keys", {}).items() if not v.get("allowed", True)
        }
        for key in uns_def["keys"]:
            key_def = uns_def["keys"][key]
            if "add_labels" in key_def:
                # Skip uns keys explicitly forbidden in pre-analysis mode (e.g. default_embedding)
                if key in uns_pa_forbidden_keys:
                    continue
                label_type = key_def["add_labels"][0]["type"]
                if label_type == "curie":
                    label_to_write = key_def["add_labels"][0]["to_key"]
                    term_id = self.adata.uns[key]
                    allowed_ontologies = key_def["curie_constraints"]["ontologies"]
                    self.adata.uns[label_to_write] = self._get_ontology_term_label(
                        term_id=term_id, allowed_ontologies=allowed_ontologies
                    )
                else:
                    raise TypeError(f"'{label_type}' is not supported with uns 'add-labels'")

    def _remove_categories_with_zero_values(self):
        df = self.adata.obs
        for column in df.columns:
            col = df[column]
            if col.dtype == "category":
                df[column] = col.cat.remove_unused_categories()

    def _get_perturbation_types_for_obs(self, ec_val: Optional[str], gp_val: Optional[str]) -> str:
        """
        Derive the perturbation_types value for a single observation.

        :param ec_val: value of experimental_condition_ontology_term_id, or None if absent
        :param gp_val: value of genetic_perturbation_id, or None if absent
        :return: "no perturbations" or sorted " || "-delimited perturbation type set
        """
        DIET_ONTO_TERM = "EFO:0002755"
        TEMPERATURE_ONTO_TERM = "EFO:0001702"

        ec_is_na = ec_val is None or ec_val == "na"
        gp_is_na = gp_val is None or gp_val == "na"

        if ec_is_na and gp_is_na:
            return "no perturbations"

        types: set = set()
        delimiter = " || "

        if not ec_is_na:
            for term in ec_val.split(delimiter):
                if term.startswith("CHEBI:"):
                    types.add("chemical")
                elif term.startswith("uniprot:") or term.startswith("anti-uniprot:"):
                    types.add("protein")
                elif term == TEMPERATURE_ONTO_TERM:
                    types.add("temperature")
                elif term.startswith("EFO:"):
                    ancestors = ONTOLOGY_PARSER.get_term_ancestors(term, include_self=True)
                    if DIET_ONTO_TERM in ancestors:
                        types.add("diet")

        if not gp_is_na:
            types.add("genetic")

        return delimiter.join(sorted(types)) if types else "no perturbations"

    def _annotate_perturbation_types(self):
        """
        Add obs['perturbation_types'] when experimental_condition_ontology_term_id or
        genetic_perturbation_id is present in obs.

        perturbation_types is a cross-column derivation (depends on both
        experimental_condition_ontology_term_id and genetic_perturbation_id) that cannot be expressed
        as a single-source add_labels entry, so it is handled here explicitly.
        """
        obs = self.adata.obs
        has_ec = "experimental_condition_ontology_term_id" in obs.columns
        has_gp = "genetic_perturbation_id" in obs.columns

        if not has_ec and not has_gp:
            return

        ec_vals = obs["experimental_condition_ontology_term_id"].astype(str) if has_ec else None
        gp_vals = obs["genetic_perturbation_id"].astype(str) if has_gp else None

        # Only compute once per unique (ec_val, gp_val) pair for efficiency
        iter_ec = ec_vals if ec_vals is not None else (["na"] * len(obs))
        iter_gp = gp_vals if gp_vals is not None else (["na"] * len(obs))
        unique_combos = set(zip(iter_ec, iter_gp))
        combo_to_type = {(ec, gp): self._get_perturbation_types_for_obs(ec, gp) for ec, gp in unique_combos}

        pt_values = [combo_to_type[(ec, gp)] for ec, gp in zip(iter_ec, iter_gp)]
        obs["perturbation_types"] = pd.Categorical(pt_values)

    def _annotate_genetic_perturbations(self):
        """
        Annotate genetic perturbations with genomic locations and target genes if genetic_perturbations exists in uns.
        """
        try:
            self.adata = annotate_perturbations_in_h5ad(self.adata, index_path=None)
            logger.info("Genetic perturbations annotation completed")
        except RuntimeError as e:
            # guidescan2 not installed or other runtime errors
            error_msg = f"Genetic perturbations annotation skipped: {e}"
            logger.warning(error_msg)
            self.errors.append((error_msg, None))
        except Exception as e:
            # Other errors (organism not supported, etc.)
            error_msg = f"Genetic perturbations annotation failed: {e}"
            logger.warning(error_msg)
            import traceback

            tb = traceback.format_exc()
            self.errors.append((error_msg, tb))

    def _build_schema_reference_url(self, schema_version: str):
        return f"{SCHEMA_REFERENCE_BASE_URL}/{schema_version}/{SCHEMA_REFERENCE_FILE_NAME}"

    def write_labels(self, add_labels_file: str):
        """
        From a valid (per cellxgene's schema) h5ad, this function writes a new h5ad file with ontology/gene labels added
        to adata.obs  and adata.var respectively

        :param str add_labels_file: Path to new h5ad file with ontology/gene labels added

        :rtype None
        """
        logger.info("Writing labels")

        # Add columns to dataset dataframes based on values in other columns, as defined in schema definition yaml
        self._add_labels()

        # Annotate perturbation_types (cross-column derivation from ec and gp fields)
        self._annotate_perturbation_types()

        # Annotate genetic perturbations if present
        self._annotate_genetic_perturbations()

        # Remove unused categories
        self._remove_categories_with_zero_values()

        # Annotate Reserved Columns

        self.adata.uns["schema_version"] = self.schema_version
        self.adata.uns["schema_reference"] = self._build_schema_reference_url(self.schema_version)
        self.adata.uns["is_pre_analysis"] = self.pre_analysis
        self.adata.obs["observation_joinid"] = get_hash_digest_column(self.adata.obs)
        logger.info(f"Labels have been added. Writing to {add_labels_file}")
        # Write file
        try:
            self.adata.write_h5ad(add_labels_file, compression="gzip")
        except Exception as e:
            tb = traceback.format_exc()
            self.errors.append((f"Writing h5ad was unsuccessful, got exception '{e}'.", tb))

        # Print errors if any
        if self.errors:
            for e, tb in self.errors:
                logger.error(e, extra={"exec_info": tb})
            return False
        else:
            return True
