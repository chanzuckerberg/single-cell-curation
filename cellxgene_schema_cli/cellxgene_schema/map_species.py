import logging

import anndata as ad
from cellxgene_ontology_guide.ontology_parser import OntologyParser

logger = logging.getLogger(__name__)


def map_species(input_file, output_file):
    ontology_parser = OntologyParser("6.0.0")
    adata = ad.read_h5ad(input_file, backed="r")
    map_columns = [
        ("organism_cell_type_ontology_term_id", "cell_type_ontology_term_id", "CL"),
        ("organism_tissue_ontology_term_id", "tissue_ontology_term_id", "UBERON"),
    ]
    for map_from, map_to, ontology in map_columns:
        logger.info(f"Mapping {ontology} terms to {map_to}")
        if not hasattr(adata.obs, map_to):
            logger.info(f"No existing {map_to} column detected, creating one and defaulting to empty strings.")
            adata.obs[map_to] = ""
        for category, group in adata.obs.groupby(map_from):
            try:
                closest_terms = ontology_parser.get_closest_bridge_term_ids(category, ontology)
            except Exception:
                logger.error(f"{category} in {map_from} is not a valid ontology term ID, cannot find {ontology} map.")
                continue
            if len(closest_terms) < 1:
                logger.warning(f"{category} has no closest match for {ontology} in its ancestry tree.")
            elif len(closest_terms) == 1:
                adata.obs.loc[group.index, map_to] = closest_terms[0]
                logger.info(f"Setting {category} {map_to} rows to single closest match: {closest_terms[0]}")
            else:
                logger.warning(
                    f"{category} has multiple closest matches for {ontology} in its ancestry tree: "
                    f"{closest_terms}. Pick one based on biological context and assign manually."
                )
        adata.obs[map_to] = adata.obs[map_to].astype("category")

    logger.info(f"Mappings complete, writing to {output_file}")
    adata.write_h5ad(output_file, compression="gzip")
