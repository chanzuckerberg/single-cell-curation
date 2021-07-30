import pandas as pd
import anndata
import os

SCHEMA_VERSION = "2.0.0"
FIXTURES_ROOT = os.path.join(os.path.dirname(__file__))

# Pre-made examples
h5ad_dir = os.path.join(FIXTURES_ROOT, "h5ads")
h5ad_valid = os.path.join(h5ad_dir, "example_valid.h5ad")
h5ad_invalid = [
    "example_invalid_CL.h5ad",
    "example_invalid_assay.h5ad",
    "example_invalid_disease.h5ad",
    "example_invalid_organism.h5ad",
    "example_invalid_primary_data.h5ad",
    "example_invalid_sex.h5ad",
    "example_invalid_tissue.h5ad",
    "example_invalid_ethnicity.h5ad",
    "example_invalid_development_stage.h5ad",
]
h5ad_invalid = [os.path.join(h5ad_dir, i) for i in h5ad_invalid]

# Manual minimal examples
good_obs = pd.DataFrame(
    [
        [
            "CL:0000066",
            "EFO:0009899",
            "MONDO:0100096",
            "NCBITaxon:9606",
            "PATO:0000383",
            "UBERON:0002048",
            True,
            "HANCESTRO:0575",
            "HsapDv:0000003",
        ],
        [
            "CL:0000192",
            "EFO:0010183 (sci-plex)",
            "PATO:0000461",
            "NCBITaxon:10090",
            "unknown",
            "CL:0000192 (cell culture)",
            False,
            "na",
            "MmusDv:0000003",
        ],
    ],
    index=["X", "Y"],
    columns=[
        "cell_type_ontology_term_id",
        "assay_ontology_term_id",
        "disease_ontology_term_id",
        "organism_ontology_term_id",
        "sex_ontology_term_id",
        "tissue_ontology_term_id",
        "is_primary_data",
        "ethnicity_ontology_term_id",
        "development_stage_ontology_term_id",
    ],
)

obs_expected = pd.DataFrame(
    [
        [
            "epithelial cell",
            "10x 3' v2",
            "COVID-19",
            "Homo sapiens",
            "female",
            "lung",
            "Yoruban",
            "Carnegie stage 01",
        ],
        [
            "smooth muscle cell",
            "single cell library construction (sci-plex)",
            "normal",
            "Mus musculus",
            "unknown",
            "smooth muscle cell (cell culture)",
            "na",
            "Theiler stage 01",
        ],
    ],
    index=["X", "Y"],
    columns=[
        "cell_type",
        "assay",
        "disease",
        "organism",
        "sex",
        "tissue",
        "ethnicity",
        "development_stage",
    ],
)

bad_obs = pd.DataFrame(
    [
        [
            "CL:NO_TERM",
            "EFO:NO_TERM",
            "MONDO:NO_TERM",
            "NCBITaxon:NO_TERM",
            "PATO:NO_TERM",
            "UBERON:NO_TERM",
            "True",
        ],
        [
            "CL:0000182",
            "EFO:00212",
            "MONDO:0324",
            "NCBITaxon:00324",
            "PATO:2003",
            "UBERON:3203",
            "False",
        ],
    ],
    index=["X", "Y"],
    columns=[
        "cell_type_ontology_term_id",
        "assay_ontology_term_id",
        "disease_ontology_term_id",
        "organism_ontology_term_id",
        "sex_ontology_term_id",
        "tissue_ontology_term_id",
        "is_primary_data",
    ],
)

good_uns = {"schema_version": SCHEMA_VERSION}

X = pd.DataFrame(
    [[0] * good_obs.shape[1], [0] * good_obs.shape[1]],
    index=["X", "Y"],
)

adata = anndata.AnnData(X=X, obs=good_obs, uns=good_uns)
adata_with_labels = anndata.AnnData(
    X=X, obs=pd.concat([good_obs, obs_expected], axis=1), uns=good_uns
)
adata_empty = anndata.AnnData(X=X, uns=good_uns)
