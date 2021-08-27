# flake8: noqa
import pandas as pd
import numpy
import anndata
import os
from scipy import sparse

# -----------------------------------------------------------------#
# General example information
SCHEMA_VERSION = "2.0.0"
FIXTURES_ROOT = os.path.join(os.path.dirname(__file__))

# -----------------------------------------------------------------#
# Pre-made example files
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
    "example_invalid_uns.h5ad",
]
h5ad_invalid = [os.path.join(h5ad_dir, i) for i in h5ad_invalid]


# -----------------------------------------------------------------#
# Manually creating minimal anndata objects.
#
# This process entails:
# 1. Creating individual obs components: valid, invalid, and with labels (extra columns that are supposed
#   to be added by validator)
# 2. Creating individual var components: valid, invalid, and with labels
# 3. Creating individual uns components: valid and invalid
# 4. Creating expression matrices
# 5. Creating valid obsm
# 6. Putting all the components created in the previous steps into minimal anndata that used for testing in
#   the unittests

# Valid obs per schema
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

# Expected obs, this is what the obs above should look like after adding the necessary columns with the validator,
# these columns are defined in the schema
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

# invalid obs, all fields are designed to fail
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

# ---
# 2. Creating individual var components: valid, invalid, and with labels

# Valid var per schema
good_var = pd.DataFrame(
    [
        ["spike-in", False],
        ["gene", False],
        ["gene", False],
        ["gene", False],
    ],
    index=["ERCC-00002", "ENSG00000127603", "ENSMUSG00000059552", "ENSSASG00005000004"],
    columns=["feature_biotype", "feature_is_filtered"],
)
good_var.loc[:, ["feature_biotype"]] = good_var.astype("category")

# Expected var, this is what the obs above should look like after adding the necessary columns with the validator,
# these columns are defined in the schema
var_expected = pd.DataFrame(
    [
        ["spike-in", False, "ERCC-00002 spike-in control", "NCBITaxon:32630"],
        ["gene", False, "MACF1", "NCBITaxon:9606"],
        ["gene", False, "Trp53", "NCBITaxon:10090"],
        ["gene", False, "S", "NCBITaxon:2697049"],
    ],
    index=["ERCC-00002", "ENSG00000127603", "ENSMUSG00000059552", "ENSSASG00005000004"],
    columns=[
        "feature_biotype",
        "feature_is_filtered",
        "feature_name",
        "feature_reference",
    ],
)
var_expected.loc[:, ["feature_biotype"]] = var_expected.astype("category")

# invalid var, all fields are designed to fail
bad_var = pd.DataFrame(
    [
        ["gene", "False"],  # should be spike in
        ["spike-in", "False"],  # should be gene
        ["other", "False"],  # incorrect
        ["gene", "False"],
    ],
    index=["ERCC-00002", "ENSG00000127603", "ENSMUSG00000059552", "NO_GENE"],
    columns=["feature_biotype", "feature_is_filtered"],
)
bad_var.loc[:, ["feature_biotype"]] = bad_var.astype("category")


# ---
# 3. Creating individual uns components: valid and invalid
good_uns = {
    "schema_version": SCHEMA_VERSION,
    "title": "A title",
    "default_embedding": "X_umap",
    "X_normalization": "none",
    "X_approximate_distribution": "normal",
    "batch_condition": ["is_primary_data"],
}

bad_uns = {
    "schema_version": "2.0.0",
    "title": 1,
    "default_embedding": "X_PCA_1",
    "X_normalization": 1,
    "X_approximate_distribution": "CPM",
    "batch_condition": ["batchD", "batchE"],
}


# ---
# 4. Creating expression matrix,
# X has integer values and non_raw_X has real values
X = numpy.zeros([good_obs.shape[0], good_var.shape[0]])
non_raw_X = X.copy()
non_raw_X[0, 0] = 1.5

# ---
# 5.Creating valid obsm
good_obsm = {"X_umap": numpy.zeros([X.shape[0], 2])}


# ---
# 6. Putting all the components created in the previous steps into minimal anndata that used for testing in
#   the unittests

# Valid anndata
adata = anndata.AnnData(
    X=sparse.csr_matrix(X), obs=good_obs, uns=good_uns, obsm=good_obsm, var=good_var
)

# Anndata with no obs nor var
adata_empty = anndata.AnnData(X=sparse.csr_matrix(X), uns=good_uns, obsm=good_obsm)

# Anndata with a expression matrix that is not raw
adata_non_raw = anndata.AnnData(
    X=sparse.csr_matrix(non_raw_X),
    obs=good_obs,
    uns=good_uns,
    obsm=good_obsm,
    var=good_var,
)

# Expected anndata with labels that the validator must write in obs and var
adata_with_labels = anndata.AnnData(
    X=sparse.csr_matrix(X),
    obs=pd.concat([good_obs, obs_expected], axis=1),
    var=pd.concat([good_var, var_expected], axis=1),
    uns=good_uns,
    obsm=good_obsm,
)
