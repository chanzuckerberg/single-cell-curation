# flake8: noqa
import pandas as pd
import numpy
import anndata
import os
from scipy import sparse
from cellxgene_schema.utils import get_hash_digest_column

# -----------------------------------------------------------------#
# General example information
FIXTURES_ROOT = os.path.join(os.path.dirname(__file__))

# -----------------------------------------------------------------#
# Pre-made example files
h5ad_dir = os.path.join(FIXTURES_ROOT, "h5ads")
h5ad_valid = os.path.join(h5ad_dir, "example_valid.h5ad")
h5ad_invalid = os.path.join(h5ad_dir, "example_invalid_CL.h5ad")

# -----------------------------------------------------------------#
# Manually creating minimal anndata objects.
#
# The valid objects mentioned below contain all valid cases covered in the schema, including multiple examples for
# fields that allow multiple valid options.
#
# This process entails:
# 1. Creating individual obs components: one valid dataframe, and one with labels (extra columns that are supposed
#   to be added by validator)
# 2. Creating individual var components: valid, and one with labels
# 3. Creating individual uns valid component
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
            "tissue",
            True,
            "HANCESTRO:0575",
            "HsapDv:0000003",
            "donor_1",
            "nucleus",
        ],
        [
            "CL:0000192",
            "EFO:0009918",
            "PATO:0000461",
            "NCBITaxon:10090",
            "unknown",
            "CL:0000192",
            "cell culture",
            False,
            "na",
            "MmusDv:0000003",
            "donor_2",
            "na",
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
        "tissue_type",
        "is_primary_data",
        "self_reported_ethnicity_ontology_term_id",
        "development_stage_ontology_term_id",
        "donor_id",
        "suspension_type",
    ],
)

good_obs.loc[:, ["donor_id"]] = good_obs.astype("category")
good_obs.loc[:, ["suspension_type"]] = good_obs.astype("category")
good_obs.loc[:, ["tissue_type"]] = good_obs.astype("category")

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
            "smFISH",
            "normal",
            "Mus musculus",
            "unknown",
            "smooth muscle cell",
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
        "self_reported_ethnicity",
        "development_stage",
    ],
)

obs_expected["observation_joinid"] = get_hash_digest_column(obs_expected)

# Valid spatial obs per schema
good_obs_visium = pd.DataFrame(
    [
        [
            1,
            1,
            "unknown",
            "EFO:0010961",
            "MONDO:0100096",
            "NCBITaxon:9606",
            "PATO:0000383",
            "UBERON:0002048",
            "tissue",
            True,
            "HANCESTRO:0575",
            "HsapDv:0000003",
            "donor_1",
            "na",
            0,
        ],
        [
            2,
            2,
            "CL:0000192",
            "EFO:0010961",
            "PATO:0000461",
            "NCBITaxon:10090",
            "unknown",
            "CL:0000192",
            "cell culture",
            False,
            "na",
            "MmusDv:0000003",
            "donor_2",
            "na",
            1,
        ],
    ],
    index=["X", "Y"],
    columns=[
        "array_col",
        "array_row",
        "cell_type_ontology_term_id",
        "assay_ontology_term_id",
        "disease_ontology_term_id",
        "organism_ontology_term_id",
        "sex_ontology_term_id",
        "tissue_ontology_term_id",
        "tissue_type",
        "is_primary_data",
        "self_reported_ethnicity_ontology_term_id",
        "development_stage_ontology_term_id",
        "donor_id",
        "suspension_type",
        "in_tissue",
    ],
)

good_obs_visium.loc[:, ["donor_id"]] = good_obs_visium.astype("category")
good_obs_visium.loc[:, ["suspension_type"]] = good_obs_visium.astype("category")
good_obs_visium.loc[:, ["tissue_type"]] = good_obs_visium.astype("category")

# Valid spatial obs per schema
good_obs_slide_seqv2 = pd.DataFrame(
    [
        [
            "CL:0000066",
            "EFO:0030062",
            "MONDO:0100096",
            "NCBITaxon:9606",
            "PATO:0000383",
            "UBERON:0002048",
            "tissue",
            True,
            "HANCESTRO:0575",
            "HsapDv:0000003",
            "donor_1",
            "na",
        ],
        [
            "CL:0000192",
            "EFO:0030062",
            "PATO:0000461",
            "NCBITaxon:10090",
            "unknown",
            "CL:0000192",
            "cell culture",
            False,
            "na",
            "MmusDv:0000003",
            "donor_2",
            "na",
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
        "tissue_type",
        "is_primary_data",
        "self_reported_ethnicity_ontology_term_id",
        "development_stage_ontology_term_id",
        "donor_id",
        "suspension_type",
    ],
)

good_obs_slide_seqv2.loc[:, ["donor_id"]] = good_obs_slide_seqv2.astype("category")
good_obs_slide_seqv2.loc[:, ["suspension_type"]] = good_obs_slide_seqv2.astype("category")
good_obs_slide_seqv2.loc[:, ["tissue_type"]] = good_obs.astype("category")

good_obs_visium_is_single_false = pd.DataFrame(
    [
        [
            "CL:0000066",
            "EFO:0010961",
            "MONDO:0100096",
            "NCBITaxon:9606",
            "PATO:0000383",
            "UBERON:0002048",
            "tissue",
            False,
            "HANCESTRO:0575",
            "HsapDv:0000003",
            "donor_1",
            "na",
        ],
        [
            "CL:0000192",
            "EFO:0010961",
            "PATO:0000461",
            "NCBITaxon:10090",
            "unknown",
            "CL:0000192",
            "cell culture",
            False,
            "na",
            "MmusDv:0000003",
            "donor_2",
            "na",
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
        "tissue_type",
        "is_primary_data",
        "self_reported_ethnicity_ontology_term_id",
        "development_stage_ontology_term_id",
        "donor_id",
        "suspension_type",
    ],
)

good_obs_visium_is_single_false.loc[:, ["donor_id"]] = good_obs_visium_is_single_false.astype("category")
good_obs_visium_is_single_false.loc[:, ["suspension_type"]] = good_obs_visium_is_single_false.astype("category")
good_obs_visium_is_single_false.loc[:, ["tissue_type"]] = good_obs_visium_is_single_false.astype("category")

# ---
# 2. Creating individual var components: valid object and valid object and with labels

# Valid var per schema
good_var = pd.DataFrame(
    [
        [False],
        [False],
        [False],
        [False],
    ],
    index=["ERCC-00002", "ENSG00000127603", "ENSMUSG00000059552", "ENSSASG00005000004"],
    columns=["feature_is_filtered"],
)

# Expected var, this is what the obs above should look like after adding the necessary columns with the validator,
# these columns are defined in the schema
var_expected = pd.DataFrame(
    [
        ["spike-in", False, "ERCC-00002 (spike-in control)", "NCBITaxon:32630", 0],
        ["gene", False, "MACF1", "NCBITaxon:9606", 70573],
        ["gene", False, "Trp53", "NCBITaxon:10090", 4045],
        ["gene", False, "S", "NCBITaxon:2697049", 3822],
    ],
    index=["ERCC-00002", "ENSG00000127603", "ENSMUSG00000059552", "ENSSASG00005000004"],
    columns=[
        "feature_biotype",
        "feature_is_filtered",
        "feature_name",
        "feature_reference",
        "feature_length",
    ],
)

# ---
# 3. Creating individual uns component
good_uns = {
    "title": "A title",
    "default_embedding": "X_umap",
    "X_approximate_distribution": "normal",
    "batch_condition": ["is_primary_data"],
}

good_uns_with_labels = {
    "schema_version": "4.0.0",
    "schema_reference": "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/4.0.0/schema.md",
    "citation": "Publication: <doi> Dataset Version: "
    "https://datasets.cellxgene.cziscience.com/<dataset_version_id>.h5ad curated and distributed by CZ "
    "CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/<collection_id>",
    "title": "A title",
    "default_embedding": "X_umap",
    "X_approximate_distribution": "normal",
    "batch_condition": ["is_primary_data"],
}

good_uns_with_colors = {
    "title": "A title",
    "default_embedding": "X_umap",
    "X_approximate_distribution": "normal",
    "batch_condition": ["is_primary_data"],
    "suspension_type_colors": numpy.array(["red", "blue"]),
    "donor_id_colors": numpy.array(["#000000", "#ffffff"]),
    "tissue_type_colors": numpy.array(["black", "pink"]),
}

visium_library_id = "Proj2023_Lung_C001"

good_uns_with_visium_spatial = {
    "title": "A title",
    "default_embedding": "X_umap",
    "X_approximate_distribution": "normal",
    "batch_condition": ["is_primary_data"],
    "spatial": {
        "is_single": numpy.bool_(True),
        visium_library_id: {
            "images": {
                "hires": numpy.zeros((2000, 100, 3), dtype=numpy.uint8),
                "fullres": numpy.array([[[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]]], dtype=numpy.uint8),
            },
            "scalefactors": {
                "spot_diameter_fullres": 1.0,
                "tissue_hires_scalef": 1.0,
            },
        },
    },
}

good_uns_with_is_single_false = {
    "title": "A title",
    "default_embedding": "X_umap",
    "X_approximate_distribution": "normal",
    "batch_condition": ["is_primary_data"],
    "spatial": {"is_single": False},
}

good_uns_with_slide_seqV2_spatial = {
    "title": "A title",
    "default_embedding": "X_umap",
    "X_approximate_distribution": "normal",
    "batch_condition": ["is_primary_data"],
    "spatial": {"is_single": True},
}

# ---
# 4. Creating expression matrix,
# X has integer values and non_raw_X has real values
X = numpy.ones([good_obs.shape[0], good_var.shape[0]], dtype=numpy.float32)
non_raw_X = sparse.csr_matrix(X.copy())
non_raw_X[0, 0] = 1.5

# ---
# 5.Creating valid obsm
good_obsm = {"X_umap": numpy.zeros([X.shape[0], 2])}
good_obsm_spatial = {"X_umap": numpy.zeros([X.shape[0], 2]), "spatial": numpy.zeros([X.shape[0], 2])}

# ---
# 6. Putting all the components created in the previous steps into minimal anndata that used for testing in
#   the unittests

# Valid anndata
adata = anndata.AnnData(X=sparse.csr_matrix(X), obs=good_obs, uns=good_uns, obsm=good_obsm, var=good_var)
adata.raw = adata.copy()
adata.X = non_raw_X
adata.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# Anndata with "X" and "raw.X" but neither has actual raw values
adata_no_raw_values = anndata.AnnData(
    X=sparse.csr_matrix(non_raw_X),
    obs=good_obs,
    uns=good_uns,
    obsm=good_obsm,
    var=good_var,
)
adata_no_raw_values.raw = adata_no_raw_values.copy()
adata_no_raw_values.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# Anndata with no obs nor var
adata_minimal = anndata.AnnData(X=sparse.csr_matrix(X), uns=good_uns, obsm=good_obsm)

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
    var=var_expected,
    uns=good_uns_with_labels,
    obsm=good_obsm,
)

# Expected anndata with colors for categorical obs fields
adata_with_colors = anndata.AnnData(
    X=sparse.csr_matrix(X), obs=good_obs, uns=good_uns_with_colors, obsm=good_obsm, var=good_var
)

# Expected anndata with Visium spatial data
adata_visium = anndata.AnnData(
    X=sparse.csr_matrix(X), obs=good_obs_visium, uns=good_uns_with_visium_spatial, obsm=good_obsm_spatial, var=good_var
)
adata_visium.raw = adata_visium.copy()
adata_visium.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# Expected anndata with Slide-seqV2 spatial data
adata_slide_seqv2 = anndata.AnnData(
    X=sparse.csr_matrix(X),
    obs=good_obs_slide_seqv2,
    uns=good_uns_with_slide_seqV2_spatial,
    obsm=good_obsm_spatial,
    var=good_var,
)

adata_spatial_is_single_false = anndata.AnnData(
    X=sparse.csr_matrix(X),
    obs=good_obs_visium_is_single_false,
    uns=good_uns_with_is_single_false,
    obsm=good_obsm_spatial,
    var=good_var,
)

# anndata for testing migration
unmigrated_obs = pd.DataFrame(
    [
        [
            "cell_type:1",
            "assay:1",
            "disease:1",
            "organism:1",
            "sex:1",
            "tissue:1",
            "sre:1",
            "development_stage:1",
        ],
        [
            "cell_type:1",
            "assay:1",
            "disease:1",
            "organism:1",
            "sex:1",
            "tissue:1",
            "sre:1",
            "development_stage:1",
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
        "self_reported_ethnicity_ontology_term_id",
        "development_stage_ontology_term_id",
    ],
)

var_unmigrated = pd.DataFrame(
    [
        [False],
        [False],
    ],
    index=["ENSSASG00005000004", "DUMMY"],
    columns=[
        "feature_is_filtered",
    ],
)

unmigrated_X = numpy.zeros([unmigrated_obs.shape[0], var_unmigrated.shape[0]], dtype=numpy.float32)
adata_with_labels_unmigrated = anndata.AnnData(
    X=sparse.csr_matrix(unmigrated_X),
    obs=unmigrated_obs,
    uns=good_uns_with_labels,
    var=var_unmigrated,
    obsm={"X_umap": numpy.zeros([unmigrated_X.shape[0], 2])},
)
adata_with_labels_unmigrated.raw = adata_with_labels_unmigrated.copy()
