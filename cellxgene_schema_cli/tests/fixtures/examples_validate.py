# flake8: noqa
import pandas as pd
import numpy
import anndata
import os
from scipy import sparse
from cellxgene_schema.utils import get_hash_digest_column
from dask.array import from_array

# -----------------------------------------------------------------#
# General example information
FIXTURES_ROOT = os.path.join(os.path.dirname(__file__))

# -----------------------------------------------------------------#
# Pre-made example files
h5ad_dir = os.path.join(FIXTURES_ROOT, "h5ads")
h5ad_valid = os.path.join(h5ad_dir, "example_valid.h5ad")
h5ad_invalid = os.path.join(h5ad_dir, "example_invalid_CL.h5ad")
h5ad_perturbations = os.path.join(h5ad_dir, "small_perturbations.h5ad")

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
            "PATO:0000383",
            "UBERON:0002048",
            "tissue",
            True,
            "HANCESTRO:0019",
            "HsapDv:0000003",
            "donor_1",
            "nucleus",
        ],
        [
            "CL:0000066",
            "EFO:0009899",
            "MONDO:0100096",
            "PATO:0000383",
            "UBERON:0002048",
            "tissue",
            True,
            "HANCESTRO:0019",
            "HsapDv:0000003",
            "donor_1",
            "nucleus",
        ],
    ],
    index=["X", "Y"],
    columns=[
        "cell_type_ontology_term_id",
        "assay_ontology_term_id",
        "disease_ontology_term_id",
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

good_obs["donor_id"] = good_obs["donor_id"].astype("category")
good_obs["suspension_type"] = good_obs["suspension_type"].astype("category")
good_obs["tissue_type"] = good_obs["tissue_type"].astype("category")
good_obs["tissue_type"] = good_obs["tissue_type"].cat.add_categories(["primary cell culture", "organoid", "cell line"])

# Expected obs, this is what the obs above should look like after adding the necessary columns with the validator,
# these columns are defined in the schema
obs_expected = pd.DataFrame(
    [
        [
            "epithelial cell",
            "10x 3' v2",
            "COVID-19",
            "female",
            "lung",
            "Japanese",
            "Carnegie stage 01",
        ],
        [
            "epithelial cell",
            "10x 3' v2",
            "COVID-19",
            "female",
            "lung",
            "Japanese",
            "Carnegie stage 01",
        ],
    ],
    index=["X", "Y"],
    columns=[
        "cell_type",
        "assay",
        "disease",
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
            "EFO:0022859",
            "MONDO:0100096",
            "PATO:0000383",
            "UBERON:0002048",
            "tissue",
            True,
            "HANCESTRO:0019",
            "HsapDv:0000003",
            "donor_1",
            "na",
            0,
        ],
        [
            1,
            1,
            "unknown",
            "EFO:0022859",
            "MONDO:0100096",
            "PATO:0000383",
            "UBERON:0002048",
            "tissue",
            True,
            "HANCESTRO:0019",
            "HsapDv:0000003",
            "donor_1",
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

good_obs_visium["donor_id"] = good_obs_visium["donor_id"].astype("category")
good_obs_visium["suspension_type"] = good_obs_visium["suspension_type"].astype("category")
good_obs_visium["tissue_type"] = good_obs_visium["tissue_type"].astype("category")
good_obs_visium["tissue_type"] = good_obs_visium["tissue_type"].cat.add_categories(
    ["primary cell culture", "organoid", "cell line"]
)

# Valid spatial obs per schema
good_obs_slide_seqv2 = pd.DataFrame(
    [
        [
            "CL:0000066",
            "EFO:0030062",
            "MONDO:0100096",
            "PATO:0000383",
            "UBERON:0002048",
            "tissue",
            True,
            "HANCESTRO:0019",
            "HsapDv:0000003",
            "donor_1",
            "na",
        ],
        [
            "CL:0000066",
            "EFO:0030062",
            "MONDO:0100096",
            "PATO:0000383",
            "UBERON:0002048",
            "tissue",
            True,
            "HANCESTRO:0019",
            "HsapDv:0000003",
            "donor_1",
            "na",
        ],
    ],
    index=["X", "Y"],
    columns=[
        "cell_type_ontology_term_id",
        "assay_ontology_term_id",
        "disease_ontology_term_id",
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

good_obs_slide_seqv2["donor_id"] = good_obs_slide_seqv2["donor_id"].astype("category")
good_obs_slide_seqv2["suspension_type"] = good_obs_slide_seqv2["suspension_type"].astype("category")
good_obs_slide_seqv2["tissue_type"] = good_obs_slide_seqv2["tissue_type"].astype("category")
good_obs_slide_seqv2["tissue_type"] = good_obs_slide_seqv2["tissue_type"].cat.add_categories(
    ["primary cell culture", "organoid", "cell line"]
)

good_obs_visium_is_single_false = pd.DataFrame(
    [
        [
            "CL:0000066",
            "EFO:0022859",
            "MONDO:0100096",
            "PATO:0000383",
            "UBERON:0002048",
            "tissue",
            False,
            "HANCESTRO:0019",
            "HsapDv:0000003",
            "donor_1",
            "na",
        ],
        [
            "CL:0000066",
            "EFO:0022859",
            "MONDO:0100096",
            "PATO:0000383",
            "UBERON:0002048",
            "tissue",
            False,
            "HANCESTRO:0019",
            "HsapDv:0000003",
            "donor_1",
            "na",
        ],
    ],
    index=["X", "Y"],
    columns=[
        "cell_type_ontology_term_id",
        "assay_ontology_term_id",
        "disease_ontology_term_id",
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

good_obs_visium_is_single_false["donor_id"] = good_obs_visium_is_single_false["donor_id"].astype("category")
good_obs_visium_is_single_false["suspension_type"] = good_obs_visium_is_single_false["suspension_type"].astype(
    "category"
)
good_obs_visium_is_single_false["tissue_type"] = good_obs_visium_is_single_false["tissue_type"].astype("category")
good_obs_visium_is_single_false["tissue_type"] = good_obs_visium_is_single_false["tissue_type"].cat.add_categories(
    ["primary cell culture", "organoid", "cell line"]
)

good_obs_mouse = pd.DataFrame(
    [
        [
            "CL:0000192",
            "EFO:0008992",
            "PATO:0000461",
            "unknown",
            "CL:0000192",
            "primary cell culture",
            False,
            "na",
            "MmusDv:0000003",
            "donor_2",
            "na",
        ],
        [
            "CL:0000192",
            "EFO:0008992",
            "PATO:0000461",
            "unknown",
            "CL:0000192",
            "primary cell culture",
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

good_obs_mouse["donor_id"] = good_obs_mouse["donor_id"].astype("category")
good_obs_mouse["suspension_type"] = good_obs_mouse["suspension_type"].astype("category")
good_obs_mouse["tissue_type"] = good_obs_mouse["tissue_type"].astype("category")
good_obs_mouse["tissue_type"] = good_obs_mouse["tissue_type"].cat.add_categories(["tissue", "organoid", "cell line"])

# Creating a cell line obs by copying good_obs and changing the necessary fields
good_obs_cell_line = good_obs.copy()
good_obs_cell_line.loc[:, "tissue_type"] = "cell line"
good_obs_cell_line.loc[:, "tissue_ontology_term_id"] = "CVCL_0001"
good_obs_cell_line.loc[:, "development_stage_ontology_term_id"] = "na"
good_obs_cell_line.loc[:, "sex_ontology_term_id"] = "na"
good_obs_cell_line.loc[:, "self_reported_ethnicity_ontology_term_id"] = "na"
good_obs_cell_line.loc[:, "donor_id"] = "na"
good_obs_cell_line["donor_id"] = good_obs_cell_line["donor_id"].astype("category")

# ---
# 2. Creating individual var components: valid object and valid object and with labels

# Valid var per schema
good_var = pd.DataFrame(
    [[False]],
    index=[
        "ENSG00000127603",
        "ENSG00000141510",
        "ENSG00000012048",
        "ENSG00000139618",
        "ENSG00000002330",
        "ENSG00000000005",
        "ENSG00000000419",
    ],
    columns=["feature_is_filtered"],
)
good_var_mouse = pd.DataFrame(
    [[False]],
    index=[
        "ENSMUSG00000102693",
        "ENSMUSG00000064842",
        "ENSMUSG00000051951",
        "ENSMUSG00000102851",
        "ENSMUSG00000103377",
        "ENSMUSG00000104017",
        "ENSMUSG00000103025",
    ],
    columns=["feature_is_filtered"],
)

# Expected var, this is what the obs above should look like after adding the necessary columns with the validator,
# these columns are defined in the schema
var_expected = pd.DataFrame(
    [
        ["gene", False, "MACF1", "NCBITaxon:9606", 2821, "protein_coding"],
        ["gene", False, "TP53", "NCBITaxon:9606", 2426, "protein_coding"],
        ["gene", False, "BRCA1", "NCBITaxon:9606", 3757, "protein_coding"],
        ["gene", False, "BRCA2", "NCBITaxon:9606", 11428, "protein_coding"],
        ["gene", False, "BAD", "NCBITaxon:9606", 552, "protein_coding"],
        ["gene", False, "TNMD", "NCBITaxon:9606", 873, "protein_coding"],
        ["gene", False, "DPM1", "NCBITaxon:9606", 1262, "protein_coding"],
    ],
    index=[
        "ENSG00000127603",
        "ENSG00000141510",
        "ENSG00000012048",
        "ENSG00000139618",
        "ENSG00000002330",
        "ENSG00000000005",
        "ENSG00000000419",
    ],
    columns=[
        "feature_biotype",
        "feature_is_filtered",
        "feature_name",
        "feature_reference",
        "feature_length",
        "feature_type",
    ],
)

NUMBER_OF_GENES = len(var_expected.index)

# ---
# 3. Creating individual uns component
good_uns = {
    "organism_ontology_term_id": "NCBITaxon:9606",
    "title": "A title",
    "default_embedding": "X_umap",
    "X_approximate_distribution": "normal",
    "batch_condition": ["is_primary_data"],
}

good_uns_mouse = {
    "organism_ontology_term_id": "NCBITaxon:10090",
    "title": "A title",
    "default_embedding": "X_umap",
    "X_approximate_distribution": "normal",
    "batch_condition": ["is_primary_data"],
}

good_uns_with_labels = {
    "organism_ontology_term_id": "NCBITaxon:9606",
    "organism": "Homo sapiens",
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
    "organism_ontology_term_id": "NCBITaxon:9606",
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
    "organism_ontology_term_id": "NCBITaxon:9606",
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
    "organism_ontology_term_id": "NCBITaxon:9606",
    "title": "A title",
    "default_embedding": "X_umap",
    "X_approximate_distribution": "normal",
    "batch_condition": ["is_primary_data"],
    "spatial": {"is_single": False},
}

good_uns_with_slide_seqV2_spatial = {
    "organism_ontology_term_id": "NCBITaxon:9606",
    "title": "A title",
    "default_embedding": "X_umap",
    "X_approximate_distribution": "normal",
    "batch_condition": ["is_primary_data"],
    "spatial": {"is_single": True},
}

# ---
# 4. Creating expression matrix,
# X has integer values and non_raw_X has real values
X = from_array(sparse.csr_matrix((good_obs.shape[0], good_var.shape[0]), dtype=numpy.float32))
for i in range(good_obs.shape[0]):
    for j in range(good_var.shape[0]):
        X[i, j] = i + j
non_raw_X = X.copy()
non_raw_X[0, 0] = 1.5

# ---
# 5.Creating valid obsm
good_obsm = {"X_umap": numpy.zeros([X.shape[0], 2])}
good_obsm_spatial = {"X_umap": numpy.zeros([X.shape[0], 2]), "spatial": numpy.zeros([X.shape[0], 2])}

# ---
# 6. Putting all the components created in the previous steps into minimal anndata that used for testing in
#   the unittests

# Valid anndata
adata = anndata.AnnData(X=X.copy(), obs=good_obs, uns=good_uns, obsm=good_obsm, var=good_var)
adata.raw = adata.copy()
adata.X = non_raw_X
adata.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# Anndata with "X" and "raw.X" but neither has actual raw values
adata_no_raw_values = anndata.AnnData(
    X=non_raw_X.copy(),
    obs=good_obs,
    uns=good_uns,
    obsm=good_obsm,
    var=good_var,
)
adata_no_raw_values.raw = adata_no_raw_values.copy()
adata_no_raw_values.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# Anndata with no obs nor var
adata_minimal = anndata.AnnData(X=X.copy(), uns=good_uns, obsm=good_obsm)

# Anndata with a expression matrix that is not raw g
adata_non_raw = anndata.AnnData(
    X=non_raw_X.copy(),
    obs=good_obs,
    uns=good_uns,
    obsm=good_obsm,
    var=good_var,
)

# Expected anndata with labels that the validator must write in obs and var
adata_with_labels = anndata.AnnData(
    X=X.copy(),
    obs=pd.concat([good_obs, obs_expected], axis=1),
    var=var_expected,
    uns=good_uns_with_labels,
    obsm=good_obsm,
)

# Expected anndata with colors for categorical obs fields
adata_with_colors = anndata.AnnData(X=X.copy(), obs=good_obs, uns=good_uns_with_colors, obsm=good_obsm, var=good_var)

# Expected anndata with Visium spatial data
adata_visium = anndata.AnnData(
    X=X.copy(), obs=good_obs_visium, uns=good_uns_with_visium_spatial, obsm=good_obsm_spatial, var=good_var
)
adata_visium.raw = adata_visium.copy()
adata_visium.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# Expected anndata with Slide-seqV2 spatial data
adata_slide_seqv2 = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_slide_seqv2,
    uns=good_uns_with_slide_seqV2_spatial,
    obsm=good_obsm_spatial,
    var=good_var,
)

adata_spatial_is_single_false = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_visium_is_single_false,
    uns=good_uns_with_is_single_false,
    obsm=good_obsm_spatial,
    var=good_var,
)

adata_mouse = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_mouse,
    uns=good_uns_mouse,
    obsm=good_obsm,
    var=good_var_mouse,
)

# Expected anndata with cell line data
adata_with_cell_line = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_cell_line,
    uns=good_uns,  # Use human organism since we test with human development stage terms
    obsm=good_obsm,
    var=good_var,  # Use human feature IDs since we're using human organism
)

# anndata for testing migration
unmigrated_obs = pd.DataFrame(
    [
        ["cell_type:1", "assay:1", "disease:1", "sex:1", "tissue:1", "sre:1", "development_stage:1", "tissue"],
        ["cell_type:1", "assay:1", "disease:1", "sex:1", "tissue:1", "sre:1", "development_stage:1", "tissue"],
    ],
    index=["X", "Y"],
    columns=[
        "cell_type_ontology_term_id",
        "assay_ontology_term_id",
        "disease_ontology_term_id",
        "sex_ontology_term_id",
        "tissue_ontology_term_id",
        "self_reported_ethnicity_ontology_term_id",
        "development_stage_ontology_term_id",
        "tissue_type",
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

unmigrated_X = sparse.csr_matrix(numpy.zeros([unmigrated_obs.shape[0], var_unmigrated.shape[0]], dtype=numpy.float32))

adata_with_labels_unmigrated = anndata.AnnData(
    X=unmigrated_X.copy(),
    obs=unmigrated_obs,
    uns=good_uns_with_labels,
    var=var_unmigrated,
    obsm={"X_umap": numpy.zeros([unmigrated_X.shape[0], 2])},
)
adata_with_labels_unmigrated.raw = adata_with_labels_unmigrated.copy()

# -----------------------------------------------------------------#
# Genetic perturbation fixtures (schema 7.1.0)

gp_id_1 = "CERS6-2"
gp_id_2 = "KARS-1"

# -------------------------
# Valid: non-control perturbations
# - No obs row may contain "na"
# - One row uses a single ID, the other uses two IDs joined with " || " in lexical order
good_obs_gene_perturbations = good_obs.copy()
good_obs_gene_perturbations["genetic_perturbation_id"] = pd.Categorical(
    [gp_id_1, f"{gp_id_1} || {gp_id_2}"],
    categories=[gp_id_1, gp_id_2, f"{gp_id_1} || {gp_id_2}"],
)
good_obs_gene_perturbations["genetic_perturbation_strategy"] = pd.Categorical(
    ["CRISPR knockout screen", "CRISPR knockout screen"],
    categories=[
        "control",
        "CRISPR activation screen",
        "CRISPR interference screen",
        "CRISPR knockout mutant",
        "CRISPR knockout screen",
        "no perturbations",
    ],
)

good_uns_with_gene_perturbations_curator = {
    **good_uns,
    "genetic_perturbations": {
        gp_id_1: {
            "role": "targeting",
            "protospacer_sequence": "CAGAGGGAGGAGAGAACCG",
            "protospacer_adjacent_motif": "3' NGG",
        },
        gp_id_2: {
            "role": "targeting",
            "protospacer_sequence": "GGGCCCTCCGGGAAGATGG",
            "protospacer_adjacent_motif": "3' NGG",
        },
    },
}

adata_gene_perturbations = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations,
    uns=good_uns_with_gene_perturbations_curator,
    obsm=good_obsm,
    var=good_var,
)
adata_gene_perturbations.raw = adata_gene_perturbations.copy()
adata_gene_perturbations.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Valid: control perturbations (strategy == "control" => role MUST be "control")
gp_control_id = "CTRL-1"
good_obs_gene_perturbations_control = good_obs.copy()
good_obs_gene_perturbations_control["genetic_perturbation_id"] = pd.Categorical(
    [gp_control_id, gp_control_id],
    categories=[gp_control_id],
)
good_obs_gene_perturbations_control["genetic_perturbation_strategy"] = pd.Categorical(
    ["control", "control"],
    categories=[
        "control",
        "CRISPR activation screen",
        "CRISPR interference screen",
        "CRISPR knockout mutant",
        "CRISPR knockout screen",
        "no perturbations",
    ],
)

good_uns_with_gene_perturbations_control = {
    **good_uns,
    "genetic_perturbations": {
        gp_control_id: {
            "role": "control",
            "protospacer_sequence": "TTTTTTTTTTTTTTTTTTTT",  # 20bp
            "protospacer_adjacent_motif": "3' NGG",
        }
    },
}

adata_gene_perturbations_control = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations_control,
    uns=good_uns_with_gene_perturbations_control,
    obsm=good_obsm,
    var=good_var,
)
adata_gene_perturbations_control.raw = adata_gene_perturbations_control.copy()
adata_gene_perturbations_control.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Invalid: strategy is "no perturbations" when id is not "na"
bad_obs_gene_perturbations_bad_strategy = good_obs.copy()
bad_obs_gene_perturbations_bad_strategy["genetic_perturbation_id"] = pd.Categorical(
    [gp_id_1, gp_id_1],
    categories=[gp_id_1],
)
bad_obs_gene_perturbations_bad_strategy["genetic_perturbation_strategy"] = pd.Categorical(
    ["no perturbations", "no perturbations"],
    categories=[
        "control",
        "CRISPR activation screen",
        "CRISPR interference screen",
        "CRISPR knockout mutant",
        "CRISPR knockout screen",
        "no perturbations",
    ],
)
adata_gene_perturbations_invalid_bad_strategy = anndata.AnnData(
    X=X.copy(),
    obs=bad_obs_gene_perturbations_bad_strategy,
    uns=good_uns_with_gene_perturbations_curator,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: multi-id not in ascending lexical order OR contains duplicates
bad_obs_gene_perturbations_bad_multi = good_obs.copy()
bad_obs_gene_perturbations_bad_multi["genetic_perturbation_id"] = pd.Categorical(
    [f"{gp_id_2} || {gp_id_1}", f"{gp_id_1} || {gp_id_1}"],  # wrong order, duplicate
    categories=[gp_id_1, gp_id_2, f"{gp_id_2} || {gp_id_1}", f"{gp_id_1} || {gp_id_1}"],
)
bad_obs_gene_perturbations_bad_multi["genetic_perturbation_strategy"] = pd.Categorical(
    ["CRISPR knockout screen", "CRISPR knockout screen"],
    categories=[
        "control",
        "CRISPR activation screen",
        "CRISPR interference screen",
        "CRISPR knockout mutant",
        "CRISPR knockout screen",
        "no perturbations",
    ],
)
adata_gene_perturbations_invalid_bad_multi = anndata.AnnData(
    X=X.copy(),
    obs=bad_obs_gene_perturbations_bad_multi,
    uns=good_uns_with_gene_perturbations_curator,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: obs references an ID missing from uns['genetic_perturbations']
bad_obs_gene_perturbations_missing_key = good_obs.copy()
bad_obs_gene_perturbations_missing_key["genetic_perturbation_id"] = pd.Categorical(
    ["MISSING_ID", "MISSING_ID"],
    categories=["MISSING_ID"],
)
bad_obs_gene_perturbations_missing_key["genetic_perturbation_strategy"] = pd.Categorical(
    ["CRISPR knockout screen", "CRISPR knockout screen"],
    categories=[
        "control",
        "CRISPR activation screen",
        "CRISPR interference screen",
        "CRISPR knockout mutant",
        "CRISPR knockout screen",
        "no perturbations",
    ],
)
adata_gene_perturbations_invalid_missing_key = anndata.AnnData(
    X=X.copy(),
    obs=bad_obs_gene_perturbations_missing_key,
    uns=good_uns_with_gene_perturbations_curator,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: strategy == "control" but role is not "control"
bad_uns_gene_perturbations_control_role_mismatch = {
    **good_uns,
    "genetic_perturbations": {
        gp_control_id: {
            "role": "targeting",  # invalid when strategy is "control"
            "protospacer_sequence": "TTTTTTTTTTTTTTTTTTTT",
            "protospacer_adjacent_motif": "3' NGG",
        }
    },
}
adata_gene_perturbations_invalid_control_role_mismatch = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations_control,
    uns=bad_uns_gene_perturbations_control_role_mismatch,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: Discover-only key present in curator submission
bad_uns_gene_perturbations_contains_derived = {
    **good_uns,
    "genetic_perturbations": {
        gp_id_1: {
            "role": "targeting",
            "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
            "protospacer_adjacent_motif": "3' NGG",
            "derived_genomic_regions": ["16:75647615-75647633(-)"],  # Discover-only
        }
    },
}
adata_gene_perturbations_invalid_contains_derived = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_gene_perturbations_contains_derived,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: uns has genetic_perturbations but obs is missing required columns
# This tests the cross-validation check that requires obs columns when uns has genetic_perturbations
adata_gene_perturbations_invalid_missing_obs_columns = anndata.AnnData(
    X=X.copy(),
    obs=good_obs.copy(),  # good_obs doesn't have genetic_perturbation_id or genetic_perturbation_strategy
    uns=good_uns_with_gene_perturbations_curator,  # uns has genetic_perturbations
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: key name contains whitespace (violates key_pattern)
bad_uns_key_with_whitespace = {
    **good_uns,
    "genetic_perturbations": {
        "guide 1": {  # Contains space
            "role": "targeting",
            "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
            "protospacer_adjacent_motif": "3' NGG",
        }
    },
}
adata_gene_perturbations_invalid_key_whitespace = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_key_with_whitespace,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: key name contains slash (violates key_pattern)
bad_uns_key_with_slash = {
    **good_uns,
    "genetic_perturbations": {
        "guide/1": {  # Contains slash
            "role": "targeting",
            "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
            "protospacer_adjacent_motif": "3' NGG",
        }
    },
}
adata_gene_perturbations_invalid_key_slash = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_key_with_slash,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: key name contains comma (violates key_pattern)
bad_uns_key_with_comma = {
    **good_uns,
    "genetic_perturbations": {
        "guide,1": {  # Contains comma
            "role": "targeting",
            "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
            "protospacer_adjacent_motif": "3' NGG",
        }
    },
}
adata_gene_perturbations_invalid_key_comma = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_key_with_comma,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: key name contains single quote (violates key_pattern)
bad_uns_key_with_quote = {
    **good_uns,
    "genetic_perturbations": {
        "guide'1": {  # Contains single quote
            "role": "targeting",
            "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
            "protospacer_adjacent_motif": "3' NGG",
        }
    },
}
adata_gene_perturbations_invalid_key_quote = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_key_with_quote,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: protospacer_sequence too short (< 14bp)
bad_uns_protospacer_too_short = {
    **good_uns,
    "genetic_perturbations": {
        gp_id_1: {
            "role": "targeting",
            "protospacer_sequence": "ACGTACGTACGT",  # Only 12bp, need 14-22
            "protospacer_adjacent_motif": "3' NGG",
        }
    },
}
adata_gene_perturbations_invalid_protospacer_short = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_protospacer_too_short,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: protospacer_sequence too long (> 22bp)
bad_uns_protospacer_too_long = {
    **good_uns,
    "genetic_perturbations": {
        gp_id_1: {
            "role": "targeting",
            "protospacer_sequence": "ACGTACGTACGTACGTACGTACG",  # 23bp, need 14-22
            "protospacer_adjacent_motif": "3' NGG",
        }
    },
}
adata_gene_perturbations_invalid_protospacer_long = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_protospacer_too_long,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: protospacer_sequence contains non-ACGT characters
bad_uns_protospacer_invalid_chars = {
    **good_uns,
    "genetic_perturbations": {
        gp_id_1: {
            "role": "targeting",
            "protospacer_sequence": "ACGTACGTACGTACGN",  # Contains 'N', need only ACGT
            "protospacer_adjacent_motif": "3' NGG",
        }
    },
}
adata_gene_perturbations_invalid_protospacer_chars = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_protospacer_invalid_chars,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: protospacer_adjacent_motif wrong format (doesn't start with "3' ")
bad_uns_pam_wrong_format = {
    **good_uns,
    "genetic_perturbations": {
        gp_id_1: {
            "role": "targeting",
            "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
            "protospacer_adjacent_motif": "NGG",  # Missing "3' " prefix
        }
    },
}
adata_gene_perturbations_invalid_pam_format = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_pam_wrong_format,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: missing required field 'role'
bad_uns_missing_role = {
    **good_uns,
    "genetic_perturbations": {
        gp_id_1: {
            # "role": "targeting",  # Missing required field
            "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
            "protospacer_adjacent_motif": "3' NGG",
        }
    },
}
adata_gene_perturbations_invalid_missing_role = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_missing_role,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: missing required field 'protospacer_sequence'
bad_uns_missing_protospacer = {
    **good_uns,
    "genetic_perturbations": {
        gp_id_1: {
            "role": "targeting",
            # "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",  # Missing required field
            "protospacer_adjacent_motif": "3' NGG",
        }
    },
}
adata_gene_perturbations_invalid_missing_protospacer = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_missing_protospacer,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: missing required field 'protospacer_adjacent_motif'
bad_uns_missing_pam = {
    **good_uns,
    "genetic_perturbations": {
        gp_id_1: {
            "role": "targeting",
            "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
            # "protospacer_adjacent_motif": "3' NGG",  # Missing required field
        }
    },
}
adata_gene_perturbations_invalid_missing_pam = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_missing_pam,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: role has invalid enum value
bad_uns_invalid_role = {
    **good_uns,
    "genetic_perturbations": {
        gp_id_1: {
            "role": "unknown",  # Invalid, must be 'control' or 'targeting'
            "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
            "protospacer_adjacent_motif": "3' NGG",
        }
    },
}
adata_gene_perturbations_invalid_role_enum = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_invalid_role,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: reserved key 'derived_features' present
bad_uns_contains_derived_features = {
    **good_uns,
    "genetic_perturbations": {
        gp_id_1: {
            "role": "targeting",
            "protospacer_sequence": "GCTGCTGCTGCTGCTGCTGA",
            "protospacer_adjacent_motif": "3' NGG",
            "derived_features": ["ENSG00000012345"],  # Reserved key
        }
    },
}
adata_gene_perturbations_invalid_contains_derived_features = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),
    uns=bad_uns_contains_derived_features,
    obsm=good_obsm,
    var=good_var,
)

# -------------------------
# Invalid: obs has genetic_perturbation columns but uns is missing genetic_perturbations
# This tests the dependency that requires uns['genetic_perturbations'] when obs columns are present
adata_gene_perturbations_invalid_obs_without_uns = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_gene_perturbations.copy(),  # Has genetic_perturbation_id and genetic_perturbation_strategy
    uns=good_uns,  # Does NOT have genetic_perturbations
    obsm=good_obsm,
    var=good_var,
)

# -----------------------------------------------------------------#
# Experimental condition fixtures (schema 7.1.0)

# -------------------------
# Valid: single CHEBI term (aspirin = chemical perturbation)
good_obs_ec_chebi = good_obs.copy()
good_obs_ec_chebi["experimental_condition_ontology_term_id"] = pd.Categorical(
    ["CHEBI:15365", "CHEBI:15365"],
    categories=["CHEBI:15365"],
)

adata_ec_valid_chebi = anndata.AnnData(X=X.copy(), obs=good_obs_ec_chebi, uns=good_uns, obsm=good_obsm, var=good_var)
adata_ec_valid_chebi.raw = adata_ec_valid_chebi.copy()
adata_ec_valid_chebi.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Valid: anti-uniprot protein term (antibody)
good_obs_ec_anti_uniprot = good_obs.copy()
good_obs_ec_anti_uniprot["experimental_condition_ontology_term_id"] = pd.Categorical(
    ["anti-uniprot:Q99467", "anti-uniprot:Q99467"],
    categories=["anti-uniprot:Q99467"],
)

adata_ec_valid_anti_uniprot = anndata.AnnData(
    X=X.copy(), obs=good_obs_ec_anti_uniprot, uns=good_uns, obsm=good_obsm, var=good_var
)
adata_ec_valid_anti_uniprot.raw = adata_ec_valid_anti_uniprot.copy()
adata_ec_valid_anti_uniprot.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Valid: uniprot protein term (no anti- prefix)
good_obs_ec_uniprot = good_obs.copy()
good_obs_ec_uniprot["experimental_condition_ontology_term_id"] = pd.Categorical(
    ["uniprot:Q99467", "uniprot:Q99467"],
    categories=["uniprot:Q99467"],
)

adata_ec_valid_uniprot = anndata.AnnData(
    X=X.copy(), obs=good_obs_ec_uniprot, uns=good_uns, obsm=good_obsm, var=good_var
)
adata_ec_valid_uniprot.raw = adata_ec_valid_uniprot.copy()
adata_ec_valid_uniprot.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Valid: EFO:0001702 temperature term
good_obs_ec_temperature = good_obs.copy()
good_obs_ec_temperature["experimental_condition_ontology_term_id"] = pd.Categorical(
    ["EFO:0001702", "EFO:0001702"],
    categories=["EFO:0001702"],
)

adata_ec_valid_temperature = anndata.AnnData(
    X=X.copy(), obs=good_obs_ec_temperature, uns=good_uns, obsm=good_obsm, var=good_var
)
adata_ec_valid_temperature.raw = adata_ec_valid_temperature.copy()
adata_ec_valid_temperature.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Valid: EFO:0002755 diet root term
good_obs_ec_diet_root = good_obs.copy()
good_obs_ec_diet_root["experimental_condition_ontology_term_id"] = pd.Categorical(
    ["EFO:0002755", "EFO:0002755"],
    categories=["EFO:0002755"],
)

adata_ec_valid_diet_root = anndata.AnnData(
    X=X.copy(), obs=good_obs_ec_diet_root, uns=good_uns, obsm=good_obsm, var=good_var
)
adata_ec_valid_diet_root.raw = adata_ec_valid_diet_root.copy()
adata_ec_valid_diet_root.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Valid: EFO:0002757 (a descendant of diet EFO:0002755)
good_obs_ec_diet_desc = good_obs.copy()
good_obs_ec_diet_desc["experimental_condition_ontology_term_id"] = pd.Categorical(
    ["EFO:0002757", "EFO:0002757"],
    categories=["EFO:0002757"],
)

adata_ec_valid_diet_descendant = anndata.AnnData(
    X=X.copy(), obs=good_obs_ec_diet_desc, uns=good_uns, obsm=good_obsm, var=good_var
)
adata_ec_valid_diet_descendant.raw = adata_ec_valid_diet_descendant.copy()
adata_ec_valid_diet_descendant.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Valid: multi-term sorted " || " (CHEBI + EFO:0002755)
good_obs_ec_multi = good_obs.copy()
good_obs_ec_multi["experimental_condition_ontology_term_id"] = pd.Categorical(
    ["CHEBI:15365 || EFO:0002755", "CHEBI:15365 || EFO:0002755"],
    categories=["CHEBI:15365 || EFO:0002755"],
)

adata_ec_valid_multi = anndata.AnnData(X=X.copy(), obs=good_obs_ec_multi, uns=good_uns, obsm=good_obsm, var=good_var)
adata_ec_valid_multi.raw = adata_ec_valid_multi.copy()
adata_ec_valid_multi.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Valid: absent (column not present) — field is not required
adata_ec_absent = adata  # reuse base adata without the column

# -------------------------
# Valid: all "na" values — column absent is preferred but this must fail
# (forbidden_when_all_na)
bad_obs_ec_all_na = good_obs.copy()
bad_obs_ec_all_na["experimental_condition_ontology_term_id"] = pd.Categorical(
    ["na", "na"],
    categories=["na"],
)

adata_ec_invalid_all_na = anndata.AnnData(X=X.copy(), obs=bad_obs_ec_all_na, uns=good_uns, obsm=good_obsm, var=good_var)
adata_ec_invalid_all_na.raw = adata_ec_invalid_all_na.copy()
adata_ec_invalid_all_na.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Invalid: forbidden CHEBI term (CHEBI:23367 = molecular entity)
bad_obs_ec_forbidden_chebi = good_obs.copy()
bad_obs_ec_forbidden_chebi["experimental_condition_ontology_term_id"] = pd.Categorical(
    ["CHEBI:23367", "CHEBI:23367"],
    categories=["CHEBI:23367"],
)

adata_ec_invalid_forbidden_chebi = anndata.AnnData(
    X=X.copy(), obs=bad_obs_ec_forbidden_chebi, uns=good_uns, obsm=good_obsm, var=good_var
)
adata_ec_invalid_forbidden_chebi.raw = adata_ec_invalid_forbidden_chebi.copy()
adata_ec_invalid_forbidden_chebi.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Invalid: EFO term that is not allowed (not EFO:0001702 / not a diet descendant)
bad_obs_ec_efo_not_allowed = good_obs.copy()
bad_obs_ec_efo_not_allowed["experimental_condition_ontology_term_id"] = pd.Categorical(
    ["EFO:0009899", "EFO:0009899"],  # 10x 3' v2 assay — valid EFO but not allowed here
    categories=["EFO:0009899"],
)

adata_ec_invalid_efo_not_allowed = anndata.AnnData(
    X=X.copy(), obs=bad_obs_ec_efo_not_allowed, uns=good_uns, obsm=good_obsm, var=good_var
)
adata_ec_invalid_efo_not_allowed.raw = adata_ec_invalid_efo_not_allowed.copy()
adata_ec_invalid_efo_not_allowed.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Invalid: multi-term unsorted
bad_obs_ec_unsorted = good_obs.copy()
bad_obs_ec_unsorted["experimental_condition_ontology_term_id"] = pd.Categorical(
    ["EFO:0002755 || CHEBI:15365", "EFO:0002755 || CHEBI:15365"],  # unsorted (should be CHEBI first)
    categories=["EFO:0002755 || CHEBI:15365"],
)

adata_ec_invalid_multi_unsorted = anndata.AnnData(
    X=X.copy(), obs=bad_obs_ec_unsorted, uns=good_uns, obsm=good_obsm, var=good_var
)
adata_ec_invalid_multi_unsorted.raw = adata_ec_invalid_multi_unsorted.copy()
adata_ec_invalid_multi_unsorted.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Invalid: reserved column 'experimental_condition' manually set by curator
bad_obs_ec_reserved = good_obs.copy()
bad_obs_ec_reserved["experimental_condition"] = pd.Categorical(
    ["some label", "some label"],
    categories=["some label"],
)

adata_ec_invalid_reserved_label_column = anndata.AnnData(
    X=X.copy(), obs=bad_obs_ec_reserved, uns=good_uns, obsm=good_obsm, var=good_var
)
adata_ec_invalid_reserved_label_column.raw = adata_ec_invalid_reserved_label_column.copy()
adata_ec_invalid_reserved_label_column.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Invalid: reserved column 'perturbation_types' manually set by curator
bad_obs_perturbation_types_reserved = good_obs.copy()
bad_obs_perturbation_types_reserved["perturbation_types"] = pd.Categorical(
    ["chemical", "chemical"],
    categories=["chemical"],
)

adata_ec_invalid_reserved_perturbation_types = anndata.AnnData(
    X=X.copy(), obs=bad_obs_perturbation_types_reserved, uns=good_uns, obsm=good_obsm, var=good_var
)
adata_ec_invalid_reserved_perturbation_types.raw = adata_ec_invalid_reserved_perturbation_types.copy()
adata_ec_invalid_reserved_perturbation_types.raw.var.drop("feature_is_filtered", axis=1, inplace=True)

# -------------------------
# Valid: experimental_condition_ontology_term_id + genetic_perturbation_id together
# (used for perturbation_types derivation: chemical + genetic)
good_obs_ec_and_gp = good_obs_gene_perturbations.copy()
good_obs_ec_and_gp["experimental_condition_ontology_term_id"] = pd.Categorical(
    ["CHEBI:15365", "CHEBI:15365"],
    categories=["CHEBI:15365"],
)

adata_ec_valid_with_gp = anndata.AnnData(
    X=X.copy(),
    obs=good_obs_ec_and_gp,
    uns=good_uns_with_gene_perturbations_curator,
    obsm=good_obsm,
    var=good_var,
)
adata_ec_valid_with_gp.raw = adata_ec_valid_with_gp.copy()
adata_ec_valid_with_gp.raw.var.drop("feature_is_filtered", axis=1, inplace=True)
