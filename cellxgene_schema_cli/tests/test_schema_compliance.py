import unittest
import pandas as pd
import numpy
from cellxgene_schema import validate
import fixtures.examples_validate as examples

# Tests for schema compliance of an AnnData object


class TestExpressionMatrix(unittest.TestCase):

    """
    Fail cases for expression matrices (anndata.X and anndata.raw.X)
    """

    def setUp(self):

        self.validator = validate.Validator()
        self.validator.adata = examples.adata.copy()

    def test_shapes(self):

        """
        All matrix layers MUST have the same shape, and have the same cell labels and gene labels.
        """

        # Creates a raw layer
        self.validator.adata.raw = self.validator.adata
        self.validator.adata.X = examples.adata_non_raw.X.copy()
        self.validator.adata.uns["X_normalization"] = "CPM"

        # remove one gene
        self.validator.adata = self.validator.adata[:, 1:]
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: Number of genes in X (3) is different than raw.X (4)"],
        )

    def test_sparsity(self):

        """
        In any layer, if a matrix has 50% or more values that are zeros, it is STRONGLY RECOMMENDED that
        the matrix be encoded as a scipy.sparse.csr_matrix
        """

        self.validator.adata.X = self.validator.adata.X.toarray()
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.warnings,
            [
                "WARNING: Sparsity of 'X' is 1.0 which is greater than 0.5, "
                "and it is not a 'scipy.sparse.csr_matrix'. It is "
                "STRONGLY RECOMMENDED to use this type of matrix for "
                "the given sparsity."
            ],
        )

    def test_raw_existance(self):

        """
        Except for ATAC-seq and methylation data, raw data is REQUIRED
        """

        self.validator.adata.uns["X_normalization"] = "CPM"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: uns['X_normalization'] is 'CPM' but raw data seems "
                "to be in X, if X is raw then uns['X_normalization'] "
                "MUST be 'none'."
            ],
        )


class TestValidAnndata(unittest.TestCase):

    """
    Tests a valid anndata. Most other tests below modify this anndata object and test for failure cases
    """

    def setUp(self):
        self.validator = validate.Validator()
        self.validator.adata = examples.adata.copy()

    def test_valid_anndata(self):
        self.validator.validate_adata()
        self.assertFalse(self.validator.errors)


class TestObs(unittest.TestCase):

    """
    Fail cases in adata.uns
    """

    def setUp(self):
        self.validator = validate.Validator()
        self.validator.adata = examples.adata.copy()

    def test_column_presence(self):
        """
        obs is a pandas.DataFrame. Curators MUST annotate the following columns in the obs dataframe.
        """

        columns = [
            "assay_ontology_term_id",
            "development_stage_ontology_term_id",
            "disease_ontology_term_id",
            "ethnicity_ontology_term_id",
            "is_primary_data",
            "sex_ontology_term_id",
            "tissue_ontology_term_id",
        ]

        for column in columns:
            self.validator.errors = []
            self.validator.adata = examples.adata.copy()

            self.validator.adata.obs.drop(column, axis=1, inplace=True)
            # Remove batch condition because it has a dependency with is_primary_data
            self.validator.adata.uns.pop("batch_condition")

            self.validator.validate_adata()
            self.assertEqual(
                self.validator.errors,
                [f"ERROR: Dataframe 'obs' is missing " f"column '{column}'."],
            )

    def test_column_presence_organism(self):
        """
        obs is a pandas.DataFrame. Curators MUST annotate the following columns in the obs dataframe.

        A separate check is need for organism_ontology_term_id because removing from anndata results in multiple
        errors given that other columns depend on its presence
        """

        self.validator.adata.obs.drop("organism_ontology_term_id", axis=1, inplace=True)
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: Dataframe 'obs' is missing column "
                "'organism_ontology_term_id'.",
                "ERROR: Checking values with dependencies failed for "
                "adata.obs['ethnicity_ontology_term_id'], this is likely due "
                "to missing dependent column in adata.obs",
                "ERROR: Checking values with dependencies failed for "
                "adata.obs['development_stage_ontology_term_id'], this is likely due "
                "to missing dependent column in adata.obs",
            ],
        )

    def test_obsolete_term_id(self):
        """
        Terms documented as obsolete in an ontology MUST NOT be used. For example, EFO:0009310
        for obsolete_10x v2 was marked as obsolete in EFO version 3.31.0 and replaced by
        EFO:0009899 for 10x 3' v2.
        """

        # Not a valid term
        self.validator.adata.obs["assay_ontology_term_id"][0] = "EFO:0009310"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0009310' in 'assay_ontology_term_id' is a deprecated term id of 'EFO'",
                "ERROR: 'EFO:0009310' in 'assay_ontology_term_id' is not a children term id "
                "of '[['EFO:0002772', 'EFO:0010183']]'",
            ],
        )

    def test_assay_ontology_term_id(self):

        """
        assay_ontology_term_id categorical with str categories. categorical with str categories.
        This MUST be an EFO term and either child of "EFO:0002772" or "EFO:0010183"
        """

        # Not a valid term
        self.validator.adata.obs["assay_ontology_term_id"][0] = "CL:000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'CL:000001' in 'assay_ontology_term_id' is not a valid "
                "ontology term id of 'EFO'",
                "ERROR: 'CL:000001' in 'assay_ontology_term_id' is not a children "
                "term id of '[['EFO:0002772', 'EFO:0010183']]'",
            ],
        )

        # Not a valid child
        self.validator.adata.obs["assay_ontology_term_id"][0] = "EFO:0000001"
        self.validator.errors = []
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'assay_ontology_term_id' is not a "
                "children term id of '[['EFO:0002772', 'EFO:0010183']]'"
            ],
        )

    def test_cell_type_ontology_term_id(self):

        """
        cell_type_ontology_term_id categorical with str categories. This MUST be a CL term.
        """

        # Not a valid term
        self.validator.adata.obs["cell_type_ontology_term_id"][0] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'cell_type_ontology_term_id' is not a valid "
                "ontology term id of 'CL'"
            ],
        )

    def test_development_stage_ontology_term_id(self):

        """
        development_stage_ontology_term_id categorical with str categories. If unavailable, this MUST be "unknown".
        If organism_ontolology_term_id is "NCBITaxon:9606" for Homo sapiens,
        this MUST be the most accurate HsapDv term.
        If organism_ontolology_term_id is "NCBITaxon:10090" for Mus musculus,
        this MUST be the most accurate MmusDv term
        All other it MUST be children of UBERON:0000105 and not UBERON:0000071
        """

        # If organism_ontolology_term_id is "NCBITaxon:9606" for Homo sapiens,
        # this MUST be the most accurate HsapDv term.
        self.validator.adata.obs["organism_ontology_term_id"][0] = "NCBITaxon:9606"
        self.validator.adata.obs["development_stage_ontology_term_id"][
            0
        ] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'development_stage_ontology_term_id' is "
                "not a valid ontology term id of 'HsapDv'"
            ],
        )

        # If organism_ontolology_term_id is "NCBITaxon:10090" for Mus musculus,
        # this MUST be the most accurate MmusDv term
        self.validator.errors = []
        self.validator.adata.obs["organism_ontology_term_id"][0] = "NCBITaxon:10090"
        self.validator.adata.obs["development_stage_ontology_term_id"][
            0
        ] = "EFO:0000001"
        self.validator.adata.obs["ethnicity_ontology_term_id"][0] = "na"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'development_stage_ontology_term_id' is "
                "not a valid ontology term id of 'MmusDv'"
            ],
        )

        # All other it MUST be children of UBERON:0000105 and not UBERON:0000071
        self.validator.errors = []
        self.validator.adata.obs["organism_ontology_term_id"][0] = "NCBITaxon:10114"
        self.validator.adata.obs["development_stage_ontology_term_id"][
            0
        ] = "EFO:0000001"
        self.validator.adata.obs["ethnicity_ontology_term_id"][0] = "na"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'development_stage_ontology_term_id' is "
                "not a valid ontology term id of 'UBERON'",
                "ERROR: 'EFO:0000001' in 'development_stage_ontology_term_id' is not "
                "a children term id of '[['UBERON:0000105']]'",
            ],
        )

    def test_disease_ontology_term_id(self):

        """
        disease_ontology_term_id categorical with str categories. This MUST be a MONDO term or
        PATO:0000461 for normal or healthy.
        """

        # If organism_ontolology_term_id is "NCBITaxon:9606" for Homo sapiens,
        # this MUST be the most accurate HsapDv term.
        self.validator.adata.obs["disease_ontology_term_id"][0] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'disease_ontology_term_id' is not a "
                "valid ontology term id of 'MONDO, PATO'"
            ],
        )

    def test_ethnicity_ontology_term_id(self):

        """
        ethnicity_ontology_term_id categorical with str categories.
        If organism_ontolology_term_id is "NCBITaxon:9606" for Homo sapiens,
        this MUST be either a HANCESTRO term or "unknown" if unavailable.
        Otherwise, for all other organisms this MUST be "na".
        """

        # If organism_ontolology_term_id is "NCBITaxon:9606" for Homo sapiens,
        # this MUST be either a HANCESTRO term or "unknown" if unavailable.
        self.validator.adata.obs["organism_ontology_term_id"][0] = "NCBITaxon:9606"
        self.validator.adata.obs["ethnicity_ontology_term_id"][0] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'ethnicity_ontology_term_id' is "
                "not a valid ontology term id of 'HANCESTRO'"
            ],
        )

        # Otherwise, for all other organisms this MUST be "na".
        self.validator.errors = []
        self.validator.adata.obs["organism_ontology_term_id"][0] = "NCBITaxon:10090"
        self.validator.adata.obs["development_stage_ontology_term_id"][
            0
        ] = "MmusDv:0000003"
        self.validator.adata.obs["ethnicity_ontology_term_id"][0] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'ethnicity_ontology_term_id' is not a "
                "valid value of 'ethnicity_ontology_term_id'"
            ],
        )

    def test_organism_ontology_term_id(self):

        """
        organism_ontology_term_id categorical with str categories. This MUST be a child of NCBITaxon:33208.
        """

        # If organism_ontolology_term_id is "NCBITaxon:9606" for Homo sapiens,
        # this MUST be either a HANCESTRO term or "unknown" if unavailable.
        self.validator.adata.obs["organism_ontology_term_id"][0] = "EFO:0000001"
        self.validator.adata.obs["development_stage_ontology_term_id"][0] = "unknown"
        self.validator.adata.obs["ethnicity_ontology_term_id"][0] = "na"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'organism_ontology_term_id' is not a valid "
                "ontology term id of 'NCBITaxon'"
            ],
        )

    def test_tissue_ontology_term_id(self):

        """
        tissue_ontology_term_id categorical with str categories. This MUST be the term that best describes the tissue
        that this cell was derived from, depending on the type of biological sample:
            Cell Culture - MUST be a CL term appended with " (cell culture)"
            Organoid - MUST be an UBERON term appended with " (organoid)"
            Enriched, Sorted,or Isolated Cells from a Tissue - MUST be an UBERON or CL term and SHOULD NOT use terms
                that do not capture the tissue of origin (e.g. In the case of CD3+ kidney cells, use "UBERON:0002113"
                for kidney instead of "CL:000084" for T cell. However, in the case of EPCAM+ cervical cells,
                use "CL:000066" for epithelial cell of the cervix.)
        """

        self.validator.adata.obs["tissue_ontology_term_id"][0] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'tissue_ontology_term_id' is not a "
                "valid ontology term id of 'UBERON, CL'"
            ],
        )

        # Cell Culture - MUST be a CL term appended with " (cell culture)"
        self.errors = []
        self.validator.adata.obs["tissue_ontology_term_id"][
            0
        ] = "CL:0000057 (CELL culture)"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'CL:0000057 (CELL culture)' in 'tissue_ontology_term_id' is "
                "not a valid ontology term id of 'UBERON, CL'"
            ],
        )

        # Organoid - MUST be an UBERON term appended with " (organoid)"
        self.errors = []
        self.validator.adata.obs["tissue_ontology_term_id"][0] = "CL:0000057 (ORGANOID)"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'CL:0000057 (ORGANOID)' in 'tissue_ontology_term_id' is "
                "not a valid ontology term id of 'UBERON, CL'"
            ],
        )

    def test_sex_ontology_term_id(self):

        """
        sex_ontology_term_id categorical with str categories.
        This MUST be a child of PATO:0001894 for phenotypic sex or "unknown" if unavailable
        """

        self.validator.adata.obs["sex_ontology_term_id"][0] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'sex_ontology_term_id' is "
                "not a valid ontology term id of 'PATO'"
            ],
        )

    def test_is_primary_data(self):

        """
        is_primary_data	bool. This MUST be True if this is the canonical instance of this cellular
        observation and False if not. This is commonly False
        for meta-analyses reusing data or for secondary views of data.
        """

        self.validator.adata.obs["is_primary_data"] = "FALSE"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: Column 'is_primary_data' in dataframe 'obs' "
                "must be boolean not 'object'"
            ],
        )


class TestVar(unittest.TestCase):

    """
    Fail cases in adata.var
    """

    def setUp(self):
        self.validator = validate.Validator()
        self.validator.adata = examples.adata.copy()

    def test_check_unique_var(self):

        """
        var.index MUST contain unique ENSEMBL gene identifiers for features.
        """

        # Duplicate 1st row in var and assigned to 2nd
        new_index = list(self.validator.adata.var_names)
        new_index[1] = new_index[0]
        self.validator.adata.var_names = new_index
        self.validator.adata.var.iloc[1, :] = self.validator.adata.var.iloc[0, :]

        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: Column 'index' in dataframe 'var' is not unique."],
        )

    def test_feature_id(self):

        """
        feature_id (var.index) str.
        If the feature_biotype is "gene" then this MUST be an ENSEMBL term.
        If the feature_biotype is "spike-in" then this MUST be an ERCC Spike-In identifier.
        """

        # If the feature_biotype is "gene" then this MUST be an ENSEMBL term.
        # First not and ENSEMBL ID
        new_index = list(self.validator.adata.var_names)
        new_index[0] = "ENSEBML_NOGENE"
        self.validator.adata.var_names = new_index
        self.validator.adata.var["feature_biotype"][0] = "gene"

        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: Could not infer organism from feature ID 'ENSEBML_NOGENE' "
                "in 'var', make sure it is a valid ID"
            ],
        )

        # Second, something that looks like ENSEBML id but it isn't
        new_index = list(self.validator.adata.var_names)
        new_index[0] = "ENSG000"
        self.validator.adata.var_names = new_index
        self.validator.adata.var["feature_biotype"][0] = "gene"

        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: 'ENSG000' is not a valid feature ID in 'var'"],
        )

        # If the feature_biotype is "spike-in" then this MUST be an ERCC Spike-In identifier.
        new_index = list(self.validator.adata.var_names)
        new_index[0] = "ERCC-000000"
        self.validator.adata.var_names = new_index
        self.validator.adata.var["feature_biotype"][0] = "spike-in"

        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: 'ERCC-000000' is not a valid feature ID in 'var'"],
        )


class TestUns(unittest.TestCase):

    """
    Fail cases in adata.uns
    """

    def setUp(self):
        self.validator = validate.Validator()
        self.validator.adata = examples.adata.copy()

    def test_required_fields_schema_version(self):

        """
        Curators MUST annotate the following keys and values in uns: schema_version, title, X_normalization
        """

        del self.validator.adata.uns["schema_version"]
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: adata has no schema definition in 'adata.uns'. "
                "Validation cannot be performed."
            ],
        )

    def test_required_fields_title(self):

        """
        Curators MUST annotate the following keys and values in uns: schema_version, title, X_normalization
        """

        del self.validator.adata.uns["title"]
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors, ["ERROR: 'title' in 'uns' is not present"]
        )

    def test_required_fields_X_normalization(self):

        """
        Curators MUST annotate the following keys and values in uns: schema_version, title, X_normalization
        """

        del self.validator.adata.uns["X_normalization"]
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors, ["ERROR: 'X_normalization' in 'uns' is not present"]
        )

    def test_leading_trailing_double_spaces_in_strings(self):

        """
        The following sequences MUST NOT appear in str types documented in the schema:
            Leading control or space separators - ”     This is an example”
            Trailing control or space separators - “This is an example     ”
            Multiple (internal) control or space separators - "This     is an example"
        """

        self.validator.adata.uns["title"] = " There is a leading space"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: ' There is a leading space' in 'uns['title']' is not valid, it contains leading spaces."
            ],
        )

        self.validator.adata.uns["title"] = "There is a trailing space "
        self.validator.errors = []
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'There is a trailing space ' in 'uns['title']' is not valid, it contains trailing spaces."
            ],
        )

        self.validator.adata.uns["title"] = "There are   double   spaces"
        self.validator.errors = []
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'There are   double   spaces' in 'uns['title']' is not valid, it contains double spaces."
            ],
        )

    def test_schema_version(self):

        """
        Schema_version, This MUST be "2.0.0".
        """

        self.validator.adata.uns["schema_version"] = "1.0.0"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: Schema version '1.0.0' is not supported. Validation "
                "cannot be performed."
            ],
        )

    def test_title(self):

        """
        Title MUST be a string
        """

        # list instead of string
        self.validator.adata.uns["title"] = ["title"]
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: '['title']' in 'uns['title']' is not valid, "
                "it must be a string."
            ],
        )

    def test_X_normalization_is_str(self):

        """
        X_normalization str.
        """

        # list instead of string
        self.validator.adata.uns["X_normalization"] = ["normalization"]
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: '['normalization']' in 'uns['X_normalization']' is "
                "not valid, it must be a string."
            ],
        )

    def test_X_normalization_not_raw(self):

        """
        X_normalization str. This SHOULD describe the method used to normalize the data stored in AnnData X.
        If data in X are raw, this SHOULD be "none".

        FAIL CASE for when X_normalization was set to "none" but X is not raw data
        """

        # Assign a real value to X while X_normalization is 'none'
        self.validator.adata.X[0, 0] = 1.5
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: Matrix 'X' seems to be the raw matrix but not all of "
                "its values are integers."
            ],
        )

    def test_batch_condition_is_list(self):

        """
        batch_condition list[str]
        """

        self.validator.adata.uns["batch_condition"] = "cell_type_ontology_term_id"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'cell_type_ontology_term_id' in 'uns['batch_condition']' "
                "is not valid, it must be a list or numpy array"
            ],
        )

    def test_batch_condition_is_column_from_obs(self):

        """
        batch_condition list[str]. str values MUST refer to cell metadata keys in obs.
        """

        self.validator.adata.uns["batch_condition"] = ["NO_COLUMN"]
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: Value 'NO_COLUMN' of list 'batch_condition' is not a "
                "column in 'adata.obs'"
            ],
        )

    def test_default_embedding_is_str(self):

        """
        Default_embedding str.
        """

        self.validator.adata.uns["default_embedding"] = ["X_umap"]
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: '['X_umap']' in 'uns['default_embedding']' is not valid, "
                "it must be a string."
            ],
        )

    def test_default_embedding_is_key_from_obsm(self):

        """
        Default_embedding str. The value MUST match a key to an embedding in obsm
        """

        self.validator.adata.uns["default_embedding"] = "X_other"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'X_other' in 'uns['default_embedding']' is not valid, "
                "it must be a key of 'adata.obsm'."
            ],
        )

    def test_X_approximate_distribution_is_str(self):

        """
        X_approximate_distribution str. The value MUST be "count" [...] or "normal"
        """

        self.validator.adata.uns["X_approximate_distribution"] = ["count"]
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: '['count']' in 'uns['X_approximate_distribution']' "
                "is not valid, it must be a string."
            ],
        )

    def test_X_approximate_distribution_is_valid(self):

        """
        X_approximate_distribution str. The value MUST be "count" [...] or "normal"
        """

        self.validator.adata.uns["X_approximate_distribution"] = "COUNT"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'COUNT' in 'uns['X_approximate_distribution']' is "
                "not valid. Allowed terms: ['count', 'normal']."
            ],
        )


class TestObsm(unittest.TestCase):

    """
    Fail cases for adata.obsm
    """

    def setUp(self):

        self.validator = validate.Validator()
        self.validator.adata = examples.adata.copy()

    def test_obsm_values_ara_numpy(self):

        """
        values in obsm MUST be a numpy.ndarray
        """

        self.validator.adata.obsm["X_tsne"] = pd.DataFrame(
            self.validator.adata.obsm["X_umap"], index=self.validator.adata.obs_names
        )
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: All embeddings have to be of 'numpy.ndarray' type, "
                "'adata.obsm['X_tsne']' is <class 'pandas.core.frame.DataFrame'>')"
            ],
        )

    def test_obsm_values_at_least_one_X(self):

        """
        At least one key for the embedding MUST be prefixed with "X_"
        """

        self.validator.adata.obsm["umap"] = self.validator.adata.obsm["X_umap"]
        self.validator.adata.uns["default_embedding"] = "umap"
        del self.validator.adata.obsm["X_umap"]
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: At least one embedding in 'obsm' has to have a "
                "key with an 'X_' prefix."
            ],
        )

    def test_obsm_shape(self):

        """
        Curators MUST annotate one or more two-dimensional (m >= 2) embeddings
        """

        # Makes 1 column array
        self.validator.adata.obsm["X_umap"] = numpy.delete(
            self.validator.adata.obsm["X_umap"], 0, 1
        )
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: All embeddings must have as many rows as cells, and "
                "at least two columns.'adata.obsm['X_umap']' has shape "
                "of '(2, 1)'"
            ],
        )


class TestAddingLabels(unittest.TestCase):

    """
    Tests the addition of labels from IDs based on schema specification. The test is done by comparing manually
    created dataframes (positive control) agains the ones produced by the validator
    """

    @classmethod
    def setUpClass(cls):

        # Manually created  data (positive control)
        cls.adata_with_labels = examples.adata_with_labels

        # Validate test data
        validator = validate.Validator()
        validator.adata = examples.adata.copy()
        validator.validate_adata()

        # Add labels through validator
        cls.label_writer = validate.LabelWriter(validator)
        cls.label_writer._add_labels()

    def test_var_added_labels(self):

        """
        When a dataset is uploaded, cellxgene Data Portal MUST automatically add the matching human-readable
        name for the corresponding feature identifier and the inferred NCBITaxon term for the reference organism
        to the var dataframe. Curators MUST NOT annotate the following columns:

            - feature_name. If the feature_biotype is "gene" then this MUST be the human-readable ENSEMBL gene
            name assigned to the feature_id. If the feature_biotype is "spike-in" then this MUST be the
            ERCC Spike-In identifier appended with " spike-in control".
            - feature_reference. This MUST be the reference organism for a feature:
                Homo sapiens	"NCBITaxon:9606"
                Mus musculus	"NCBITaxon:10090"
                SARS-CoV-2	"NCBITaxon:2697049"
                ERCC Spike-Ins	"NCBITaxon:32630"
        """

        for column in ["feature_name", "feature_reference"]:
            expected_column = self.adata_with_labels.var[column]
            obtained_column = self.label_writer.adata.var[column]

            for i, j in zip(expected_column.tolist(), obtained_column.tolist()):
                self.assertEqual(i, j)

    def test_obs_added_labels(self):

        """
        When a dataset is uploaded, the cellxgene Data Portal MUST automatically add the matching human-readable
        name for the corresponding ontology term to the obs dataframe.
        Curators MUST NOT annotate the following columns.

            - assay. categorical with str categories. This MUST be the human-readable name assigned to the value
            of assay_ontology_term_id. Any clarifying text enclosed in parentheses and appended to
            assay_ontology_term_id MUST be appended to assay.
            - cell_type. categorical with str categories. This MUST be the human-readable name assigned to the value
            of cell_type_ontology_term_id.
            - development_stage. categorical with str categories. This MUST be "unknown" if set in
            development_stage_ontology_term_id; otherwise, this MUST be the human-readable name assigned to
            the value of development_stage_ontology_term_id.
            - disease. categorical with str categories. This MUST be the human-readable name assigned to
            the value of disease_ontology_term_id.
            - ethnicity. categorical with str categories. This MUST be "na" or "unknown" if
            set in ethnicity_ontology_term_id; otherwise, this MUST be the human-readable
            name assigned to the value of ethnicity_ontology_term_id.
            - organism. categorical with str categories. This MUST be the human-readable name assigned
            to the value of organism_ontology_term_id.
            - sex. categorical with str categories. This MUST be "unknown" if set in sex_ontology_term_id;
            otherwise, this MUST be the human-readable name assigned to the value of sex_ontology_term_id.
            - tissue. categorical with str categories. This MUST be the human-readable name assigned to the
            value of tissue_ontology_term_id. " (cell culture)" or " (organoid)" MUST
            be appended if present in tissue_ontology_term_id.
        """

        for column in [
            "assay",
            "cell_type",
            "development_stage",
            "disease",
            "ethnicity",
            "organism",
            "sex",
            "tissue",
        ]:
            expected_column = self.adata_with_labels.obs[column]
            obtained_column = self.label_writer.adata.obs[column]

            for i, j in zip(expected_column.tolist(), obtained_column.tolist()):
                self.assertEqual(i, j)
