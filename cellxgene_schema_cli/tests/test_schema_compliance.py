import unittest
from unittest.mock import patch

import fixtures.examples_validate as examples
import numpy
import pandas as pd
from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender

# Tests for schema compliance of an AnnData object


class TestValidAnndata(unittest.TestCase):

    """
    Tests a valid AnnData object. Most other tests below modify this AnnData object and test for failure cases.

    The valid AnnData object has all valid cases described in the schema.
    """

    def setUp(self):
        self.validator = Validator()
        self.validator.adata = examples.adata.copy()

    def test_valid_anndata(self):
        self.validator.validate_adata()
        self.assertFalse(self.validator.errors)


class TestH5adValidation(unittest.TestCase):

    """
    Checks that validation from h5ad works, only does one invalid example as extensive testing is done in the classes
    below
    """

    def setUp(self):
        self.h5ad_valid_file = examples.h5ad_valid
        self.h5ad_invalid_file = examples.h5ad_invalid
        self.validator = Validator()

    def test_validate(self):
        # Valid h5ad
        self.assertTrue(self.validator.validate_adata(self.h5ad_valid_file))
        # Invalid h5ads
        self.assertFalse(self.validator.validate_adata(self.h5ad_invalid_file))


class BaseValidationTest(unittest.TestCase):
    def setUp(self):
        self.validator = Validator()
        self.validator.adata = examples.adata.copy()

        # Override the schema definition here
        self.validator._set_schema_def()

        # lower threshold for low gene count warning
        self.validator.schema_def["components"]["var"]["warn_if_less_than_rows"] = 1


class TestExpressionMatrix(BaseValidationTest):

    """
    Fail cases for expression matrices (anndata.X and anndata.raw.X)
    """

    def test_shapes(self):
        """
        All matrix layers MUST have the same shape, and have the same cell labels and gene labels.
        """

        # Creates a raw layer
        self.validator.adata.raw = self.validator.adata
        self.validator.adata.raw.var.drop("feature_is_filtered", axis=1, inplace=True)
        self.validator.adata.X = examples.adata_non_raw.X.copy()

        # remove one gene
        self.validator.adata = self.validator.adata[:, 1:]
        self.validator.validate_adata()
        self.assertIn(
            "ERROR: Number of genes in X (3) is different than raw.X (4).",
            self.validator.errors,
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
                "WARNING: Sparsity of 'X' is 0.875 which is greater than 0.5, "
                "and it is not a 'scipy.sparse.csr_matrix'. It is "
                "STRONGLY RECOMMENDED to use this type of matrix for "
                "the given sparsity."
            ],
        )

    def test_raw_values(self):
        """
        When both `adata.X` and `adata.raw.X` are present, but `adata.raw.X` contains non-integer values an error
        is raised.
        """

        self.validator.adata = examples.adata_no_raw_values.copy()
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: Raw data may be missing: data in 'raw.X' contains non-integer values."],
        )

    def test_raw_existence(self):
        """
        Except for ATAC-seq and methylation data, raw data is REQUIRED
        """

        # ATAC - raw layer not required
        # The assignment above makes X to not be raw: self.validator.adata.uns["X_normalization"] = "CPM"
        # The following line makes it to be scATAC-seq data (EFO:0010891)
        # Missing raw data in atac-seq data is allowed, thus the following should not return an error message
        self.validator.errors = []
        self.validator.adata.obs["assay_ontology_term_id"] = "EFO:0010891"
        self.validator.adata.obs["suspension_type"] = "nucleus"
        self.validator.adata.obs.loc[:, ["suspension_type"]] = self.validator.adata.obs.astype("category")
        self.validator.validate_adata()
        self.assertEqual(self.validator.errors, [])

    def test_final_strongly_recommended(self):
        """
        Except for ATAC-seq and methylation data, final matrix is STRONGLY RECOMMENDED
        """

        # move raw to X amd: i.e. there is no final
        self.validator.adata.X = self.validator.adata.raw.X
        del self.validator.adata.raw
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.warnings,
            [
                "WARNING: Only raw data was found, i.e. there is no 'raw.X'. "
                "It is STRONGLY RECOMMENDED that 'final' (normalized) data is provided."
            ],
        )


class TestObs(BaseValidationTest):

    """
    Fail cases in adata.uns
    """

    def test_column_presence(self):
        """
        obs is a pandas.DataFrame. Curators MUST annotate the following columns in the obs dataframe.
        """

        columns = [
            "development_stage_ontology_term_id",
            "disease_ontology_term_id",
            "self_reported_ethnicity_ontology_term_id",
            "is_primary_data",
            "sex_ontology_term_id",
            "tissue_ontology_term_id",
            "donor_id",
            "suspension_type",
        ]

        for column in columns:
            with self.subTest(column=column):
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
                "ERROR: Dataframe 'obs' is missing column " "'organism_ontology_term_id'.",
                "ERROR: Checking values with dependencies failed for "
                "adata.obs['self_reported_ethnicity_ontology_term_id'], this is likely due "
                "to missing dependent column in adata.obs.",
                "ERROR: Checking values with dependencies failed for "
                "adata.obs['development_stage_ontology_term_id'], this is likely due "
                "to missing dependent column in adata.obs.",
            ],
        )

    def test_column_presence_assay(self):
        """
        obs is a pandas.DataFrame. Curators MUST annotate the following columns in the obs dataframe.

        A separate check is need for assay_ontology_term_id because removing from anndata results in multiple
        errors given that other columns depend on its presence
        """

        self.validator.adata = examples.adata.copy()
        self.validator.adata.obs.drop("assay_ontology_term_id", axis=1, inplace=True)
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: Dataframe 'obs' is missing column " "'assay_ontology_term_id'.",
                "ERROR: Checking values with dependencies failed for "
                "adata.obs['suspension_type'], this is likely due "
                "to missing dependent column in adata.obs.",
            ],
        )

    def test_obsolete_term_id(self):
        """
        Terms documented as obsolete in an ontology MUST NOT be used. For example, EFO:0009310
        for obsolete_10x v2 was marked as obsolete in EFO version 3.31.0 and replaced by
        EFO:0009899 for 10x 3' v2.

        https://www.ebi.ac.uk/ols/ontologies/efo/terms?short_form=EFO_0009310
        """

        # Not a valid term
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "assay_ontology_term_id"] = "EFO:0009310"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: 'EFO:0009310' in 'assay_ontology_term_id' is a deprecated term id of 'EFO'."],
        )

    def test_assay_ontology_term_id(self):
        """
        assay_ontology_term_id categorical with str categories.
        This MUST be an EFO term and either child of "EFO:0002772" or "EFO:0010183"
        """

        # Not a valid term
        # self.validator.adata.obs["assay_ontology_term_id"][0] = "CL:000001"
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "assay_ontology_term_id"] = "CL:000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: 'CL:000001' in 'assay_ontology_term_id' is not a valid " "ontology term id of 'EFO'."],
        )

        # Not a valid child
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "assay_ontology_term_id"] = "EFO:0000001"
        self.validator.errors = []
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'assay_ontology_term_id' is not a "
                "child term id of '[['EFO:0002772', 'EFO:0010183']]'."
            ],
        )

        # Includes extraneous text
        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0], "assay_ontology_term_id"
        ] = "EFO:0010183 (sci-plex)"
        self.validator.errors = []
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: 'EFO:0010183 (sci-plex)' in 'assay_ontology_term_id' is not a valid ontology term id of 'EFO'."],
        )

    def test_cell_type_ontology_term_id(self):
        """
        cell_type_ontology_term_id categorical with str categories. This MUST be a CL term.
        """

        # Not a valid term
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "cell_type_ontology_term_id"] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: 'EFO:0000001' in 'cell_type_ontology_term_id' is not a valid " "ontology term id of 'CL'."],
        )

    def test_development_stage_ontology_term_id_human(self):
        """
        development_stage_ontology_term_id categorical with str categories. If unavailable, this MUST be "unknown".
        If organism_ontolology_term_id is "NCBITaxon:9606" for Homo sapiens,
        this MUST be the most accurate HsapDv term.
        """

        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "organism_ontology_term_id"] = "NCBITaxon:9606"
        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0], "development_stage_ontology_term_id"
        ] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'development_stage_ontology_term_id' is "
                "not a valid ontology term id of 'HsapDv'. When 'organism_ontology_term_id' is 'NCBITaxon:9606' "
                "(Homo sapiens), 'development_stage_ontology_term_id' MUST be a term id of 'HsapDv' or unknown."
            ],
        )

    def test_development_stage_ontology_term_id_mouse(self):
        """
        If organism_ontolology_term_id is "NCBITaxon:10090" for Mus musculus,
        this MUST be the most accurate MmusDv term
        """

        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "organism_ontology_term_id"] = "NCBITaxon:10090"
        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0], "development_stage_ontology_term_id"
        ] = "EFO:0000001"
        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "na"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'development_stage_ontology_term_id' is "
                "not a valid ontology term id of 'MmusDv'. When 'organism_ontology_term_id' is 'NCBITaxon:10090' "
                "(Mus musculus), 'development_stage_ontology_term_id' MUST be a term id of 'MmusDv' or unknown."
            ],
        )

    def test_development_stage_ontology_term_id_all_species(self):
        """
        All other it MUST be children of UBERON:0000105 and not UBERON:0000071
        """

        # Fail case not an UBERON term
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "organism_ontology_term_id"] = "NCBITaxon:10114"
        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0], "development_stage_ontology_term_id"
        ] = "EFO:0000001"
        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "na"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'development_stage_ontology_term_id' is "
                "not a valid ontology term id of 'UBERON'. When 'organism_ontology_term_id' is not 'NCBITaxon:10090' "
                "nor 'NCBITaxon:9606', 'development_stage_ontology_term_id' MUST be a child term id of "
                "'UBERON:0000105' excluding 'UBERON:0000071', or unknown."
            ],
        )

        # All other it MUST be children of UBERON:0000105 and not UBERON:0000071
        # Fail case UBERON:0000071
        self.validator.errors = []
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "organism_ontology_term_id"] = "NCBITaxon:10114"
        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0], "development_stage_ontology_term_id"
        ] = "UBERON:0000071"
        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "na"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'UBERON:0000071' in 'development_stage_ontology_term_id' is not allowed'. When "
                "'organism_ontology_term_id' is not 'NCBITaxon:10090' "
                "nor 'NCBITaxon:9606', 'development_stage_ontology_term_id' MUST be a child term id of "
                "'UBERON:0000105' excluding 'UBERON:0000071', or unknown.",
            ],
        )

    def test_disease_ontology_term_id(self):
        """
        disease_ontology_term_id categorical with str categories. This MUST be a MONDO term or
        PATO:0000461 for normal or healthy.
        """

        # Invalid ontology
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "disease_ontology_term_id"] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'disease_ontology_term_id' is not a "
                "valid ontology term id of 'MONDO, PATO'. Only 'PATO:0000461' is allowed for 'PATO' term ids."
            ],
        )

        # Invalid PATO term id
        self.validator.errors = []
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "disease_ontology_term_id"] = "PATO:0001894"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'PATO:0001894' in 'disease_ontology_term_id' is not an allowed term: '[['PATO:0000461']]'. "
                "Only 'PATO:0000461' is allowed for 'PATO' term ids."
            ],
        )

    def test_self_reported_ethnicity_ontology_term_id(self):
        """
        self_reported_ethnicity_ontology_term_id categorical with str categories.
        If organism_ontolology_term_id is "NCBITaxon:9606" for Homo sapiens,
        this MUST be either a HANCESTRO term, "multiethnic", or "unknown" if unavailable.
        Otherwise, for all other organisms this MUST be "na".
        """

        # If organism_ontolology_term_id is "NCBITaxon:9606" for Homo sapiens,
        # this MUST be either a HANCESTRO term, "multiethnic", or "unknown" if unavailable.
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "organism_ontology_term_id"] = "NCBITaxon:9606"
        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'self_reported_ethnicity_ontology_term_id' is "
                "not a valid ontology term id of 'HANCESTRO'. When 'organism_ontology_term_id' is 'NCBITaxon:9606' "
                "(Homo sapiens), self_reported_ethnicity_ontology_term_id MUST be either: a term id of 'HANCESTRO', "
                "'multiethnic' if more than one ethnicity is reported, or 'unknown' if unavailable."
            ],
        )

        # Otherwise, for all other organisms this MUST be "na". Below is the test case for mouse data.
        # development_stage_ontology_term_id has to be set to an appropriate mouse term id, otherwise there
        # will be an error in that field.
        self.validator.errors = []
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "organism_ontology_term_id"] = "NCBITaxon:10090"
        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0], "development_stage_ontology_term_id"
        ] = "MmusDv:0000003"
        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'self_reported_ethnicity_ontology_term_id' is not a "
                "valid value of 'self_reported_ethnicity_ontology_term_id'. When 'organism_ontology_term_id' is NOT "
                "'NCBITaxon:9606' (Homo sapiens), self_reported_ethnicity_ontology_term_id MUST be 'na'."
            ],
        )

    def test_organism_ontology_term_id(self):
        """
        organism_ontology_term_id categorical with str categories. This MUST be a child of NCBITaxon:33208.
        """

        # Setting "organism_ontology_term_id" to "EFO:0000001" is the fail case. However since this represents neither
        # human nor mouse, then two other columns that are dependent on it need to be set appropriately to avoid
        # other error messages: "development_stage_ontology_term_id" and "self_reported_ethnicity_ontology_term_id"
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "organism_ontology_term_id"] = "EFO:0000001"
        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0], "development_stage_ontology_term_id"
        ] = "unknown"
        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0],
            "self_reported_ethnicity_ontology_term_id",
        ] = "na"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'organism_ontology_term_id' is not a valid "
                "ontology term id of 'NCBITaxon'. Only children term ids of 'NCBITaxon:33208' for metazoan are allowed."
            ],
        )

    def test_tissue_ontology_term_id_base(self):
        """
        tissue_ontology_term_id categorical with str categories. This MUST be the term that best describes the tissue
        that this cell was derived from, depending on the type of biological sample:
        """

        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "tissue_ontology_term_id"] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: 'EFO:0000001' in 'tissue_ontology_term_id' is not a " "valid ontology term id of 'UBERON, CL'."],
        )

    def test_tissue_ontology_term_id_cell_culture(self):
        """
        Cell Culture - MUST be a CL term appended with " (cell culture)"
        """

        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0], "tissue_ontology_term_id"
        ] = "CL:0000057 (CELL culture)"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'CL:0000057 (CELL culture)' in 'tissue_ontology_term_id' is "
                "not a valid ontology term id of 'UBERON, CL'."
            ],
        )

    def test_tissue_ontology_term_id_organoid(self):
        """
        Organoid - MUST be an UBERON term appended with " (organoid)"
        """

        self.validator.adata.obs.loc[
            self.validator.adata.obs.index[0], "tissue_ontology_term_id"
        ] = "CL:0000057 (ORGANOID)"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'CL:0000057 (ORGANOID)' in 'tissue_ontology_term_id' is "
                "not a valid ontology term id of 'UBERON, CL'."
            ],
        )

    def test_sex_ontology_term_id(self):
        """
        sex_ontology_term_id categorical with str categories.
        This MUST be a child of PATOPATO:0001894 for phenotypic sex or "unknown" if unavailable
        """

        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "sex_ontology_term_id"] = "EFO:0000001"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'EFO:0000001' in 'sex_ontology_term_id' is "
                "not a valid ontology term id of 'PATO'. Only 'PATO:0000383', 'PATO:0000384', 'PATO:0001340', "
                "or 'unknown' are allowed."
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
            ["ERROR: Column 'is_primary_data' in dataframe 'obs' " "must be boolean, not 'object'."],
        )

    def test_donor_id_must_be_categorical(self):
        """
        donor_id categorical with str categories. This MUST be free-text that identifies
        a unique individual that data were derived from. It is STRONGLY RECOMMENDED
        that this identifier be designed so that it is unique to:
        - a given individual within the collection of datasets that includes this dataset
        - a given individual across all collections in the cellxgene Data Portal
        """

        self.validator.adata.obs["donor_id"] = "NA"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: Column 'donor_id' in dataframe 'obs' " "must be categorical, not object."],
        )

    def test_donor_id_must_not_be_empty(self):
        self.validator.adata.obs["donor_id"] = self.validator.adata.obs["donor_id"].cat.add_categories("")
        self.validator.adata.obs["donor_id"].iloc[0] = ""
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: Column 'donor_id' in dataframe 'obs' " "must not contain empty values."],
        )

    def test_donor_id_must_not_be_nan(self):
        self.validator.adata.obs["donor_id"][0] = numpy.nan
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: Column 'donor_id' in dataframe 'obs' " "must not contain NaN values."],
        )

    def test_suspension_type(self):
        """
        suspension_id categorical with str categories. This field MUST be "cell", "nucleus", or "na". The allowed
        values depend on the assay_ontology_term_id. MUST  fail if the corresponding assay is present in the table, but
        the value of the suspension_type does not match the required value(s) in the table.
        """
        match_assays = {
            "EFO:0010010": ["cell", "nucleus"],
            "EFO:0008720": ["nucleus"],
            "EFO:0008722": ["cell", "nucleus"],
            "EFO:0030002": ["cell"],
            "EFO:0008853": ["cell"],
            "EFO:0030026": ["nucleus"],
            "EFO:0010550": ["cell", "nucleus"],
            "EFO:0008919": ["cell"],
            "EFO:0008939": ["nucleus"],
            "EFO:0030027": ["nucleus"],
        }

        for assay, suspension_types in match_assays.items():
            with self.subTest(assay=assay):
                # Resetting validator
                self.validator.adata = examples.adata.copy()
                self.validator.errors = []
                self.validator.warnings = []

                invalid_suspension_type = "na"
                self.validator.adata.obs.loc[
                    self.validator.adata.obs.index[1], "suspension_type"
                ] = invalid_suspension_type
                self.validator.adata.obs.loc[self.validator.adata.obs.index[1], "assay_ontology_term_id"] = assay
                self.validator.validate_adata()
                self.assertEqual(
                    self.validator.errors,
                    [
                        f"ERROR: Column 'suspension_type' in dataframe 'obs' contains invalid values "
                        f"'['{invalid_suspension_type}']'. Values must be one of {suspension_types} when "
                        f"'assay_ontology_term_id' is {assay}"
                    ],
                )

    def test_suspension_type_ancestors_inclusive(self):
        """
        suspension_id categorical with str categories. This field MUST be "cell", "nucleus", or "na". The allowed
        values depend on the assay_ontology_term_id. MUST  fail if the corresponding assay is present in the table, but
        the value of the suspension_type does not match the required value(s) in the table.
        """
        match_assays_or_children = {
            "EFO:0030080": ["cell", "nucleus"],
            "EFO:0007045": ["nucleus"],
            "EFO:0009294": ["cell"],
            "EFO:0010184": ["cell", "nucleus"],
            "EFO:0009918": ["na"],
            "EFO:0700000": ["na"],
            "EFO:0030005": ["na"],
        }
        for assay, suspension_types in match_assays_or_children.items():
            with self.subTest(assay=assay):
                # Resetting validator
                self.validator.adata = examples.adata.copy()
                self.validator.errors = []
                self.validator.warnings = []

                invalid_suspension_type = "na"
                if assay in {"EFO:0009918", "EFO:0700000", "EFO:0030005"}:
                    invalid_suspension_type = "nucleus"
                    self.validator.adata.obs["suspension_type"] = self.validator.adata.obs[
                        "suspension_type"
                    ].cat.remove_unused_categories()
                self.validator.adata.obs.loc[self.validator.adata.obs.index[1], "assay_ontology_term_id"] = assay
                self.validator.adata.obs.loc[
                    self.validator.adata.obs.index[1], "suspension_type"
                ] = invalid_suspension_type
                self.validator.validate_adata()
                self.assertEqual(
                    self.validator.errors,
                    [
                        f"ERROR: Column 'suspension_type' in dataframe 'obs' contains invalid values "
                        f"'['{invalid_suspension_type}']'. Values must be one of {suspension_types} when "
                        f"'assay_ontology_term_id' is {assay} or its children"
                    ],
                )

    def test_suspension_type_with_child_term_id(self):
        """
        suspension_id categorical with str categories. This field MUST be "cell", "nucleus", or "na". The allowed
        values depend on the assay_ontology_term_id. MUST support matching against ancestor term rules if specified.
        """
        with self.subTest("failure"):
            self.validator.adata.obs.loc[
                self.validator.adata.obs.index[0], "assay_ontology_term_id"
            ] = "EFO:0030008"  # child of EFO:0009294
            self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "suspension_type"] = "nucleus"

            self.validator.validate_adata()
            self.assertEqual(
                self.validator.errors,
                [
                    "ERROR: Column 'suspension_type' in dataframe 'obs' contains invalid values "
                    "'['nucleus']'. Values must be one of ['cell'] when "
                    "'assay_ontology_term_id' is EFO:0009294 or its children"
                ],
            )

        with self.subTest("success"):
            self.validator.adata.obs.loc[
                self.validator.adata.obs.index[0], "assay_ontology_term_id"
            ] = "EFO:0008904"  # child of EFO:0007045
            self.validator.adata.obs["suspension_type"][0] = "nucleus"

            self.validator.validate_adata()
            self.assertEqual(
                self.validator.errors,
                [],
            )

    def test_suspension_type_unrecognized_assay(self):
        """
        suspension_id categorical with str categories. This field MUST be "cell", "nucleus", or "na". The allowed
        values depend on the assay_ontology_term_id. MUST warn if the corresponding assay is not recognized.
        """
        self.validator.adata.obs.loc[self.validator.adata.obs.index[1], "assay_ontology_term_id"] = "EFO:0010183"
        self.validator.validate_adata()
        with self.subTest("no errors"):
            self.assertEqual(self.validator.errors, [])

        with self.subTest("warnings"):
            self.assertEqual(
                self.validator.warnings,
                [
                    "WARNING: Data contains assay(s) that are not represented in the 'suspension_type' schema "
                    "definition table. Ensure you have selected the most appropriate value for the assay(s) between "
                    "'cell', 'nucleus', and 'na'. Please contact cellxgene@chanzuckerberg.com "
                    "during submission so that the assay(s) can be added to the schema definition document."
                ],
            )

    def test_categories_with_zero_values_warn(self):
        modified_donor_id = self.validator.adata.obs["donor_id"].cat.add_categories("donor_3")
        self.validator.adata.obs["donor_id"] = modified_donor_id
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.warnings,
            [
                "WARNING: Column 'donor_id' in dataframe 'obs' "
                "contains a category 'donor_3' with zero observations. "
                "These categories will be removed when `--add-labels` "
                "flag is present."
            ],
        )

    def test_deprecated_fields(self):
        self.validator.adata.obs["ethnicity"] = "test"
        self.validator.adata.obs["ethnicity_ontology_term_id"] = "unknown"

        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: The field 'ethnicity' is present in 'obs', but it is deprecated.",
                "ERROR: The field 'ethnicity_ontology_term_id' is present in 'obs', but it is deprecated.",
            ],
        )

    def test_fields_with_double_underscore_fail(self):
        self.validator.adata.obs["__test_field"] = "test"

        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: The field '__test_field' in 'obs' is invalid. Fields that start with '__' are reserved.",
            ],
        )

    def test_nan_values_must_be_rejected(self):
        """
        NaN values should not be allowed in dataframes
        """
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "tissue_ontology_term_id"] = numpy.nan
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: Column 'tissue_ontology_term_id' in dataframe 'obs' must not contain NaN values."],
        )


class TestVar(BaseValidationTest):

    """
    Fail cases in adata.var and adata.raw.var
    """

    def test_var_and_raw_var_same_index(self):
        """
        var.index MUST contain unique identifiers for features. raw.var.index MUST be identical to var.index.
        """

        # Swap first row for second one
        var = Validator.getattr_anndata(self.validator.adata, "var")

        # First swap the index
        new_index = list(var.index)
        tmp = new_index[0]
        new_index[0] = new_index[1]
        new_index[1] = tmp
        var.set_index(pd.Index(new_index), inplace=True)

        # Then swap the actual rows
        tmp = var.iloc[0, :].copy()
        var.iloc[0, :] = var.iloc[1, :].copy()
        var.iloc[1, :] = tmp

        self.validator.validate_adata()
        print("FOO", self.validator.errors)
        self.assertEqual(
            self.validator.errors,
            ["ERROR: Index of 'raw.var' is not identical to index of 'var'."],
        )

    def test_check_unique_var(self):
        """
        var.index MUST contain unique ENSEMBL gene identifiers for features.
        """

        for component_name in ["var", "raw.var"]:
            with self.subTest(component_name=component_name):
                # Resetting validator
                self.validator.adata = examples.adata.copy()
                self.validator.errors = []

                # Duplicate 1st row in var and assign it to 2nd
                component = Validator.getattr_anndata(self.validator.adata, component_name)
                new_index = list(component.index)
                new_index[1] = new_index[0]
                component.set_index(pd.Index(new_index), inplace=True)
                component.iloc[1, :] = component.iloc[0, :]

                self.validator.validate_adata()
                self.assertEqual(
                    self.validator.errors,
                    [f"ERROR: Column 'index' in dataframe '{component_name}' is not unique."],
                )

    def test_column_presence(self):
        """
        var is a pandas.DataFrame. Curators MUST annotate the following columns in the obs dataframe.
        feature_is_filtered must not be in raw.var, and it's only checked in var
        """

        columns = ["feature_is_filtered"]

        for component_name in ["var", "raw.var"]:
            for column in columns:
                if column == "feature_is_filtered" and component_name == "raw.var":
                    continue
                with self.subTest(component_name=component_name, column=column):
                    # Resetting validator
                    self.validator.errors = []
                    self.validator.adata = examples.adata.copy()

                    component = Validator.getattr_anndata(self.validator.adata, component_name)
                    component.drop(column, axis=1, inplace=True)

                    self.validator.validate_adata()
                    self.assertEqual(
                        self.validator.errors,
                        [f"ERROR: Dataframe '{component_name}' is missing " f"column '{column}'."],
                    )

    def test_feature_is_filtered(self):
        """
        feature_is_filtered bool. This MUST be True if the feature was filtered out in the final matrix (X)
        but is present in the raw matrix (raw.X). The value for all cells of the given feature in the
        final matrix MUST be 0.

        Otherwise, this MUST be False.
        """

        # Duplicate 1st row in var and assigned to 2nd
        self.validator.adata.var["feature_is_filtered"][0] = True
        for i in range(self.validator.adata.X.shape[0]):
            self.validator.adata.X[i, 0] = 0
        self.validator.adata.X[0, 0] = 1

        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: Some features are 'True' in 'feature_is_filtered' of dataframe 'var', "
                "but there are 1 non-zero values in the corresponding columns of the matrix 'X'. "
                "All values for these features must be 0."
            ],
        )

        # Test that feature_is_filtered is a bool and not a string
        self.validator.adata.var["feature_is_filtered"] = "string"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: Column 'feature_is_filtered' in dataframe 'var' must be boolean, not 'object'."],
        )

    def test_columns_not_in_raw_var(self):
        """
        Curators MUST annotate the following column only in the var dataframe.
        This column MUST NOT be present in raw.var:
            feature_is_filtered
        """

        self.validator.adata.raw = self.validator.adata
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: Column 'feature_is_filtered' must not be present in 'raw.var'."],
        )

    def test_feature_id_wrong_format(self):
        """
        feature_id (var.index) str.

        This tests the case of an ID with an incorrect format "ENSEBML_NOGENE"
        """

        for component_name in ["var", "raw.var"]:
            with self.subTest(component_name=component_name):
                # Resetting validator
                self.validator.adata = examples.adata.copy()
                self.validator.errors = []

                component = Validator.getattr_anndata(self.validator.adata, component_name)

                new_index = list(component.index)
                new_index[0] = "ENSEBML_NOGENE"
                component.set_index(pd.Index(new_index), inplace=True)

                self.validator.validate_adata()
                self.assertEqual(
                    self.validator.errors,
                    [
                        f"ERROR: Could not infer organism from feature ID 'ENSEBML_NOGENE' "
                        f"in '{component_name}', make sure it is a valid ID."
                    ],
                )

    def test_feature_id_non_existent_ensembl(self):
        """
        feature_id (var.index) str.

        This tests the case of an ENSEMBL ID that has the right format but doesn't exist
        """

        for component_name in ["var", "raw.var"]:
            with self.subTest(component_name=component_name):
                # Resetting validator
                self.validator.adata = examples.adata.copy()
                self.validator.errors = []

                component = Validator.getattr_anndata(self.validator.adata, component_name)

                new_index = list(component.index)
                new_index[0] = "ENSG000"
                component.set_index(pd.Index(new_index), inplace=True)

                self.validator.validate_adata()
                self.assertEqual(
                    self.validator.errors,
                    [f"ERROR: 'ENSG000' is not a valid feature ID in '{component_name}'."],
                )

    def test_feature_id_non_existent_ercc(self):
        """
        feature_id (var.index) str.

        This tests the case of an ERCC ID that has the right format but doesn't exist
        """

        for component_name in ["var", "raw.var"]:
            with self.subTest(component_name=component_name):
                # Resetting validator
                self.validator.adata = examples.adata.copy()
                self.validator.errors = []

                component = Validator.getattr_anndata(self.validator.adata, component_name)

                new_index = list(component.index)
                new_index[0] = "ERCC-000000"
                component.set_index(pd.Index(new_index), inplace=True)

                self.validator.validate_adata()
                self.assertEqual(
                    self.validator.errors,
                    [f"ERROR: 'ERCC-000000' is not a valid feature ID in '{component_name}'."],
                )

    def test_should_warn_for_low_gene_count(self):
        """
        Raise a warning if there are too few genes
        """
        self.validator.schema_def["components"]["var"]["warn_if_less_than_rows"] = 100
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.warnings,
            ["WARNING: Dataframe 'var' only has 4 rows. Features SHOULD NOT be filtered from expression matrix."],
        )


class TestUns(BaseValidationTest):

    """
    Fail cases in adata.uns
    """

    def test_optional_fields_schema_version_is_old(self):
        latest_version = "1.0.0"
        old_version = "0.0.1"
        self.validator.adata.uns["schema_version"] = old_version
        self.validator.schema_version = None
        with self.assertLogs(level="WARNING") as logs, patch(
            "cellxgene_schema.validate.schema.get_schema_versions_supported", return_value=[latest_version]
        ) as mock_supported_versions, patch(
            "cellxgene_schema.validate.schema.get_current_schema_version", return_value=latest_version
        ) as mock_get_current_schema, patch(
            "cellxgene_schema.validate.schema.get_schema_definition"
        ) as mock_get_schema_definition:
            self.validator._set_schema_def()
        self.assertEqual(
            logs.output,
            [
                f"WARNING:cellxgene_schema.validate:Schema version '{old_version}' is not supported. Current supported "
                f"versions: '{latest_version}'. Validating with latest version '{latest_version}'."
            ],
        )
        mock_supported_versions.assert_called_once()
        mock_get_schema_definition.assert_called_with(latest_version)

    def test_optional_fields_schema_version_is_missing(self):
        latest_version = "1.0.0"
        self.validator.schema_version = None
        del self.validator.adata.uns["schema_version"]
        with patch(
            "cellxgene_schema.validate.schema.get_schema_versions_supported", return_value=[latest_version]
        ) as mock_supported_versions, patch(
            "cellxgene_schema.validate.schema.get_current_schema_version", return_value=latest_version
        ) as mock_get_current_schema, patch(
            "cellxgene_schema.validate.schema.get_schema_definition"
        ) as mock_get_schema_definition:
            self.validator._set_schema_def()
        mock_supported_versions.assert_called_once()
        mock_get_current_schema.assert_called_once()
        mock_get_schema_definition.assert_called_with(latest_version)
        self.assertEqual(self.validator.schema_version, latest_version)

    def test_optional_fields_schema_version_is_latest(self):
        latest_version = "1.0.0"
        self.validator.schema_version = None
        self.validator.adata.uns["schema_version"] = latest_version
        with patch(
            "cellxgene_schema.validate.schema.get_schema_versions_supported", return_value=[latest_version]
        ) as mock_supported_versions, patch(
            "cellxgene_schema.validate.schema.get_current_schema_version", return_value=latest_version
        ) as mock_get_current_schema, patch(
            "cellxgene_schema.validate.schema.get_schema_definition"
        ) as mock_get_schema_definition:
            self.validator._set_schema_def()
        mock_supported_versions.assert_called_once()
        mock_get_current_schema.assert_called_once()
        mock_get_schema_definition.assert_called_with(latest_version)
        self.assertEqual(self.validator.schema_version, latest_version)

    def test_optional_fields_schema_version_with_multiple_supported_versions(self):
        latest_version = "1.0.0"
        old_version = "0.0.1"
        self.validator.schema_version = None
        self.validator.adata.uns["schema_version"] = old_version
        with self.assertLogs(level="WARNING") as logs, patch(
            "cellxgene_schema.validate.schema.get_schema_versions_supported", return_value=[old_version, latest_version]
        ) as mock_supported_versions, patch(
            "cellxgene_schema.validate.schema.get_current_schema_version", return_value=latest_version
        ) as mock_get_current_schema, patch(
            "cellxgene_schema.validate.schema.get_schema_definition"
        ) as mock_get_schema_definition:
            self.validator._set_schema_def()
        self.assertEqual(
            logs.output,
            [
                f"WARNING:cellxgene_schema.validate:Schema version '{old_version}' is not supported. Current supported "
                f"versions: '{latest_version}'. Validating with latest version '{latest_version}'."
            ],
        )
        mock_supported_versions.assert_called_once()
        mock_get_current_schema.assert_called_once()
        mock_get_schema_definition.assert_called_with(latest_version)
        self.assertEqual(self.validator.schema_version, latest_version)

    def test_required_fields_title(self):
        """
        Curators MUST annotate `schema_version` and values in uns (title)
        """

        del self.validator.adata.uns["title"]
        self.validator.validate_adata()
        self.assertEqual(self.validator.errors, ["ERROR: 'title' in 'uns' is not present."])

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
            ["ERROR: ' There is a leading space' in 'uns['title']' is not valid, it contains leading spaces."],
        )

        self.validator.adata.uns["title"] = "There is a trailing space "
        self.validator.errors = []
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: 'There is a trailing space ' in 'uns['title']' is not valid, it contains trailing spaces."],
        )

        self.validator.adata.uns["title"] = "There are   double   spaces"
        self.validator.errors = []
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: 'There are   double   spaces' in 'uns['title']' is not valid, it contains double spaces."],
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
            ["ERROR: '['title']' in 'uns['title']' is not valid, " "it must be a string."],
        )

    def test_batch_condition_is_list(self):
        """
        batch_condition list[str]
        """

        # Check valid case of numpy array which is interchangeable with lists
        self.validator.adata.uns["batch_condition"] = numpy.array(self.validator.adata.uns["batch_condition"])
        self.validator.validate_adata()
        self.assertEqual(self.validator.errors, [])

        # Check fail case: not a list nor numpy array
        self.validator.adata.uns["batch_condition"] = "cell_type_ontology_term_id"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: 'cell_type_ontology_term_id' in 'uns['batch_condition']' "
                "is not valid, it must be a list or numpy array."
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
            ["ERROR: Value 'NO_COLUMN' of list 'batch_condition' is not a " "column in 'adata.obs'."],
        )

    def test_default_embedding_is_str(self):
        """
        Default_embedding str.
        """

        self.validator.adata.uns["default_embedding"] = ["X_umap"]
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: '['X_umap']' in 'uns['default_embedding']' is not valid, " "it must be a string."],
        )

    def test_default_embedding_is_key_from_obsm(self):
        """
        Default_embedding str. The value MUST match a key to an embedding in obsm
        """

        self.validator.adata.uns["default_embedding"] = "X_other"
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: 'X_other' in 'uns['default_embedding']' is not valid, " "it must be a key of 'adata.obsm'."],
        )

    def test_X_approximate_distribution_is_str(self):
        """
        X_approximate_distribution str. The value MUST be "count" [...] or "normal".
        Note that `normal` is tested in the happy path test case using `good_uns`.
        """

        # Check valid case of "count" which is not included in valid object
        self.validator.adata.uns["X_approximate_distribution"] = "count"
        self.validator.validate_adata()
        self.assertEqual(self.validator.errors, [])

        # Invalid type: list
        self.validator.adata.uns["X_approximate_distribution"] = ["count"]
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            ["ERROR: '['count']' in 'uns['X_approximate_distribution']' " "is not valid, it must be a string."],
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

    def test_deprecated_fields(self):
        self.validator.adata.uns["X_normalization"] = "test_value"
        self.validator.adata.uns["default_field"] = "test_value"
        self.validator.adata.uns["layer_descriptions"] = "test_value"
        self.validator.adata.uns["tags"] = "test_value"
        self.validator.adata.uns["version"] = "test_value"
        self.validator.adata.uns["contributors"] = "test_value"
        self.validator.adata.uns["preprint_doi"] = "test_value"
        self.validator.adata.uns["project_description"] = "test_value"
        self.validator.adata.uns["project_links"] = "test_value"
        self.validator.adata.uns["project_name"] = "test_value"
        self.validator.adata.uns["publication_doi"] = "test_value"

        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: The field 'X_normalization' is present in 'uns', but it is deprecated.",
                "ERROR: The field 'default_field' is present in 'uns', but it is deprecated.",
                "ERROR: The field 'layer_descriptions' is present in 'uns', but it is deprecated.",
                "ERROR: The field 'tags' is present in 'uns', but it is deprecated.",
                "ERROR: The field 'version' is present in 'uns', but it is deprecated.",
                "ERROR: The field 'contributors' is present in 'uns', but it is deprecated.",
                "ERROR: The field 'preprint_doi' is present in 'uns', but it is deprecated.",
                "ERROR: The field 'project_description' is present in 'uns', but it is deprecated.",
                "ERROR: The field 'project_links' is present in 'uns', but it is deprecated.",
                "ERROR: The field 'project_name' is present in 'uns', but it is deprecated.",
                "ERROR: The field 'publication_doi' is present in 'uns', but it is deprecated.",
            ],
        )


class TestObsm(BaseValidationTest):

    """
    Fail cases for adata.obsm
    """

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
                "'adata.obsm['X_tsne']' is <class 'pandas.core.frame.DataFrame'>')."
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
            ["ERROR: At least one embedding in 'obsm' has to have a " "key with an 'X_' prefix."],
        )

    def test_obsm_shape(self):
        """
        Curators MUST annotate one or more two-dimensional (m >= 2) embeddings
        """

        # Makes 1 column array
        self.validator.adata.obsm["X_umap"] = numpy.delete(self.validator.adata.obsm["X_umap"], 0, 1)
        self.validator.validate_adata()
        self.assertEqual(
            self.validator.errors,
            [
                "ERROR: All embeddings must have as many rows as cells, and "
                "at least two columns.'adata.obsm['X_umap']' has shape "
                "of '(2, 1)'."
            ],
        )


class TestAddingLabels(unittest.TestCase):

    """
    Tests the addition of labels from IDs based on schema specification. The test is done by comparing manually
    created dataframes (positive control) against the ones produced by the validator
    """

    @classmethod
    def setUpClass(cls):
        # Manually created  data (positive control)
        cls.adata_with_labels = examples.adata_with_labels

        # Validate test data
        validator = Validator()
        validator.adata = examples.adata.copy()
        validator.validate_adata()

        # Add labels through validator
        cls.label_writer = AnnDataLabelAppender(validator)
        cls.label_writer._add_labels()

    def test_var_added_labels(self):
        """
        When a dataset is uploaded, cellxgene Data Portal MUST automatically add the matching human-readable
        name for the corresponding feature identifier and the inferred NCBITaxon term for the reference organism
        to the var dataframe. Curators MUST NOT annotate the following columns:

            - feature_name. this MUST be a human-readable ENSEMBL gene name or a ERCC Spike-In identifier
            appended with " spike-in control", corresponding to the feature_id
            - feature_reference. This MUST be the reference organism for a feature:
                Homo sapiens	"NCBITaxon:9606"
                Mus musculus	"NCBITaxon:10090"
                SARS-CoV-2	"NCBITaxon:2697049"
                ERCC Spike-Ins	"NCBITaxon:32630"
            - feature_biotype. This MUST be "gene" if the feature_id is an ENSEMBL gene, or "spike-in" if the feature_id
            is an ERCC Spike-In identifier.
        """

        for column in ["feature_name", "feature_reference", "feature_biotype"]:
            expected_column = self.adata_with_labels.var[column]
            obtained_column = self.label_writer.adata.var[column]

            for i, j in zip(expected_column.tolist(), obtained_column.tolist()):
                with self.subTest(i=i, j=j):
                    self.assertEqual(i, j)

    def test_obs_added_labels(self):
        """
        When a dataset is uploaded, the cellxgene Data Portal MUST automatically add the matching human-readable
        name for the corresponding ontology term to the obs dataframe.
        Curators MUST NOT annotate the following columns.

            - assay. categorical with str categories. This MUST be the human-readable name assigned to the value
            of assay_ontology_term_id.
            - cell_type. categorical with str categories. This MUST be the human-readable name assigned to the value
            of cell_type_ontology_term_id.
            - development_stage. categorical with str categories. This MUST be "unknown" if set in
            development_stage_ontology_term_id; otherwise, this MUST be the human-readable name assigned to
            the value of development_stage_ontology_term_id.
            - disease. categorical with str categories. This MUST be the human-readable name assigned to
            the value of disease_ontology_term_id.
            - self_reported_ethnicity. categorical with str categories. This MUST be "na" or "unknown" if
            set in self_reported_ethnicity_ontology_term_id; otherwise, this MUST be the human-readable
            name assigned to the value of self_reported_ethnicity_ontology_term_id.
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
            "self_reported_ethnicity",
            "organism",
            "sex",
            "tissue",
        ]:
            expected_column = self.adata_with_labels.obs[column]
            obtained_column = self.label_writer.adata.obs[column]

            for i, j in zip(expected_column.tolist(), obtained_column.tolist()):
                with self.subTest(i=i, j=j):
                    self.assertEqual(i, j)

    def test_remove_unused_categories(self):
        modified_donor_id = self.label_writer.adata.obs["donor_id"].cat.add_categories("donor_3")
        self.label_writer.adata.obs["donor_id"] = modified_donor_id

        self.assertCountEqual(
            self.label_writer.adata.obs["donor_id"].dtype.categories,
            ["donor_1", "donor_2", "donor_3"],
        )
        self.label_writer._remove_categories_with_zero_values()

        print(self.label_writer.adata.obs["donor_id"].dtype.categories)

        self.assertCountEqual(
            self.label_writer.adata.obs["donor_id"].dtype.categories,
            ["donor_1", "donor_2"],
        )
