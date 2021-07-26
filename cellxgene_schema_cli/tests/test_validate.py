import unittest
import pandas as pd
import os
import anndata
from cellxgene_schema import ontology
from cellxgene_schema import validate

SCHEMA_VERSION = "2.0.0"
FIXTURES_ROOT = os.path.join(os.path.dirname(__file__), "fixtures")


class ExampleData:

    good_obs = pd.DataFrame(
        [
            [
                "CL:0000066", "EFO:0009899", "MONDO:0100096", "NCBITaxon:9606", "PATO:0000383", "UBERON:0002048", True
            ],
            [
                "CL:0000192", "EFO:0010183 (sci-plex)", "PATO:0000461", "NCBITaxon:10090", "unknown",
                "CL:0000192 (cell culture)", False
            ]
        ],
        index=["X", "Y"],
        columns=["cell_type_ontology_term_id", "assay_ontology_term_id", "disease_ontology_term_id",
                 "organism_ontology_term_id", "sex_ontology_term_id", "tissue_ontology_term_id", "is_primary_data"]
    )

    obs_expected = pd.DataFrame(
        [
            ["epithelial cell", "10x 3' v2", "COVID-19", "Homo sapiens", "female", "lung"],
            ["smooth muscle cell", "single cell library construction (sci-plex)", "normal", "Mus musculus", "unknown",
             "smooth muscle cell (cell culture)"]
        ],
        index=["X", "Y"],
        columns=["cell_type", "assay", "disease", "organism", "sex", "tissue"]
    )

    bad_obs = pd.DataFrame(
        [
            [
                "CL:NO_TERM", "EFO:NO_TERM", "MONDO:NO_TERM", "NCBITaxon:NO_TERM", "PATO:NO_TERM", "UBERON:NO_TERM",
                "True"
            ],
            [
                "CL:0000182", "EFO:00212", "MONDO:0324", "NCBITaxon:00324", "PATO:2003", "UBERON:3203", "False"
            ]
        ],
        index=["X", "Y"],

        columns=["cell_type_ontology_term_id", "assay_ontology_term_id", "disease_ontology_term_id",
                 "organism_ontology_term_id", "sex_ontology_term_id", "tissue_ontology_term_id", "is_primary_data"]
    )

    good_uns = {"schema_version": SCHEMA_VERSION}

    X = pd.DataFrame(
        [[0] * good_obs.shape[1],
         [0] * good_obs.shape[1]
         ],
        index=["X", "Y"],
    )

    adata = anndata.AnnData(X=X, obs=good_obs, uns=good_uns)
    adata_with_labels = anndata.AnnData(X=X, obs=pd.concat([good_obs, obs_expected], axis=1), uns=good_uns)
    adata_empty = anndata.AnnData(X=X, uns=good_uns)


class TestFieldValidation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.schema_def = validate._get_schema_definition(SCHEMA_VERSION)
        cls.OntologyChecker = ontology.OntologyChecker()

    def setUp(self):
        self.validator = validate.Validator()
        self.validator.adata = ExampleData.adata_empty
        self.column_name = "cell_type_ontology_term_id"
        self.column_schema = self.schema_def["components"]["obs"]["columns"][self.column_name]
        self.curie_constraints = self.schema_def["components"]["obs"]["columns"][self.column_name]["curie_constraints"]

    def test_schema_defintion(self):
        """
        Tests that the definition of schema is correct
        """

        self.assertIsInstance(self.schema_def["components"], dict)
        self.assertIsInstance(self.schema_def["components"]["obs"], dict)
        self.assertIsInstance(self.schema_def["components"]["obs"]["columns"], dict)

        # Check that any columns in obs that are "curie" have "curie_constraints" and "ontologies" under the constraints
        for i in self.schema_def["components"]["obs"]["columns"]:
            self.assertTrue("type" in self.schema_def["components"]["obs"]["columns"][i])
            if i == "curie":
                self.assertIsInstance(self.schema_def["components"]["obs"]["columns"][i]["curie_constrains"], dict)
                self.assertIsInstance(self.schema_def["components"]["obs"]["columns"][i]["curie_constrains"]["ontolgies"], list)

                # Check that the allowed ontologies are in the ontology checker
                for ontology in self.schema_def["components"]["obs"]["columns"][i]["curie_constrains"]["ontolgies"]:
                    self.assertTrue(self.OntologyChecker.is_valid_ontology(ontology))

    def test_validate_ontology_good(self):
        self.validator._validate_curie("CL:0000066", self.column_name, self.curie_constraints)
        self.validator._validate_curie("CL:0000192", self.column_name, self.curie_constraints)
        self.assertFalse(self.validator.errors)

    def test_validate_ontology_wrong_ontology(self):
        self.validator._validate_curie("EFO:0009899", self.column_name, self.curie_constraints)
        self.assertTrue(self.validator.errors)

    def test_validate_ontology_wrong_term(self):
        self.validator._validate_curie("NO_TERM2", self.column_name, self.curie_constraints)
        self.assertTrue(self.validator.errors)


class TestColumnValidation(unittest.TestCase):

    def setUp(self):
        self.validator = validate.Validator()
        self.validator.adata = ExampleData.adata_empty

        self.unique = pd.DataFrame(
            [["abc", "def"], ["ghi", "jkl"], ["mnop", "qrs"]],
            index=["X", "Y", "Z"],
            columns=["col1", "col2"],
        )

        self.duped = pd.DataFrame(
            [["abc", "def"], ["ghi", "qrs"], ["abc", "qrs"]],
            index=["X", "Y", "X"],
            columns=["col1", "col2"],
        )

        self.column_def_uniq = {"unique": True}
        self.column_def_not_uniq = {"unique": False}

    def test_validate_unique_good(self):

        self.validator._validate_column(self.unique.index, "index", "unique_df", self.column_def_uniq)
        self.assertFalse(self.validator.errors)

        self.validator._validate_column(self.unique["col1"], "col1", "unique_df", self.column_def_uniq)
        self.assertFalse(self.validator.errors)

    def test_validate_unique_valid_duped(self):

        self.validator._validate_column(self.duped["col1"], "col1", "duped_df", self.column_def_not_uniq)
        self.assertFalse(self.validator.errors)

    def test_validate_unique_invalid_duped(self):

        self.validator._validate_column(self.duped.index, "index", "duped_df", self.column_def_uniq)
        self.assertTrue(self.validator.errors)

    def test_validate_unique_invalid_duped2(self):

        self.validator._validate_column(self.duped["col1"], "col1", "duped_df", self.column_def_uniq)
        self.assertTrue(self.validator.errors)


    def test_ontology_columns(self):

        columns_def = validate._get_schema_definition(SCHEMA_VERSION)["components"]["obs"]["columns"]

        # Correct example
        good_df = ExampleData.good_obs

        for column in good_df.columns:
            self.validator.errors = [] # Reset errors
            self.validator._validate_column(good_df[column], column, "obs", columns_def[column],)
            self.assertFalse(self.validator.errors)

        # Bad columns, do each individually
        bad_df = ExampleData.bad_obs
        for column in bad_df.columns:
            self.validator.errors = [] # Reset errors
            self.validator._validate_column(bad_df[column], column, "obs", columns_def[column],)
            self.assertTrue(self.validator.errors)

class TestH5adValidation(unittest.TestCase):

    def test_validate(self):

        h5ad_dir = os.path.join(FIXTURES_ROOT, "h5ads" )
        h5ad_valid = os.path.join(h5ad_dir, "example_valid.h5ad")

        validator = validate.Validator()

        # Good
        self.assertTrue(validator.validate_adata(h5ad_valid))

        # Test all bad cases
        invalid_files = ["example_invalid_CL.h5ad", "example_invalid_assay.h5ad", "example_invalid_disease.h5ad",
                         "example_invalid_organism.h5ad", "example_invalid_primary_data.h5ad",
                         "example_invalid_sex.h5ad", "example_invalid_tissue.h5ad"]
        for i in invalid_files:
            self.assertFalse(validator.validate_adata(os.path.join(h5ad_dir, i)))


class TestAddLabelFunctions(unittest.TestCase):

    def setUp(self):

        # Set up test data
        self.test_adata = ExampleData.adata
        self.test_adata_with_labels = ExampleData.adata_with_labels
        self.schema_def = validate._get_schema_definition(SCHEMA_VERSION)

        validator = validate.Validator()
        validator.adata = self.test_adata
        validator.validate_adata()
        self.writer = validate.LabelWriter(validator)

    def test_get_dictionary_mapping(self):

        # Good
        ids = ["CL:0000066", "CL:0000192"]
        labels = ["epithelial cell", "smooth muscle cell"]
        curie_constraints = self.schema_def["components"]["obs"]["columns"]["cell_type_ontology_term_id"]["curie_constraints"]
        expected_dict = {i: j for i,j in zip(ids, labels)}
        self.assertEqual(validate._get_mapping_dict_curie(ids, curie_constraints), expected_dict)

        ids = ["EFO:0009899", "EFO:0009922"]
        labels = ["10x 3' v2", "10x 3' v3"]
        curie_constraints = self.schema_def["components"]["obs"]["columns"]["assay_ontology_term_id"]["curie_constraints"]
        expected_dict = {i: j for i,j in zip(ids, labels)}
        self.assertEqual(validate._get_mapping_dict_curie(ids, curie_constraints), expected_dict)

        ids = ["MONDO:0100096"]
        labels = ["COVID-19"]
        curie_constraints = self.schema_def["components"]["obs"]["columns"]["disease_ontology_term_id"]["curie_constraints"]
        expected_dict = {i: j for i,j in zip(ids, labels)}
        self.assertEqual(validate._get_mapping_dict_curie(ids, curie_constraints), expected_dict)

        #ids = ["MmusDv:0000062", "HsapDv:0000174"]
        #labels = ["2 month-old stage", "1 month-old human stage"]
        #curie_constraints = self.schema_def["components"]["obs"]["columns"]["development_stage_ontology_term_id"]
        #expected_dict = {i: j for i,j in zip(ids, labels)}
        #self.assertEqual(validate._get_mapping_dict_curie(ids, curie_constraints), expected_dict)

        # Bad
        curie_constraints = self.schema_def["components"]["obs"]["columns"]["cell_type_ontology_term_id"]["curie_constraints"]

        ids = ["CL:0000066", "CL:0000192_FOO"]
        with self.assertRaises(ValueError):
            validate._get_mapping_dict_curie(ids, curie_constraints)

        ids = ["CL:0000066", "CL:0000192", "UBERON:0002048"]
        with self.assertRaises(ValueError):
            validate._get_mapping_dict_curie(ids, curie_constraints)

        ids = ["CL:NO_TERM"]
        with self.assertRaises(ValueError):
            validate._get_mapping_dict_curie(ids, curie_constraints)

    def test_get_new_labels(self):

        # Test getting a column with labels based on ids for adata.obs
        component = "obs"
        for column, column_definition in self.schema_def["components"]["obs"]["columns"].items():
            if "add_labels" in column_definition:
                expected_column = self.test_adata_with_labels.obs[column_definition["add_labels"]["to"]]
                obtained_column = self.writer._get_labels(component, column, column_definition)
                for i, j in zip(expected_column.tolist(), obtained_column.tolist()):
                    self.assertEqual(i, j)

    def test_get_new_adata(self):

        # Test getting a column with labels based on ids
        expected_adata = self.test_adata_with_labels
        self.writer._add_labels()
        obtained_adata = self.writer.adata
        self.assertTrue(all(expected_adata.obs == obtained_adata.obs))
