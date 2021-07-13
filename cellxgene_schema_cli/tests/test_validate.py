import unittest
import pandas as pd
import os
import anndata
from cellxgene_schema import ontology
from cellxgene_schema import validate

SCHEMA_VERSION = "2.0.0"
FIXTURES_ROOT = os.path.join(os.path.dirname(__file__), "fixtures")


class TestFieldValidation(unittest.TestCase):

    def setUp(self):
        self.schema_def = validate.get_schema_definition(SCHEMA_VERSION)
        self.ontologyChecker = ontology.ontologyChecker()

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
                    self.assertTrue(self.ontologyChecker.is_valid_ontology(ontology))

    def test_cell_type_ontology(self):
        column_name = "cell_type_ontology_term_id"
        column_schema = self.schema_def["components"]["obs"]["columns"][column_name]
        curie_constraints = self.schema_def["components"]["obs"]["columns"][column_name]["curie_constraints"]

        self.assertEqual(column_schema["type"], "curie")
        self.assertEqual(curie_constraints["ontologies"], ["CL"])

        # Good curies should be an empty list of errors
        errors = validate._validate_curie("CL:0000066", column_name, curie_constraints)
        self.assertFalse(errors)

        errors= validate._validate_curie("CL:0000192", column_name, curie_constraints)
        self.assertFalse(errors)

        # Bad curies should be a non-empty list of errors
        errors= validate._validate_curie("EFO:0009899", column_name, curie_constraints)
        self.assertTrue(errors)

        errors= validate._validate_curie("NO_TERM2", column_name, curie_constraints)
        self.assertTrue(errors)


class TestColumnValidation(unittest.TestCase):

    def test_validate_unique(self):
        unique = pd.DataFrame(
            [["abc", "def"], ["ghi", "jkl"], ["mnop", "qrs"]],
            index=["X", "Y", "Z"],
            columns=["col1", "col2"],
        )
        duped = pd.DataFrame(
            [["abc", "def"], ["ghi", "qrs"], ["abc", "qrs"]],
            index=["X", "Y", "X"],
            columns=["col1", "col2"],
        )

        schema_def = {"unique": True}

        errors = validate._validate_column(
            unique.index, "index", "unique_df", schema_def
        )
        self.assertFalse(errors)

        errors = validate._validate_column(duped.index, "index", "duped_df", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("is not unique", errors[0])

        errors = validate._validate_column(
            unique["col1"], "col1", "unique_df", schema_def
        )
        self.assertFalse(errors)

        errors = validate._validate_column(
            duped["col1"], "col1", "duped_df", schema_def
        )
        self.assertEqual(len(errors), 1)
        self.assertIn("is not unique", errors[0])

        schema_def = {"unique": False}
        errors = validate._validate_column(
            duped["col1"], "col1", "duped_df", schema_def
        )
        self.assertFalse(errors)

    def test_ontology_columns(self):

        good_df = pd.DataFrame(
            [
                [
                    "CL:0000066",
                ],
                [
                    "CL:0000192"
                ]
            ],
            index=["X", "Y"],
            columns=["cell_type_ontology_term_id"]
        )

        bad_df = pd.DataFrame(
            [
                [
                    "EFO:0009899", #Bad term of CL
                ],
                [
                    "NO_TERM"
                ]
            ],
            index=["X", "Y"],
            columns=["cell_type_ontology_term_id"]
        )

        columns_def = validate.get_schema_definition(SCHEMA_VERSION)["components"]["obs"]["columns"]
        for column in good_df.columns:
            errors = validate._validate_column(good_df[column], column, "obs", columns_def[column],)
            self.assertFalse(errors)


class TestH5adValidation(unittest.TestCase):

    def setUp(self):
        self.anndata_valid = anndata.read(os.path.join(FIXTURES_ROOT, "h5ads", "example_valid.h5ad"))
        self.anndata_invalid_CL = anndata.read(os.path.join(FIXTURES_ROOT, "h5ads", "example_invalid_CL.h5ad"))
        self.schema_def = validate.get_schema_definition(SCHEMA_VERSION)


    def test_validate(self):

        # Good
        self.assertTrue(validate.validate_adata(self.anndata_valid, self.schema_def))

        # Bad
        self.assertFalse(validate.validate_adata(self.anndata_invalid_CL, self.schema_def))


class TestAddLabelFunctions(unittest.TestCase):

    def setUp(self):

        # Set up test data with cell_type as example
        X = pd.DataFrame(
            [[0],
             [0]],
            index=["X", "Y"],
        )

        obs = pd.DataFrame(
            [["CL:0000066"],
             ["CL:0000192"]],
            index=["X", "Y"],
            columns=["cell_type_ontology_term_id"]
        )

        self.cell_type_column = pd.Categorical(["epithelial cell", "smooth muscle cell"])
        self.test_adata = anndata.AnnData(X=X, obs=obs)
        self.test_component = "obs"
        self.test_column = "cell_type_ontology_term_id"
        self.schema_def = {
            "title": "Corpora schema version 2.0.0",
            "type": "anndata",
            "components": {
                "obs": {
                    "columns": {
                        "cell_type_ontology_term_id": {
                            "type": "curie",
                            "curie_constraints": {"ontologies": ["CL"]},
                            "add_labels": {"type": "curie", "to": "cell_type"}
                        }
                    }
                }
            }
        }

        self.test_column_definition_cell_type = self.schema_def["components"]["obs"]["columns"]["cell_type_ontology_term_id"]
        self.test_adata_with_labels = self.test_adata.copy()
        self.test_adata_with_labels.obs["cell_type"] = self.cell_type_column

    def test_get_dictionary_mapping(self):
        # Good
        ids = ["CL:0000066", "CL:0000192"]
        labels = ["epithelial cell", "smooth muscle cell"]
        ontologies = ["CL"]
        expected_dict = {i: j for i,j in zip(ids, labels)}
        self.assertEqual(validate._get_mapping_dict_curie(ids, ontologies), expected_dict)

        ids = ["EFO:0009899", "EFO:0009922"]
        labels = ["10x 3' v2", "10x 3' v3"]
        ontologies = ["EFO"]
        expected_dict = {i: j for i,j in zip(ids, labels)}
        self.assertEqual(validate._get_mapping_dict_curie(ids, ontologies), expected_dict)

        ids = ["MONDO:0100096"]
        labels = ["COVID-19"]
        ontologies = ["MONDO"]
        expected_dict = {i: j for i,j in zip(ids, labels)}
        self.assertEqual(validate._get_mapping_dict_curie(ids, ontologies), expected_dict)

        ids = ["MmusDv:0000062", "HsapDv:0000174"]
        labels = ["2 month-old stage", "1 month-old human stage"]
        ontologies = ["MmusDv", "HsapDv"]
        expected_dict = {i: j for i,j in zip(ids, labels)}
        self.assertEqual(validate._get_mapping_dict_curie(ids, ontologies), expected_dict)

        # Bad
        ids = ["CL:0000066", "CL:0000192asdf"]
        ontologies = ["UBERON"]
        with self.assertRaises(ValueError):
            validate._get_mapping_dict_curie(ids, ontologies)

        ids = ["CL:0000066", "CL:0000192", "UBERON:0002048"]
        ontologies = ["UBERON"]
        with self.assertRaises(ValueError):
            validate._get_mapping_dict_curie(ids, ontologies)

        ids = ["CL:NO_TERM"]
        ontologies = ["CL"]
        with self.assertRaises(ValueError):
            validate._get_mapping_dict_curie(ids, ontologies)

    def test_get_new_labels(self):

        # Test getting a column with labels based on ids
        expected_column = self.cell_type_column
        obtained_column = validate._get_labels(self.test_adata, self.test_component, self.test_column, self.test_column_definition_cell_type)
        for i,j in zip(expected_column.tolist(), obtained_column.tolist()):
            self.assertEqual(i, j)

    def test_get_new_adata(self):

        # Test getting a column with labels based on ids
        expected_adata = self.test_adata_with_labels
        obtained_adata = validate._add_labels(self.test_adata, self.schema_def)
        self.assertTrue(all(expected_adata.obs == obtained_adata.obs))
