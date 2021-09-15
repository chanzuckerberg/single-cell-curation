import unittest

from cellxgene_schema import ontology
from cellxgene_schema import validate
from cellxgene_schema import schema
import fixtures.examples_validate as examples

# Tests for internal functions of the Validator and LabelWriter classes.


class TestFieldValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema_def = schema.get_schema_definition(examples.SCHEMA_VERSION)
        cls.OntologyChecker = ontology.OntologyChecker()

    def setUp(self):
        self.validator = validate.Validator()
        self.validator.adata = examples.adata_minimal
        self.column_name = "cell_type_ontology_term_id"
        self.column_schema = self.schema_def["components"]["obs"]["columns"][
            self.column_name
        ]
        self.curie_constraints = self.schema_def["components"]["obs"]["columns"][
            self.column_name
        ]["curie_constraints"]

    def test_schema_definition(self):

        """
        Tests that the definition of schema is well-defined
        """

        self.assertIsInstance(self.schema_def["components"], dict)
        self.assertIsInstance(self.schema_def["components"]["obs"], dict)
        self.assertIsInstance(self.schema_def["components"]["obs"]["columns"], dict)

        # Check that any columns in obs that are "curie" have "curie_constraints" and "ontologies" under the constraints
        for i in self.schema_def["components"]["obs"]["columns"]:
            self.assertTrue(
                "type" in self.schema_def["components"]["obs"]["columns"][i]
            )
            if i == "curie":
                self.assertIsInstance(
                    self.schema_def["components"]["obs"]["columns"][i][
                        "curie_constrains"
                    ],
                    dict,
                )
                self.assertIsInstance(
                    self.schema_def["components"]["obs"]["columns"][i][
                        "curie_constrains"
                    ]["ontolgies"],
                    list,
                )

                # Check that the allowed ontologies are in the ontology checker
                for ontology_name in self.schema_def["components"]["obs"]["columns"][i][
                    "curie_constrains"
                ]["ontolgies"]:
                    self.assertTrue(
                        self.OntologyChecker.is_valid_ontology(ontology_name)
                    )

    def test_validate_ontology_good(self):
        self.validator._validate_curie(
            "CL:0000066", self.column_name, self.curie_constraints
        )
        self.validator._validate_curie(
            "CL:0000192", self.column_name, self.curie_constraints
        )
        self.assertFalse(self.validator.errors)

    def test_validate_ontology_wrong_ontology(self):
        self.validator._validate_curie(
            "EFO:0009899", self.column_name, self.curie_constraints
        )
        self.assertTrue(self.validator.errors)

    def test_validate_ontology_wrong_term(self):
        self.validator._validate_curie(
            "NO_TERM2", self.column_name, self.curie_constraints
        )
        self.assertTrue(self.validator.errors)


class TestAddLabelFunctions(unittest.TestCase):
    def setUp(self):

        # Set up test data
        self.test_adata = examples.adata
        self.test_adata_with_labels = examples.adata_with_labels
        self.schema_def = schema.get_schema_definition(examples.SCHEMA_VERSION)

        validator = validate.Validator()
        validator.adata = self.test_adata
        validator.validate_adata()
        self.writer = validate.AnnDataLabelAppender(validator)

    def test_get_dictionary_mapping_feature_id(self):

        # Good
        ids = [
            "ERCC-00002",
            "ENSG00000127603",
            "ENSMUSG00000059552",
            "ENSSASG00005000004",
        ]
        labels = [
            "ERCC-00002 (spike-in control)",
            "MACF1",
            "Trp53",
            "S_ENSSASG00005000004",
        ]
        expected_dict = {i: j for i, j in zip(ids, labels)}
        self.assertEqual(self.writer._get_mapping_dict_feature_id(ids), expected_dict)

        # Bad
        ids = ["NO_GENE"]
        expected_dict = {i: j for i, j in zip(ids, labels)}
        with self.assertRaises(KeyError):
            self.writer._get_mapping_dict_feature_id(ids)

    def test_get_dictionary_mapping_feature_reference(self):

        # Good
        ids = [
            "ERCC-00002",
            "ENSG00000127603",
            "ENSMUSG00000059552",
            "ENSSASG00005000004",
        ]
        labels = [
            "NCBITaxon:32630",
            "NCBITaxon:9606",
            "NCBITaxon:10090",
            "NCBITaxon:2697049",
        ]
        expected_dict = {i: j for i, j in zip(ids, labels)}
        self.assertEqual(
            self.writer._get_mapping_dict_feature_reference(ids), expected_dict
        )

        # Bad
        ids = ["NO_GENE"]
        expected_dict = {i: j for i, j in zip(ids, labels)}
        with self.assertRaises(KeyError):
            self.writer._get_mapping_dict_feature_id(ids)

    def test_get_dictionary_mapping_curie(self):

        # Good
        ids = ["CL:0000066", "CL:0000192"]
        labels = ["epithelial cell", "smooth muscle cell"]
        curie_constraints = self.schema_def["components"]["obs"]["columns"][
            "cell_type_ontology_term_id"
        ]["curie_constraints"]
        expected_dict = {i: j for i, j in zip(ids, labels)}
        self.assertEqual(
            self.writer._get_mapping_dict_curie(ids, curie_constraints), expected_dict
        )

        ids = ["EFO:0009899", "EFO:0009922"]
        labels = ["10x 3' v2", "10x 3' v3"]
        curie_constraints = self.schema_def["components"]["obs"]["columns"][
            "assay_ontology_term_id"
        ]["curie_constraints"]
        expected_dict = {i: j for i, j in zip(ids, labels)}
        self.assertEqual(
            self.writer._get_mapping_dict_curie(ids, curie_constraints), expected_dict
        )

        ids = ["MONDO:0100096"]
        labels = ["COVID-19"]
        curie_constraints = self.schema_def["components"]["obs"]["columns"][
            "disease_ontology_term_id"
        ]["curie_constraints"]
        expected_dict = {i: j for i, j in zip(ids, labels)}
        self.assertEqual(
            self.writer._get_mapping_dict_curie(ids, curie_constraints), expected_dict
        )

        # Bad
        curie_constraints = self.schema_def["components"]["obs"]["columns"][
            "cell_type_ontology_term_id"
        ]["curie_constraints"]

        ids = ["CL:0000066", "CL:0000192_FOO"]
        with self.assertRaises(ValueError):
            self.writer._get_mapping_dict_curie(ids, curie_constraints)

        ids = ["CL:0000066", "CL:0000192", "UBERON:0002048"]
        with self.assertRaises(ValueError):
            self.writer._get_mapping_dict_curie(ids, curie_constraints)

        ids = ["CL:NO_TERM"]
        with self.assertRaises(ValueError):
            self.writer._get_mapping_dict_curie(ids, curie_constraints)
