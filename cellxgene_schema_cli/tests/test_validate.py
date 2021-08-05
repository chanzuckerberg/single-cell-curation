import unittest
import pandas as pd
from cellxgene_schema import ontology
from cellxgene_schema import validate
import numpy
import fixtures.examples_validate as examples
from scipy import sparse


class TestFieldValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema_def = validate._get_schema_definition(examples.SCHEMA_VERSION)
        cls.OntologyChecker = ontology.OntologyChecker()

    def setUp(self):
        self.validator = validate.Validator()
        self.validator.adata = examples.adata_empty
        self.column_name = "cell_type_ontology_term_id"
        self.column_schema = self.schema_def["components"]["obs"]["columns"][
            self.column_name
        ]
        self.curie_constraints = self.schema_def["components"]["obs"]["columns"][
            self.column_name
        ]["curie_constraints"]

    def test_schema_defintion(self):
        """
        Tests that the definition of schema is correct
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


class TestSparsity(unittest.TestCase):

    def setUp(self):
        self.validator = validate.Validator()
        self.validator.adata = examples.adata.copy()
        self.validator._set_schema_def()

    def test_sparsity(self):

        # If a matrix is sparse and not a sparse csr matrix there should be a warning
        self.validator.adata.X = numpy.zeros(self.validator.adata.X.shape)
        self.validator._validate_sparsity()
        self.assertTrue(self.validator.warnings)

        self.validator.warnings = []
        self.validator.adata.X = sparse.csc_matrix(self.validator.adata.X)
        self.validator._validate_sparsity()
        self.assertTrue(self.validator.warnings)

        # Correct sparsity
        self.validator.warnings = []
        self.validator.adata.X = sparse.csr_matrix(self.validator.adata.X)
        self.validator._validate_sparsity()
        self.assertFalse(self.validator.warnings)


class TestObsmValidation(unittest.TestCase):

    def setUp(self):
        self.validator = validate.Validator()
        self.validator.adata = examples.adata.copy()
        self.schema_def = validate._get_schema_definition(examples.SCHEMA_VERSION)
        self.good_uns = examples.good_uns

    def test_validate_good_obsm(self):
        self.validator._validate_embedding_dict()
        self.assertFalse(self.validator.errors)

    def test_validate_bad_obsm(self):

        # Wrong dimension, delete one column
        key = list(self.validator.adata.obsm.keys())[0]
        self.validator.adata.obsm[key] = numpy.delete(self.validator.adata.obsm[key], 0, axis=1)

        self.validator._validate_embedding_dict()
        self.assertTrue(self.validator.errors)

    def test_validate_bad_obsm_2(self):

        # Wrong type
        key = list(self.validator.adata.obsm.keys())[0]
        self.validator.adata.obsm[key] = pd.DataFrame(self.validator.adata.obsm[key],
                                                      index = self.validator.adata.obs_names)

        self.validator._validate_embedding_dict()
        self.assertTrue(self.validator.errors)


class TestUnsValidation(unittest.TestCase):

    def setUp(self):
        self.validator = validate.Validator()
        self.validator.adata = examples.adata
        self.schema_def = validate._get_schema_definition(examples.SCHEMA_VERSION)
        self.good_uns = examples.good_uns
        self.bad_uns = examples.bad_uns

    def test_validate_good_uns(self):
        self.validator._validate_dict(
            self.good_uns, "uns", self.schema_def["components"]["uns"]
        )
        self.assertFalse(self.validator.errors)

    def test_validate_bad_uns(self):

        # Do one key at a time but skip schema_version as that's not checked here
        for key, value in self.bad_uns.items():

            if key == "schema_version":
                continue

            current_dict = self.good_uns.copy()
            current_dict[key] = value

            self.validator.errors = []
            self.validator._validate_dict(
                current_dict, "uns", self.schema_def["components"]["uns"]
            )
            self.assertTrue(self.validator.errors)


class TestColumnValidation(unittest.TestCase):
    def setUp(self):
        self.validator = validate.Validator()
        self.validator.adata = examples.adata_empty
        self.good_var = examples.good_var
        self.good_obs = examples.good_obs
        self.bad_var = examples.bad_var
        self.bad_obs = examples.bad_obs
        self.schema_def = validate._get_schema_definition(examples.SCHEMA_VERSION)

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

        self.validator._validate_column(
            self.unique.index, "index", "unique_df", self.column_def_uniq
        )
        self.assertFalse(self.validator.errors)

        self.validator._validate_column(
            self.unique["col1"], "col1", "unique_df", self.column_def_uniq
        )
        self.assertFalse(self.validator.errors)

    def test_validate_unique_valid_duped(self):

        self.validator._validate_column(
            self.duped["col1"], "col1", "duped_df", self.column_def_not_uniq
        )
        self.assertFalse(self.validator.errors)

    def test_validate_unique_invalid_duped(self):

        self.validator._validate_column(
            self.duped.index, "index", "duped_df", self.column_def_uniq
        )
        self.assertTrue(self.validator.errors)

    def test_validate_unique_invalid_duped2(self):

        self.validator._validate_column(
            self.duped["col1"], "col1", "duped_df", self.column_def_uniq
        )
        self.assertTrue(self.validator.errors)

    def test_ontology_columns(self):


        # Correct example. This only tests the column validation process and therefore the tests excludes those columns
        # that have dependencies with other columns and need the entire dataframe for validation
        dfs = [self.good_var, self.good_obs]
        components = ["var", "obs"]

        for good_df, component in zip(dfs, components):
            for column in good_df.columns:
                columns_def = self.schema_def["components"][component]["columns"]
                if "dependencies" not in columns_def[column]:
                    self.validator.errors = []  # Reset errors
                    self.validator._validate_column(
                        good_df[column], column, "obs", columns_def[column]
                    )
                    self.assertFalse(self.validator.errors)

        # Bad columns, do each individually
        dfs = [self.bad_var, self.bad_obs]
        components = ["var", "obs"]
        for bad_df, component in zip(dfs, components):
            for column in bad_df.columns:
                columns_def = self.schema_def["components"][component]["columns"]
                if "dependencies" not in columns_def[column]:
                    self.validator.errors = []  # Reset errors
                    self.validator._validate_column(
                        bad_df[column], column, "obs", columns_def[column]
                    )
                    self.assertTrue(self.validator.errors)



class TestH5adValidation(unittest.TestCase):
    def setUp(self):
        self.h5ad_valid = examples.h5ad_valid
        self.invalid_files = examples.h5ad_invalid
        self.validator = validate.Validator()

    def test_validate(self):

        # Valid h5ad
        self.assertTrue(self.validator.validate_adata(self.h5ad_valid))

        # Invalid h5ads
        for i in self.invalid_files:
            self.assertFalse(self.validator.validate_adata(i))


class TestAddLabelFunctions(unittest.TestCase):
    def setUp(self):

        # Set up test data
        self.test_adata = examples.adata
        self.test_adata_with_labels = examples.adata_with_labels
        self.schema_def = validate._get_schema_definition(examples.SCHEMA_VERSION)

        validator = validate.Validator()
        validator.adata = self.test_adata
        validator.validate_adata()
        self.writer = validate.LabelWriter(validator)

    def test_get_dictionary_mapping_feature_id(self):

        # Good
        ids = ["ERCC-00002", "ENSG00000127603", "ENSMUSG00000059552", "ENSSASG00005000004"]
        labels = ["ERCC-00002 spike-in control", "MACF1", "Trp53", "S"]
        expected_dict = {i: j for i, j in zip(ids, labels)}
        self.assertEqual(
            self.writer._get_mapping_dict_feature_id(ids), expected_dict
        )

        # Bad
        ids = ["NO_GENE"]
        expected_dict = {i: j for i, j in zip(ids, labels)}
        with self.assertRaises(KeyError):
            self.writer._get_mapping_dict_feature_id(ids)

    def test_get_dictionary_mapping_feature_reference(self):

        # Good
        ids = ["ERCC-00002", "ENSG00000127603", "ENSMUSG00000059552", "ENSSASG00005000004"]
        labels = ["NCBITaxon:32630", "NCBITaxon:9606", "NCBITaxon:10090", "NCBITaxon:2697049"]
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

    def test_get_new_labels(self):

        # Test getting a column with labels based on ids for adata.obs
        for component in ["obs", "var"]:
            for column, column_definition in self.schema_def["components"][component][
                "columns"
            ].items():
                if "add_labels" in column_definition:
                    for label_def in column_definition["add_labels"]:
                        expected_column = self.test_adata_with_labels.obs[label_def["to"]]
                        obtained_column = self.writer._get_labels(
                            component, column, column_definition, label_def["type"]
                        )
                        for i, j in zip(expected_column.tolist(), obtained_column.tolist()):
                            self.assertEqual(i, j)

    def test_get_new_adata(self):

        # Test getting a column with labels based on ids
        expected_adata = self.test_adata_with_labels
        self.writer._add_labels()
        obtained_adata = self.writer.adata
        self.assertTrue(all(expected_adata.obs == obtained_adata.obs))
