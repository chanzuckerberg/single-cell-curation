import unittest
from cellxgene_schema import ontology
import fixtures.examples_ontology_test as examples

# Tests for internal functions of the OntologyChecker and GeneChecker classes


class TestGeneChecker(unittest.TestCase):
    def setUp(self):
        self.valid_species = ontology.SupportedOrganisms
        self.invalid_species = examples.invalid_species
        self.valid_genes = examples.valid_genes
        self.invalid_genes = examples.invalid_genes

    def test_species_validity(self):
        for species in self.valid_species:
            self.assertIsInstance(ontology.GeneChecker(species), ontology.GeneChecker)
        for species in self.invalid_species:
            with self.assertRaises(ValueError):
                ontology.GeneChecker(species)

    def test_valid_genes(self):
        for species in self.valid_genes:
            geneChecker = ontology.GeneChecker(species)
            for gene_id in self.valid_genes[species]:
                gene_label = self.valid_genes[species][gene_id]

                self.assertTrue(geneChecker.is_valid_id(gene_id))
                self.assertEqual(geneChecker.get_symbol(gene_id), gene_label)

    def test_invalid_genes(self):
        for species in self.invalid_genes:
            geneChecker = ontology.GeneChecker(species)
            for gene_id in self.invalid_genes[species]:

                self.assertFalse(geneChecker.is_valid_id(gene_id))
                with self.assertRaises(ValueError):
                    geneChecker.get_symbol(gene_id)


class TestOntologyChecker(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ontologyChecker = ontology.OntologyChecker()

    def setUp(self):
        self.valid_ontologies = examples.valid_ontologies
        self.invalid_ontologies = examples.invalid_ontologies
        self.valid_terms = examples.valid_terms
        self.invalid_terms = examples.invalid_terms

    def test_ontology_validity(self):
        for i in self.valid_ontologies:
            self.assertTrue(self.ontologyChecker.is_valid_ontology(i))
            self.assertIsNone(self.ontologyChecker.assert_ontology(i))

        for i in self.invalid_ontologies:
            self.assertFalse(self.ontologyChecker.is_valid_ontology(i))
            with self.assertRaises(ValueError):
                self.ontologyChecker.assert_ontology(i)

    def test_valid_term_id(self):
        for ontology_id in self.valid_terms:
            for term_id in self.valid_terms[ontology_id]:
                term_label = self.valid_terms[ontology_id][term_id]

                self.assertTrue(
                    self.ontologyChecker.is_valid_term_id(ontology_id, term_id)
                )
                self.assertIsNone(
                    self.ontologyChecker.assert_term_id(ontology_id, term_id)
                )
                self.assertEqual(
                    self.ontologyChecker.get_term_label(ontology_id, term_id),
                    term_label,
                )

    def test_invalid_term_ids(self):
        for ontology_id in self.invalid_terms:
            for term_id in self.invalid_terms[ontology_id]:
                self.assertFalse(
                    self.ontologyChecker.is_valid_term_id(ontology_id, term_id)
                )
                with self.assertRaises(ValueError):
                    self.ontologyChecker.assert_term_id(ontology_id, term_id)
