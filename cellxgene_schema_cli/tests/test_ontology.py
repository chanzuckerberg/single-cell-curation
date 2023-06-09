import unittest
from unittest.mock import patch

import pytest
from cellxgene_schema import ontology
from fixtures.examples_ontology_test import (
    invalid_genes,
    invalid_ontologies,
    invalid_species,
    invalid_terms,
    valid_genes,
    valid_ontologies,
    valid_terms,
)

# Tests for internal functions of the OntologyChecker and GeneChecker classes


class TestGeneChecker(unittest.TestCase):
    def setUp(self):
        self.valid_species = ontology.SupportedOrganisms
        self.invalid_species = invalid_species
        self.valid_genes = valid_genes
        self.invalid_genes = invalid_genes

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
        self.valid_ontologies = valid_ontologies
        self.invalid_ontologies = invalid_ontologies
        self.valid_terms = valid_terms
        self.invalid_terms = invalid_terms

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

                self.assertTrue(self.ontologyChecker.is_valid_term_id(ontology_id, term_id))
                self.assertIsNone(self.ontologyChecker.assert_term_id(ontology_id, term_id))
                self.assertEqual(
                    self.ontologyChecker.get_term_label(ontology_id, term_id),
                    term_label,
                )

    def test_invalid_term_ids(self):
        for ontology_id in self.invalid_terms:
            for term_id in self.invalid_terms[ontology_id]:
                self.assertFalse(self.ontologyChecker.is_valid_term_id(ontology_id, term_id))
                with self.assertRaises(ValueError):
                    self.ontologyChecker.assert_term_id(ontology_id, term_id)


@pytest.fixture
def organisms():
    return ["apple", "dog", "mouse"]


def test_get_deprecated_features(tmp_path, organisms):
    expected_deprecated_feature_ids = []
    for organism in organisms:
        with open(f"{tmp_path}/{organism}_diff.txt", "w") as fp:
            organism_feature_ids = [f"{organism}:{i}" for i in range(4)]
            for feature_id in organism_feature_ids:
                fp.write(feature_id + "\n")
            expected_deprecated_feature_ids.extend(organism_feature_ids)
    with patch("cellxgene_schema.utils.env.ONTOLOGY_DIR", tmp_path):
        actual_deprecated_features = ontology.get_deprecated_feature_ids()
    expected_deprecated_feature_ids.sort()
    actual_deprecated_features.sort()
    assert expected_deprecated_feature_ids == actual_deprecated_features


def test_get_deprecated_features__no_files(tmp_path):
    with patch("cellxgene_schema.utils.env.ONTOLOGY_DIR", tmp_path):
        actual_deprecated_features = ontology.get_deprecated_feature_ids()
    assert actual_deprecated_features == []


def test_get_deprecated_features__empty_feature_files(tmp_path, organisms):
    for organism in organisms:
        with open(f"{tmp_path}/{organism}_diff.txt", "w") as fp:
            fp.write("")
    with patch("cellxgene_schema.utils.env.ONTOLOGY_DIR", tmp_path):
        actual_deprecated_features = ontology.get_deprecated_feature_ids()
    assert actual_deprecated_features == []
