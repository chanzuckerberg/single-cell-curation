import unittest

from cellxgene_schema import ontology


class TestGeneChecker(unittest.TestCase):
    def setUp(self):
        self.valid_species = ["Homo sapiens", "Mus musculus"]
        self.invalid_species = ["Caenorhabditis elegans"]

        self.valid_genes = {
            "Homo sapiens": {
                "ENSG00000141510": "TP53"
            },
            "Mus musculus": {
                "ENSMUSG00000059552": "Trp53"
            }
        }

        self.invalid_genes = {
            "Homo sapiens": ["ENSMUSG00000059552", "GENE"],
            "Mus musculus": ["ENSG00000141510", "GENE"]

        }

    def test_species_validity(self):
        for species in self.valid_species:
            self.assertIsInstance(ontology.geneChecker(species), ontology.geneChecker)
        for species in self.invalid_species:
            with self.assertRaises(ValueError):
                ontology.geneChecker(species)

    def test_valid_genes(self):
        for species in self.valid_genes:
            geneChecker = ontology.geneChecker(species)
            for gene_id in self.valid_genes[species]:
                gene_label = self.valid_genes[species][gene_id]

                self.assertTrue(geneChecker.is_valid_id(gene_id))
                self.assertEqual(geneChecker.get_symbol(gene_id), gene_label)

    def test_invalid_genes(self):
        for species in self.invalid_genes:
            geneChecker = ontology.geneChecker(species)
            for gene_id in self.invalid_genes[species]:

                self.assertFalse(geneChecker.is_valid_id(gene_id))
                with self.assertRaises(ValueError):
                    geneChecker.get_symbol(gene_id)


class TestOntologyChecker(unittest.TestCase):

    def setUp(self):
        self.ontologyChecker = ontology.ontologyChecker()

        self.valid_ontologies = ["CL", "EFO", "MONDO", "MmusDv", "NCBITaxon", "UBERON", "PATO"]
        self.invalid_ontologies = ["NOT_ONTOLOGY", "OOO"]

        self.valid_terms = {
            "CL": {
                "CL:0000066": "epithelial cell",
                "CL:0000192": "smooth muscle cell"
            },
            "EFO": {
                "EFO:0009899": "10x 3' v2",
                "EFO:0009922": "10x 3' v3",
                "EFO:0011025": "10x 5' v1",
                "EFO:0008930": "Smart-seq",
                "EFO:0008931": "Smart-seq2"
            },
            "MONDO": {
                "MONDO:0100096": "COVID-19"
            },
            "MmusDv": {
                "MmusDv:0000062": "2 month-old stage",
                "MmusDv:0000003": "Theiler stage 01"
            },
            "HsapDv": {
                "HsapDv:0000174": "1 month-old human stage"
            },
            "NCBITaxon": {
                "NCBITaxon:9606": "Homo sapiens",
                "NCBITaxon:10090": "Mus musculus",
            },

            "UBERON": {
                "UBERON:0002048": "lung"
            },

            "PATO": {
                "PATO:0000461": "normal"
            }
        }

        self.invalid_terms = {
            "CL": ["EFO:0009899", "NO_TERM"],
            "EFO": ["UBERON:0002048", "NO_TERM"],
            "MONDO": ["EFO:0009899", "NO_TERM"],
            "MmusDv": ["EFO:0009899", "NO_TERM"],
            "NCBITaxon": ["EFO:0009899", "NO_TERM"],
            "UBERON": ["EFO:0009899", "NO_TERM"],
            "PATO": ["EFO:0009899", "NO_TERM"]
        }

    def test_ontology_validity(self):
        for i in self.valid_ontologies:
            self.assertTrue(self.ontologyChecker.is_valid_ontology(i))
            self.assertIsNone(self.ontologyChecker.assert_ontology(i))

        for i in self.invalid_ontologies:
            self.assertFalse(self.ontologyChecker.is_valid_ontology(i))
            with self.assertRaises(ValueError):
                self.ontologyChecker.assert_ontology(i)

    def test_valid_term_id(self):
        for ontology in self.valid_terms:
            for term_id in self.valid_terms[ontology]:
                term_label = self.valid_terms[ontology][term_id]

                self.assertTrue(self.ontologyChecker.is_valid_term_id(ontology, term_id))
                self.assertIsNone(self.ontologyChecker.assert_term_id(ontology, term_id))
                self.assertEqual(self.ontologyChecker.get_term_label(ontology, term_id), term_label)

    def test_invalid_term_ids(self):
        for ontology in self.invalid_terms:
            for term_id in self.invalid_terms[ontology]:
                self.assertFalse(self.ontologyChecker.is_valid_term_id(ontology, term_id))
                with self.assertRaises(ValueError):
                    self.ontologyChecker.assert_term_id(ontology, term_id)