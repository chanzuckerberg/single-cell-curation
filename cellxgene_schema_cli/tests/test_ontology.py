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
    valid_genes_same_name_diff_species,
    valid_genes_same_name_and_species,
)

# Tests for internal functions of the OntologyChecker and GeneChecker classes


class TestGeneChecker:
    @pytest.mark.parametrize("valid_species", ontology.SupportedOrganisms)
    def test_species_valid(self, valid_species):
        assert isinstance(ontology.GeneChecker(valid_species), ontology.GeneChecker)

    @pytest.mark.parametrize("invalid_species", invalid_species)
    def test_species_invalid(self, invalid_species):
        with pytest.raises(ValueError):
            ontology.GeneChecker(invalid_species)

    @pytest.mark.parametrize("species,valid_genes", valid_genes.items())
    def test_valid_genes(self, species, valid_genes):
        geneChecker = ontology.GeneChecker(species)
        for gene_id in valid_genes:
            gene_label = valid_genes[gene_id][0]
            gene_length = valid_genes[gene_id][1]

            assert geneChecker.is_valid_id(gene_id)
            assert geneChecker.get_symbol(gene_id) == gene_label
            assert geneChecker.get_length(gene_id) == gene_length

    @pytest.mark.parametrize("species,invalid_genes", invalid_genes.items())
    def test_invalid_genes(self, species, invalid_genes):
        geneChecker = ontology.GeneChecker(species)
        for gene_id in invalid_genes:
            assert not geneChecker.is_valid_id(gene_id)
            with pytest.raises(ValueError):
                geneChecker.get_symbol(gene_id)
            with pytest.raises(ValueError):
                geneChecker.get_length(gene_id)
    
    @pytest.mark.parametrize("species,valid_genes_same_name_diff_species", valid_genes_same_name_diff_species.items())
    def test_valid_genes_same_name_diff_species(self, species, valid_genes_same_name_diff_species):
        geneChecker = ontology.GeneChecker(species)
        for gene_id in valid_genes_same_name_diff_species:
            gene_label = valid_genes_same_name_diff_species[gene_id][0]
            gene_length = valid_genes_same_name_diff_species[gene_id][1]

            assert geneChecker.is_valid_id(gene_id)
            assert geneChecker.get_symbol(gene_id) == gene_label
            assert geneChecker.get_length(gene_id) == gene_length

    @pytest.mark.parametrize("species,valid_genes_same_name_and_species", valid_genes_same_name_and_species.items())
    def test_valid_genes_same_name_and_species(self, species, valid_genes_same_name_and_species):
        geneChecker = ontology.GeneChecker(species)
        for gene_id in valid_genes_same_name_and_species:
            gene_label = valid_genes_same_name_and_species[gene_id][0]
            gene_length = valid_genes_same_name_and_species[gene_id][1]

            assert geneChecker.is_valid_id(gene_id)
            assert geneChecker.get_symbol(gene_id) == gene_label
            assert geneChecker.get_length(gene_id) == gene_length


@pytest.fixture(scope="class")
def ontologyChecker() -> ontology.OntologyChecker:
    return ontology.OntologyChecker()


class TestOntologyChecker:
    @pytest.mark.parametrize("valid_ontology", valid_ontologies)
    def test_ontology_valid(self, ontologyChecker, valid_ontology):
        assert ontologyChecker.is_valid_ontology(valid_ontology)
        assert ontologyChecker.assert_ontology(valid_ontology) is None

    @pytest.mark.parametrize("invalid_ontology", invalid_ontologies)
    def test_ontology_invalid(self, ontologyChecker, invalid_ontology):
        assert not ontologyChecker.is_valid_ontology(invalid_ontology)
        with pytest.raises(ValueError):
            ontologyChecker.assert_ontology(invalid_ontology)

    @pytest.mark.parametrize(
        "ontology_id,term_id",
        [(ontology_id, term_id) for ontology_id in valid_terms for term_id in valid_terms[ontology_id]],
    )
    def test_valid_term_id(self, ontologyChecker, ontology_id, term_id):
        term_label = valid_terms[ontology_id][term_id]

        assert ontologyChecker.is_valid_term_id(ontology_id, term_id)
        assert ontologyChecker.assert_term_id(ontology_id, term_id) is None
        assert ontologyChecker.get_term_label(ontology_id, term_id) == term_label

    @pytest.mark.parametrize(
        "ontology_id,term_id",
        [(ontology_id, term_id) for ontology_id in invalid_terms for term_id in invalid_terms[ontology_id]],
    )
    def test_invalid_term_ids(self, ontologyChecker, ontology_id, term_id):
        assert not ontologyChecker.is_valid_term_id(ontology_id, term_id)
        with pytest.raises(ValueError):
            ontologyChecker.assert_term_id(ontology_id, term_id)
