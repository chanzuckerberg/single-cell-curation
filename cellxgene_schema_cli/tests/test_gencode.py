import pytest
from cellxgene_schema import gencode
from fixtures.examples_ontology_test import (
    invalid_genes,
    invalid_species,
    valid_genes,
    valid_genes_same_name_and_species,
    valid_genes_same_name_diff_species,
)

# Tests for internal functions of the OntologyChecker and GeneChecker classes


class TestGeneChecker:
    @pytest.mark.parametrize("valid_species", gencode.SupportedOrganisms)
    def test_species_valid(self, valid_species):
        assert isinstance(gencode.GeneChecker(valid_species), gencode.GeneChecker)

    @pytest.mark.parametrize("invalid_species", invalid_species)
    def test_species_invalid(self, invalid_species):
        with pytest.raises(ValueError):
            gencode.GeneChecker(invalid_species)

    @pytest.mark.parametrize("species,valid_genes", valid_genes.items())
    def test_valid_genes(self, species, valid_genes):
        geneChecker = gencode.GeneChecker(species)
        for gene_id in valid_genes:
            gene_label = valid_genes[gene_id][0]
            gene_length = valid_genes[gene_id][1]

            assert geneChecker.is_valid_id(gene_id)
            assert geneChecker.get_symbol(gene_id) == gene_label
            assert geneChecker.get_length(gene_id) == gene_length

    @pytest.mark.parametrize("species,invalid_genes", invalid_genes.items())
    def test_invalid_genes(self, species, invalid_genes):
        geneChecker = gencode.GeneChecker(species)
        for gene_id in invalid_genes:
            assert not geneChecker.is_valid_id(gene_id)
            with pytest.raises(ValueError):
                geneChecker.get_symbol(gene_id)
            with pytest.raises(ValueError):
                geneChecker.get_length(gene_id)

    @pytest.mark.parametrize("species,valid_genes_same_name_diff_species", valid_genes_same_name_diff_species.items())
    def test_valid_genes_same_name_diff_species(self, species, valid_genes_same_name_diff_species):
        geneChecker = gencode.GeneChecker(species)
        for gene_id in valid_genes_same_name_diff_species:
            gene_label = valid_genes_same_name_diff_species[gene_id][0]
            gene_length = valid_genes_same_name_diff_species[gene_id][1]

            assert geneChecker.is_valid_id(gene_id)
            assert geneChecker.get_symbol(gene_id) == gene_label
            assert geneChecker.get_length(gene_id) == gene_length

    @pytest.mark.parametrize("species,valid_genes_same_name_and_species", valid_genes_same_name_and_species.items())
    def test_valid_genes_same_name_and_species(self, species, valid_genes_same_name_and_species):
        geneChecker = gencode.GeneChecker(species)
        for gene_id in valid_genes_same_name_and_species:
            gene_label = valid_genes_same_name_and_species[gene_id][0]
            gene_length = valid_genes_same_name_and_species[gene_id][1]

            assert geneChecker.is_valid_id(gene_id)
            assert geneChecker.get_symbol(gene_id) == gene_label
            assert geneChecker.get_length(gene_id) == gene_length
