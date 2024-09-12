import enum
import gzip
import os
from typing import Union

from . import env


class SupportedOrganisms(enum.Enum):
    HOMO_SAPIENS = "NCBITaxon:9606"
    MUS_MUSCULUS = "NCBITaxon:10090"
    SARS_COV_2 = "NCBITaxon:2697049"
    ERCC = "NCBITaxon:32630"


def get_organism_from_feature_id(
    feature_id: str,
) -> Union[SupportedOrganisms, None]:
    """
    Infers the organism of a feature id based on the prefix of a feature id, e.g. ENSG means Homo sapiens

    :param str feature_id: the feature id

    :rtype Union[ontology.SypportedOrganisms, None]
    :return: the organism the feature id is from
    """

    if feature_id.startswith("ENSG") or feature_id.startswith("ENST"):
        return SupportedOrganisms.HOMO_SAPIENS
    elif feature_id.startswith("ENSMUS"):
        return SupportedOrganisms.MUS_MUSCULUS
    elif feature_id.startswith("ENSSAS"):
        return SupportedOrganisms.SARS_COV_2
    elif feature_id.startswith("ERCC-"):
        return SupportedOrganisms.ERCC
    else:
        return None


class GeneChecker:
    """Handles checking gene ids, retrieves symbols"""

    GENE_FILES = {
        SupportedOrganisms.HOMO_SAPIENS: os.path.join(env.GENCODE_DIR, "genes_homo_sapiens.csv.gz"),
        SupportedOrganisms.MUS_MUSCULUS: os.path.join(env.GENCODE_DIR, "genes_mus_musculus.csv.gz"),
        SupportedOrganisms.SARS_COV_2: os.path.join(env.GENCODE_DIR, "genes_sars_cov_2.csv.gz"),
        SupportedOrganisms.ERCC: os.path.join(env.GENCODE_DIR, "genes_ercc.csv.gz"),
    }

    def __init__(self, species: SupportedOrganisms):
        """
        :param enum.Enum.SupportedSpecies species: item from SupportedOrganisms
        """

        if species not in self.GENE_FILES:
            raise ValueError(f"{species} not supported.")

        self.species = species
        self.gene_dict = {}
        with gzip.open(self.GENE_FILES[species], "rt") as genes:
            for gene in genes:
                gene = gene.rstrip().split(",")  # type: ignore
                gene_id = gene[0]
                gene_label = gene[1]
                gene_length = int(gene[3])
                gene_type = gene[4]

                self.gene_dict[gene_id] = (gene_label, gene_length, gene_type)

    def is_valid_id(self, gene_id: str) -> bool:
        """
        Checks for validity of gene id

        :param str gene_id: ENSEMBL gene id

        :rtype bool
        :return True if the gene_id is a valid ENSEMBL id, False otherwise
        """

        return gene_id in self.gene_dict

    def get_symbol(self, gene_id: str) -> str:
        """
        Gets symbol associated to the ENSEBML id

        :param str gene_id: ENSEMBL gene id

        :rtype str
        :return A gene symbol
        """

        if self.is_valid_id(gene_id):
            return self.gene_dict[gene_id][0]
        else:
            raise ValueError(f"The id '{gene_id}' is not a valid ENSEMBL id for '{self.species}'")

    def get_length(self, gene_id: str) -> int:
        """
        Gets feature length associated to the ENSEBML id

        :param str gene_id: ENSEMBL gene id

        :rtype int
        :return A gene length
        """

        if self.is_valid_id(gene_id):
            return self.gene_dict[gene_id][1]
        else:
            raise ValueError(f"The id '{gene_id}' is not a valid ENSEMBL id for '{self.species}'")

    def get_type(self, gene_id: str) -> str:
        """
        Gets feature type associated to the ENSEBML id

        :param str gene_id: ENSEMBL gene id

        :rtype str
        :return A feature type
        """

        if self.is_valid_id(gene_id):
            return self.gene_dict[gene_id][2]
        else:
            raise ValueError(f"The id '{gene_id}' is not a valid ENSEMBL id for '{self.species}'")
