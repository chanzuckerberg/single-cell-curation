import enum
import gzip
import json
import os
from typing import List, Set, Union

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
        SupportedOrganisms.HOMO_SAPIENS: os.path.join(env.ONTOLOGY_DIR, "genes_homo_sapiens.csv.gz"),
        SupportedOrganisms.MUS_MUSCULUS: os.path.join(env.ONTOLOGY_DIR, "genes_mus_musculus.csv.gz"),
        SupportedOrganisms.SARS_COV_2: os.path.join(env.ONTOLOGY_DIR, "genes_sars_cov_2.csv.gz"),
        SupportedOrganisms.ERCC: os.path.join(env.ONTOLOGY_DIR, "genes_ercc.csv.gz"),
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
            gene_labels = set()
            duplicated_gene_labels = set()

            for gene in genes:
                gene = gene.rstrip().split(",")
                gene_id = gene[0]
                gene_label = gene[1]

                self.gene_dict[gene_id] = gene_label

                # Keeps track of duplicated gene labels
                if gene_label in gene_labels:
                    duplicated_gene_labels.add(gene_label)
                else:
                    gene_labels.add(gene_label)

            # Makes gene labels unique
            for gene_id, gene_label in self.gene_dict.items():
                if gene_label in duplicated_gene_labels:
                    self.gene_dict[gene_id] = gene_label + "_" + gene_id

    def is_valid_id(self, gene_id: str) -> bool:
        """
        Checks for validity of gene id

        :param str gene_id: ENSEMBL gene id

        :rtype bool
        :return True if the gene_id is a valid ENSEMBL id, False otherwise
        """

        return gene_id in self.gene_dict

    def get_symbol(self, gene_id) -> str:
        """
        Gets symbol associated to the ENSEBML id

        :param str gene_id: ENSEMBL gene id

        :rtype str
        :return A gene symbol
        """

        if not self.is_valid_id(gene_id):
            raise ValueError(f"The id '{gene_id}' is not a valid ENSEMBL id for '{self.species}'")

        return self.gene_dict[gene_id]


class OntologyChecker:
    """Handles checking ontology term ids, retrieves ontology labels and ancestors"""

    JSON_FILE = env.PARSED_ONTOLOGIES_FILE

    def __init__(self):
        with gzip.open(self.JSON_FILE, "rt") as json_o:
            self.ontology_dict = json.load(json_o)

    def get_ontologies(self) -> List[str]:
        """
        rtype list[str]
        return: a list of ontologies available in the checker
        """

        return list(self.ontology_dict.keys())

    def get_term_dict(self, ontology: str, term_id: str) -> dict:
        """
        Returns a dictionary with all the information from a given ontology and term_id

        :param str ontology: the ontology id
        :param str term_id: the ontology term id

        :rtype dict
        :return Dictionary with all the information for the given term id
        """

        self.assert_term_id(ontology, term_id)
        return self.ontology_dict[ontology][term_id]

    def get_term_label(self, ontology: str, term_id: str) -> str:
        """
        Returns the label associated to an ontology term id

        :param str ontology: the ontology id
        :param str term_id: the ontology term id

        :rtype str
        :return Label associated to term id
        """

        self.assert_term_id(ontology, term_id)
        return self.ontology_dict[ontology][term_id]["label"]

    def get_term_ancestors(self, ontology: str, term_id: str) -> Set[str]:
        """
        Returns the ancestors of an ontology id

        :param str ontology: the ontology id
        :param str term_id: the ontology term id

        :rtype Set[str]
        :return All term ids that are ancestors of the query term id.
        """

        self.assert_term_id(ontology, term_id)
        return set(self.ontology_dict[ontology][term_id]["ancestors"])

    def is_valid_ontology(self, ontology: str) -> bool:
        """
        Returns True if the ontology is present in the ontology dict

        :param str ontology: the ontology id

        :rtype bool
        :return True if the ontology is present in the ontology dict, False otherwise
        """

        return ontology in self.ontology_dict

    def is_term_id_deprecated(self, ontology: str, term_id: str) -> bool:
        """
        Returns True if the id has been deprecated (obsolete) in the ontology

        :param str ontology: the ontology id
        :param str term_id: the ontology term id

        :rtype bool
        :return True if id has been deprecated
        """

        self.assert_term_id(ontology, term_id)

        return self.ontology_dict[ontology][term_id]["deprecated"]

    def is_valid_term_id(self, ontology: str, term_id: str) -> bool:
        """
        Returns True if term_id is a valid id  from ontology

        :param str ontology: the ontology id
        :param str term_id: the ontology term id

        :rtype bool
        :return True if term_id is a valid id  from ontology, False otherwise
        """

        self.assert_ontology(ontology)

        return term_id in self.ontology_dict[ontology]

    def is_descendent_of(self, ontology: str, query_term_id: str, target_term_id: str) -> bool:
        """
        Returns True if query_term_id is a descendent of target_term_id in a given ontology

        :param str ontology: the ontology id
        :param str query_term_id: the ontology term id
        :param str target_term_id: the ontology term id

        :rtype bool
        :return True if query_term_id is a descendent of target_term_id in a given ontology, False otherwise
        """

        self.assert_term_id(ontology, query_term_id)
        self.assert_term_id(ontology, target_term_id)

        return target_term_id in self.get_term_ancestors(ontology, query_term_id)

    def assert_ontology(self, ontology: str):
        """
        Raises error if ontology is not present in the ontology dict

        :param str ontology: the ontology id

        :rtype None
        """

        if not self.is_valid_ontology(ontology):
            raise ValueError(f"The ontology '{ontology}' is not present in the ontology checker")

    def assert_term_id(self, ontology: str, term_id: str):
        """
        Raises error if term_id is not present in ontology

        :param str ontology: the ontology id
        :param str term_id: the ontology term id

        :rtype None
        """

        if not self.is_valid_term_id(ontology, term_id):
            raise ValueError(f"The term id '{term_id}' is not present in the ontology '{ontology}'")

    def assert_descendent_of(self, ontology: str, query_term_id: str, target_term_id: str):
        """
        Raises error if query_term_id is not a descendent of target_term_id in a given ontology

        :param str ontology: the ontology id
        :param str query_term_id: the ontology term id
        :param str target_term_id: the ontology term id

        :rtype None
        """

        if not self.is_descendent_of(ontology, query_term_id, target_term_id):
            raise ValueError(
                f"The term id '{query_term_id}' is not a descendent of the term id '{target_term_id}'"
                f" in the ontology '{ontology}'"
            )
