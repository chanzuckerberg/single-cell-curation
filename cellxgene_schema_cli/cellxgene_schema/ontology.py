import os
import gzip
import json
import enum
from typing import Union, List, Set
from . import env


class SupportedOrganisms(enum.Enum):
    Homo_sapiens = 1
    Mus_musculus = 2


class GeneChecker:
    """Handles checking gene ids, retrieves symbols"""

    GENE_FILES = {
        SupportedOrganisms.Homo_sapiens: os.path.join(
            env.ONTOLOGY_DIR, "genes_homo_sapiens.csv.gz"
        ),
        SupportedOrganisms.Mus_musculus: os.path.join(
            env.ONTOLOGY_DIR, "genes_mus_musculus.csv.gz"
        ),
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
                gene = gene.rstrip().split(",")
                self.gene_dict[gene[0]] = gene[1]

    def is_valid_id(self, gene_id: str) -> bool:
        """
        Checks for validity of gene id

        :param str gene_id: ENSEMBL gene id

        :rtype bool
        :return True if the gene_id is a valid ENSEMBL id, False otherwise
        """

        return gene_id in self.gene_dict

    def assert_gene_id(self, gene_id: str):
        """
        Raises error if the gene_id is not a valid ENSEBML id

        :param str gene_id: ENSEMBL gene id

        :rtype None
        """

        if not self.is_valid_id(gene_id):
            raise ValueError(
                f"The id '{gene_id}' is not a valid ENSEMBL id for '{self.species}'"
            )

    def get_symbol(self, gene_id) -> str:
        """
        Gets symbol associated to the ENSEBML id

        :param str gene_id: ENSEMBL gene id

        :rtype str
        :return A gene symbol
        """

        self.assert_gene_id(gene_id)

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

    def get_info(self, term_id: str, mode: str = "all") -> Union[dict, str, Set[str]]:

        """
        Returns information from an ontology term id based on "mode"

        :param str term_id: the ontology term id
        :param str mode: the kind of information to retrieve: all|label|ancestors

        rtype Union[dict, str, Set[str]]:
        return: If mode is "all" retrieve dictionary with all info for the term id, if mode is "label" returns the label
        associated to term id, if mode is "ancestors" returns the ids for the ancestors associated to term_id
        """

        onto_and_term = term_id.split(":")
        ontology = onto_and_term[0]

        if len(onto_and_term) != 2:
            raise ValueError(
                f"'{term_id}' term id's format is invalid, it must be of the from 'Ontology:Id'"
            )

        if mode == "all":
            return self.get_term_dict(ontology, term_id)
        elif mode == "label":
            return self.get_term_label(ontology, term_id)
        elif mode == "ancestors":
            return self.get_term_ancestors(ontology, term_id)
        else:
            raise ValueError(f"{mode} is not a supported type")

    def print_info(self, term_id: str, mode: str = "all"):
        """
        Prints information from an ontology term id based on "mode". If mode is "all" print all info for the term id,
        if mode is "label" prints the label associated to term id, if mode is "ancestors" print the ids for the
        ancestors associated to term_id.

        :param str term_id: the ontology term id
        :param str mode: the kind of information to retrieve: all|label|ancestors

        :rtype None
        """

        if mode == "all":
            print(json.dumps(self.get_info(term_id, mode), indent=4))
        else:
            print(self.get_info(term_id, mode))

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

    def is_descendent_of(
        self, ontology: str, query_term_id: str, target_term_id: str
    ) -> bool:
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
            raise ValueError(
                f"The ontology '{ontology}' is not present in the ontology checker"
            )

    def assert_term_id(self, ontology: str, term_id: str):
        """
        Raises error if term_id is not present in ontology

        :param str ontology: the ontology id
        :param str term_id: the ontology term id

        :rtype None
        """

        if not self.is_valid_term_id(ontology, term_id):
            raise ValueError(
                f"The term id '{term_id}' is not present in the ontology '{ontology}'"
            )

    def assert_descendent_of(
        self, ontology: str, query_term_id: str, target_term_id: str
    ):
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
