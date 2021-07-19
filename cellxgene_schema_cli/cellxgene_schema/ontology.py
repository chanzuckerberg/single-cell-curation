"""Methods for working with ontologies and the OLS."""
import os
import gzip
import json


class geneChecker:
    """Handles checking gene ids, retrieves symbols"""

    GENE_FILES = {
        "Homo sapiens": os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "ontology_files/genes_homo_sapiens.csv.gz",
        ),
        "Mus musculus": os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "ontology_files/genes_mus_musculus.csv.gz",
        ),
    }

    def __init__(self, species):
        """
        :param str species: scientific name of a species, currently the only supported species are: "Homo sapiens" and "Mus musculus"
        """

        if species not in self.GENE_FILES:
            raise ValueError(f"{species} not supported.")

        self.species = species
        self.gene_dict = {}
        with gzip.open(self.GENE_FILES[species], "rt") as genes:
            for gene in genes:
                gene = gene.rstrip().split(",")
                self.gene_dict[gene[0]] = gene[1]

    def is_valid_id(self, gene_id):
        """
        Returns True if the gene_id is a valid ENSEMBL id
        """

        return gene_id in self.gene_dict

    def assert_gene_id(self, gene_id):
        """
        Raises error if the gene_id is not a valid ENSEBML id
        """

        if not self.is_valid_id(gene_id):
            raise ValueError(
                f"The id '{gene_id}' is not a valid ENSEMBL id for '{self.species}'"
            )

    def get_symbol(self, gene_id):
        """
        Returns the symbol associated to the ENSEBML id
        """

        self.assert_gene_id(gene_id)

        return self.gene_dict[gene_id]


class ontologyChecker:
    """Handles checking ontology term ids, retrieves ontology labels and ancestors"""

    JSON_FILE = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "ontology_files/all_ontology.json.gz",
    )

    def __init__(self):
        with gzip.open(self.JSON_FILE, "rt") as json_o:
            self.ontology_dict = json.load(json_o)

    def get_ontologies(self):
        """
        Returns a list of ontologies available in the checker
        """

        return list(self.ontology_dict.keys())

    def get_info(self, term_id, mode="all"):

        """
        Returns information from an ontology term id based on "mode"

        :param str term_id: the ontology term id
        :param str mode: the kind of information to retrieve: all|label|ancestors
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

    def print_info(self, term_id, mode="all"):
        """
        Prints information from an ontology term id based on "mode"

        :param str term_id: path to gzipped gtf file
        :param str mode: the kind of information to retrieve: all|label|ancestors
        """

        if mode == "all":
            print(json.dumps(self.get_info(term_id, mode), indent=4))
        else:
            print(self.get_info(term_id, mode))

    def get_term_dict(self, ontology, term_id):
        """
        Returns a dictionary with all the information from a given ontology and term_id
        """

        self.assert_term_id(ontology, term_id)
        return self.ontology_dict[ontology][term_id]

    def get_term_label(self, ontology, term_id):
        """
        Returns the label associated to an ontology term id
        """

        self.assert_term_id(ontology, term_id)
        return self.ontology_dict[ontology][term_id]["label"]

    def get_term_ancestors(self, ontology, term_id):
        """
        Returns the ancestros of an ontology id
        """

        self.assert_term_id(ontology, term_id)
        return set(self.ontology_dict[ontology][term_id]["ancestors"])

    def is_valid_ontology(self, ontology):
        """
        Returns True if the ontology is present in the ontology dict
        """

        return ontology in self.ontology_dict

    def is_valid_term_id(self, ontology, term_id):
        """
        Returns True id the id is valid from ontology
        """

        self.assert_ontology(ontology)

        return term_id in self.ontology_dict[ontology]

    def is_descendent_of(self, ontology, query_term_id, target_term_id):
        """
        Returns true if query_term_id is a descendent of target_term_id in a given ontology
        """

        self.assert_term_id(ontology, query_term_id)
        self.assert_term_id(ontology, target_term_id)

        return target_term_id in self.get_term_ancestors(ontology, query_term_id)

    def assert_ontology(self, ontology):
        """
        Raises error if ontology is not present in the ontology dict
        """

        if not self.is_valid_ontology(ontology):
            raise ValueError(
                f"The ontology '{ontology}' is not present in the ontology checker"
            )

    def assert_term_id(self, ontology, term_id):
        """
        Raises error if ontology is not present in the ontology dict
        """

        if not self.is_valid_term_id(ontology, term_id):
            raise ValueError(
                f"The term id '{term_id}' is not present in the ontology '{ontology}'"
            )

    def assert_descendent_of(self, ontology, query_term_id, target_term_id):
        """
        Raises error if query_term_id is not a descendent of target_term_id in a given ontology
        """

        if not self.is_descendent_of(ontology, query_term_id, target_term_id):
            raise ValueError(
                f"The term id '{query_term_id}' is not a descendent of the term id '{target_term_id}'"
                f" in the ontology '{ontology}'"
            )
