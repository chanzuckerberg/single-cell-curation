import os

PACKAGE_ROOT = os.path.dirname(os.path.realpath(__file__))
ONTOLOGY_DIR = os.path.join(PACKAGE_ROOT, "ontology_files")
OWL_INFO_YAML = os.path.join(ONTOLOGY_DIR, "owl_info.yml")
PARSED_ONTOLOGIES_FILE = os.path.join(ONTOLOGY_DIR, "all_ontology.json.gz")
SCHEMA_DEFINITIONS_DIR = os.path.join(PACKAGE_ROOT, "schema_definitions")
