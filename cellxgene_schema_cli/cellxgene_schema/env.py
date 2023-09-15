import os

PACKAGE_ROOT = os.path.dirname(os.path.realpath(__file__))
ONTOLOGY_DIR = os.path.join(PACKAGE_ROOT, "ontology_files")
GENE_INFO_YAML = os.path.join(ONTOLOGY_DIR, "gene_info.yml")
OWL_INFO_YAML = os.path.join(ONTOLOGY_DIR, "owl_info.yml")
PARSED_ONTOLOGIES_FILE = os.path.join(ONTOLOGY_DIR, "all_ontology.json.gz")
SCHEMA_DEFINITIONS_DIR = os.path.join(PACKAGE_ROOT, "schema_definitions")
SCHEMA_DEFINITION_FILE = os.path.join(SCHEMA_DEFINITIONS_DIR, "schema_definition.yaml")
SCHEMA_REFERENCE_BASE_URL = "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema"
SCHEMA_REFERENCE_FILE_NAME = "schema.md"
