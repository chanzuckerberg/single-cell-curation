import os

PACKAGE_ROOT = os.path.dirname(os.path.realpath(__file__))
GENCODE_DIR = os.path.join(PACKAGE_ROOT, "gencode_files")
GENE_INFO_YAML = os.path.join(GENCODE_DIR, "gene_info.yml")
SCHEMA_DEFINITIONS_DIR = os.path.join(PACKAGE_ROOT, "schema_definitions")
SCHEMA_DEFINITION_FILE = os.path.join(SCHEMA_DEFINITIONS_DIR, "schema_definition.yaml")
SCHEMA_REFERENCE_BASE_URL = "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema"
SCHEMA_REFERENCE_FILE_NAME = "schema.md"
