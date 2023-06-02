import os
import subprocess
from typing import List

from jinja2 import Template

file_path = os.path.dirname(os.path.realpath(__file__))
convert_file_path = os.path.join(file_path, "../cellxgene_schema_cli/cellxgene_schema/convert.py")


def get_current_version() -> str:
    return subprocess.run(["make", "show-current-version"], capture_output=True, text=True).stdout.strip()


def get_template() -> str:
    with open(os.path.join(file_path, "convertion_template.jinja"), "r") as fp:
        template = fp.read()
    j2_template = Template(template, trim_blocks=True, lstrip_blocks=True)
    return j2_template


def generate_script(template: Template, schema_version: str, ontology_term_map: dict, gencode_term_map: dict):
    output = template.render(schema_version=schema_version, ontology_term_map=ontology_term_map, gencode_term_map=None)

    # Overwrite the existing convert.py file
    with open(os.path.join(file_path, convert_file_path), "w") as fp:
        fp.write(output)


def get_ontology_term_map() -> dict:
    # TODO: read in the ontology term map
    return {
        "assay": {},
        "cell_type": {},
        "development_stage": {},
        "disease": {},
        "organism": {},
        "self_reported_ethnicity": {},
        "sex": {},
        "tissue": {},
    }


def get_genecode_term_map() -> List[str]:
    # TODO: read in the gencode term map
    return []


def main():
    template = get_template()
    schema_version = get_current_version()
    ontology_term_map = get_ontology_term_map()
    gencode_term_map = get_genecode_term_map()
    generate_script(template, schema_version, ontology_term_map, gencode_term_map)


if __name__ == "__main__":
    main()
