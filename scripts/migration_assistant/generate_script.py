import json
import os
from typing import List

from jinja2 import Template

file_path = os.path.dirname(os.path.realpath(__file__))
target_file = os.path.join(file_path, "../../cellxgene_schema_cli/cellxgene_schema/migrate.py")


def get_template() -> Template:
    with open(os.path.join(file_path, "migration_template.jinja"), "r") as fp:
        input_file = fp.read()
    template = Template(input_file, trim_blocks=True, lstrip_blocks=True)
    return template


def generate_script(template: Template, ontology_term_map: dict, gencode_term_map: list):
    output = template.render(ontology_term_map=ontology_term_map, gencode_term_map=None)

    # Overwrite the existing migrate.py file
    with open(target_file, "w") as fp:
        fp.write(output)


def get_ontology_term_map(term_map_filepath) -> dict:
    if os.path.exists(term_map_filepath):
        with open(term_map_filepath, "r") as fp:
            replaced_by_map = json.load(fp)
        return replaced_by_map


def get_genecode_term_map() -> List[str]:
    # TODO: read in the gencode term map
    return []


def main():
    template = get_template()
    ontology_term_map = get_ontology_term_map("replaced-by.json")
    gencode_term_map = get_genecode_term_map()
    generate_script(template, ontology_term_map, gencode_term_map)


if __name__ == "__main__":
    main()
