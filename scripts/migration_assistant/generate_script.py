import json
import os
from typing import List

from jinja2 import Template

from cellxgene_schema_cli.cellxgene_schema import env

file_path = os.path.dirname(os.path.realpath(__file__))
target_file = os.path.join(file_path, "../../cellxgene_schema_cli/cellxgene_schema/migrate.py")


def get_template() -> Template:
    with open(os.path.join(file_path, "migration_template.jinja"), "r") as fp:
        input_file = fp.read()
    template = Template(input_file, keep_trailing_newline=True, trim_blocks=True, lstrip_blocks=True)
    return template


def generate_script(template: Template, ontology_term_map: dict, deprecated_feature_ids: List[str]):  # type: ignore
    output = template.render(ontology_term_map=ontology_term_map, deprecated_feature_ids=deprecated_feature_ids)

    # Overwrite the existing migrate.py file
    with open(target_file, "w") as fp:
        fp.write(output)


def get_ontology_term_map(term_map_filepath) -> dict:  # type: ignore
    if os.path.exists(term_map_filepath):
        with open(term_map_filepath, "r") as fp:
            replaced_by_map = json.load(fp)
        return replaced_by_map  # type: ignore


def get_deprecated_feature_ids() -> List[str]:
    # return a list of deprecated feature ids.
    diff_list = []
    suffix = "_diff.txt"
    files = os.listdir(env.GENCODE_DIR)
    for file in files:
        if file.endswith(suffix):
            with open(f"{env.GENCODE_DIR}/{file}") as fp:
                lines = fp.read().splitlines()
                diff_list.extend(lines)
    return diff_list


def main():  # type: ignore
    template = get_template()
    ontology_term_map = get_ontology_term_map("replaced-by.json")
    deprecated_feature_ids = get_deprecated_feature_ids()
    generate_script(template, ontology_term_map, deprecated_feature_ids)


if __name__ == "__main__":
    main()  # type: ignore
