import json
import os

from jinja2 import Template

from cellxgene_schema_cli.cellxgene_schema.utils import get_deprecated_features

file_path = os.path.dirname(os.path.realpath(__file__))
target_file = os.path.join(file_path, "../../cellxgene_schema_cli/cellxgene_schema/migrate.py")


def get_template() -> Template:
    with open(os.path.join(file_path, "migration_template.jinja"), "r") as fp:
        input_file = fp.read()
    template = Template(input_file, trim_blocks=True, lstrip_blocks=True)
    return template


def generate_script(template: Template, ontology_term_map: dict, gencode_migrate: bool):
    output = template.render(ontology_term_map=ontology_term_map, gencode_migrate=gencode_migrate)

    # Overwrite the existing migrate.py file
    with open(target_file, "w") as fp:
        fp.write(output)


def get_ontology_term_map(term_map_filepath) -> dict:
    if os.path.exists(term_map_filepath):
        with open(term_map_filepath, "r") as fp:
            replaced_by_map = json.load(fp)
        return replaced_by_map


def migrate_gencode() -> bool:
    return len(get_deprecated_features()) > 1


def main():
    template = get_template()
    ontology_term_map = get_ontology_term_map("replaced-by.json")
    gencode_migrate = migrate_gencode()
    generate_script(template, ontology_term_map, gencode_migrate)


if __name__ == "__main__":
    main()
