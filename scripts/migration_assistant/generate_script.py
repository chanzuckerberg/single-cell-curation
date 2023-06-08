import os

from jinja2 import Template

from cellxgene_schema_cli.cellxgene_schema.ontology import get_deprecated_genecode_terms

file_path = os.path.dirname(os.path.realpath(__file__))
target_file = os.path.join(file_path, "../../cellxgene_schema_cli/cellxgene_schema/migrate.py")


def get_template() -> str:
    with open(os.path.join(file_path, "migration_template.jinja"), "r") as fp:
        input_file = fp.read()
    template = Template(input_file, trim_blocks=True, lstrip_blocks=True)
    return template


def generate_script(template: Template, ontology_term_map: dict, include_gencode_conversion: bool):
    output = template.render(ontology_term_map=ontology_term_map, include_gencode_conversion=include_gencode_conversion)

    # Overwrite the existing migrate.py file
    with open(target_file, "w") as fp:
        fp.write(output)


def get_ontology_term_map() -> bool:
    # TODO: read in the ontology term map
    return len(get_deprecated_genecode_terms()) > 1


def main():
    template = get_template()
    ontology_term_map = get_ontology_term_map()
    include_gencode_conversion = len(get_deprecated_genecode_terms()) > 1
    generate_script(template, ontology_term_map, include_gencode_conversion)


if __name__ == "__main__":
    main()
