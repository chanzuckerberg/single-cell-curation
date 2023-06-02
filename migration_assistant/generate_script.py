import os
import subprocess

from jinja2 import Template

schema_version = subprocess.run(["make", "show-current-version"], capture_output=True, text=True).stdout.strip()
ontology_term_map = {
    "assay": {},
    "cell_type": {},
    "development_stage": {},
    "disease": {},
    "organism": {},
    "self_reported_ethnicity": {},
    "sex": {},
    "tissue": {},
}


file_path = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(file_path, "convertion_template.jinja"), "r") as fp:
    template = fp.read()

j2_template = Template(template, trim_blocks=True, lstrip_blocks=True)
report = j2_template.render(schema_version=schema_version, ontology_term_map=ontology_term_map, gencode_term_map=None)

# Overwrite the existing convert.py file
convert_file_path = os.path.dirname("../cellxgene_schema_cli/cellxgene_schema/convert.py")
with open(os.path.join(file_path, "../cellxgene_schema_cli/cellxgene_schema/convert.py"), "w") as fp:
    fp.write(report)
