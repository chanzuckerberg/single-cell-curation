#!/usr/bin/env python
# coding: utf-8

# # Ontology ancestors mapping
# This notebook can be used to compute mappings between nodes in an ontology files and their ancestors.
# It will:
# 1. Load the required ontology files (human, mouse, UBERON)
# 1. For `human` and `mouse`, start from each node available in the ontology
# 1. for `UBERON`, only start from a defined, hardcoded subset of leaf nodes which is relevant for our datasets
# 1. Recursively compute ancestors using the `is_a` relationship and `part_of` restriction (defined as `BFO_0000050`)
# 1. Create mappings from each node to its ancestors
# 1. Write a json file with these mappings (a dictionary of lists)

# In[ ]:
import json
from collections import defaultdict
from typing import Dict, Iterable, Iterator, List

import yaml
from owlready2 import get_ontology  # noqa
from owlready2.entity import ThingClass
from owlready2.namespace import Ontology

# In[ ]:

RO__PART_OF = "BFO_0000050"  # RO ontology property "part of"


def get_ancestors(onto: Ontology, class_name: str) -> Iterator[str]:
    print("djh call")
    entity: ThingClass = onto.search_one(iri=f"http://purl.obolibrary.org/obo/{class_name}")

    # print(f"djh: {entity} {type(entity)}")
    # return entity
    def recurse(x: ThingClass) -> Iterator[str]:
        for e in x.is_a:
            if hasattr(e, "value") and e.property.name == RO__PART_OF:
                z = onto.search_one(iri=f"http://purl.obolibrary.org/obo/{e.value.name.replace('obo.', '')}")
                yield z.name
                yield from recurse(z)

    yield entity.name
    yield from recurse(entity)


# In[ ]:


def create_ancestors_mapping(onto: Ontology, classes: Iterable[ThingClass], prefix: str = None) -> Dict[str, List[str]]:
    ancestors: Dict[str, List[str]] = defaultdict(list)
    for cls in classes:
        class_key = cls.name.replace("_", ":")
        for name in get_ancestors(onto, cls.name):  # type: ignore
            if prefix and not name.startswith(prefix):
                continue
            ancestors[class_key].append(name.replace("_", ":"))
    return ancestors


# In[ ]:


if __name__ == "__main__":
    # Load owl.info to grab latest ontology sources
    owl_info_yml = "cellxgene_schema_cli/cellxgene_schema/ontology_files/owl_info.yml"
    with open(owl_info_yml, "r") as owl_info_handle:
        owl_info = yaml.safe_load(owl_info_handle)

    # In[ ]:

    human_latest_key = owl_info["HsapDv"]["latest"]
    human_ontology = owl_info["HsapDv"]["urls"][human_latest_key]

    human_onto = get_ontology(human_ontology)
    human_onto.load()

    m_human = create_ancestors_mapping(human_onto, human_onto.classes(), "HsapDv")  # type: ignore

    # In[ ]:

    mouse_latest_key = owl_info["MmusDv"]["latest"]
    mouse_ontology = owl_info["MmusDv"]["urls"][mouse_latest_key]

    mouse_onto = get_ontology(mouse_ontology)
    mouse_onto.load()

    m_mouse = create_ancestors_mapping(mouse_onto, mouse_onto.classes(), "MmusDv")  # type: ignore

    # In[ ]:

    uberon_latest_key = owl_info["UBERON"]["latest"]
    uberon_ontology = owl_info["UBERON"]["urls"][uberon_latest_key]

    tissue_onto = get_ontology(uberon_ontology)
    tissue_onto.load()

    # In[ ]:

    uberon_classes = [
        "UBERON_0007236",
        "UBERON_0000106",
        "UBERON_0014859",
        "UBERON_0008264",
        "UBERON_0007233",
        "UBERON_0000112",
        "UBERON_8000003",
        "UBERON_0014857",
        "UBERON_0009849",
        "UBERON_0034920",
        "UBERON_0000069",
        "UBERON_0000109",
        "UBERON_8000001",
        "UBERON_0000068",
        "UBERON_0018685",
        "UBERON_0000107",
        "UBERON_0007222",
        "UBERON_0000092",
        "UBERON_0018378",
        "UBERON_0014864",
        "UBERON_0004730",
        "UBERON_0000111",
        "UBERON_0007220",
        "UBERON_0014405",
        "UBERON_0014862",
        "UBERON_8000000",
        "UBERON_0000071",
        "UBERON_0014860",
        "UBERON_0012101",
        "UBERON_0000113",
        "UBERON_0014858",
        "UBERON_0007232",
        "UBERON_0000070",
        "UBERON_0000110",
        "UBERON_8000002",
        "UBERON_0014856",
        "UBERON_0004728",
        "UBERON_0034919",
        "UBERON_0000108",
        "UBERON_0000066",
        "UBERON_0004707",
        "UBERON_0000105",
        "UBERON_0018241",
        "UBERON_0007221",
        "UBERON_0014406",
        "UBERON_0014863",
        "UBERON_0004729",
        "UBERON_0014861",
    ]

    # In[ ]:

    uc = [tissue_onto.search_one(iri=f"http://purl.obolibrary.org/obo/{c}") for c in uberon_classes]

    m_uberon = create_ancestors_mapping(tissue_onto, uc, "UBERON")  # type: ignore

    # In[ ]:

    with open("development_stage_ontology_mapping.json", "w") as f:
        d = {}
        d.update(m_human)
        d.update(m_mouse)
        d.update(m_uberon)
        json.dump(d, f)
