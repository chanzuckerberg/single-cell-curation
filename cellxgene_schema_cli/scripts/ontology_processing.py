import gzip
import json
import os
import re
import sys
import urllib.request
from threading import Thread
from urllib.error import HTTPError, URLError

import owlready2
import yaml

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../cellxgene_schema"))
import os
from typing import List

import env


def _download_owls(owl_info_yml: str = env.OWL_INFO_YAML, output_dir: str = env.ONTOLOGY_DIR):
    """
    Downloads the ontology owl files specified in 'owl_info_yml' into 'output_dir'

    :param str owl_info_yml: path to yaml file wit OWL information
    :param str output_dir: path to writable directory where owl files will be downloaded to

    :rtype None
    """

    with open(owl_info_yml, "r") as owl_info_handle:
        owl_info = yaml.safe_load(owl_info_handle)

    def download(_ontology, _url):
        print(f"Start Downloading {_ontology}")
        # Format of owl (handles cases where they are compressed)
        download_format = _url.split(".")[-1]

        output_file = os.path.join(output_dir, _ontology + ".owl")
        if download_format == "gz":
            urllib.request.urlretrieve(_url, output_file + ".gz")
            _decompress(output_file + ".gz", output_file)
            os.remove(output_file + ".gz")
        else:
            urllib.request.urlretrieve(_url, output_file)
        print(f"Finish Downloading {_ontology}")

    threads = []
    for ontology, _ in owl_info.items():
        latest_version = owl_info[ontology]["latest"]
        url = owl_info[ontology]["urls"][latest_version]
        try:
            urllib.request.urlopen(url)
        except HTTPError as e:
            raise Exception(f"{ontology} with pinned URL {url} returns status code {e.code}") from e
        except URLError as e:
            raise Exception(f"{ontology} with pinned URL {url} fails due to {e.reason}") from e

        t = Thread(target=download, args=(ontology, url))
        t.start()
        threads.append(t)

    for t in threads:
        t.join()


def _decompress(infile: str, tofile: str):
    """
    Decompresses a gziped file

    :param str infile: path gziped file
    :param str tofile: path to output decompressed file

    :rtype None
    """
    with open(infile, "rb") as inf, open(tofile, "w", encoding="utf8") as tof:
        decom_str = gzip.decompress(inf.read()).decode("utf-8")
        tof.write(decom_str)


def _parse_owls(
    working_dir: str = env.ONTOLOGY_DIR,
    owl_info_yml: str = env.OWL_INFO_YAML,
    output_json_file: str = env.PARSED_ONTOLOGIES_FILE,
):
    """
    Parser all owl files in working_dir. Extracts information from all classes in the owl file.
    The extracted information is written into a gzipped a json file with the following structure:
    {
        "ontology_name":
            {
            "term_id": {
                "label": "..."
                "deprecated": True
                "ancestors": [
                    "ancestor1_term_id_1",
                    "ancestor2_term_id_2"
                    ]
                }
            }

            "term_id2": {
                ...
            }

            ...
            }
    }

    :param str working_dir: path to folder with owl files
    :param str owl_info_yml: path to writable directory where owl files will be downloaded to
    :param str owl_info_yml: path to yaml file wit owl information
    :param str output_json_file: path to output jsaon file

    :rtype None
    """

    with open(owl_info_yml, "r") as owl_info_handle:
        owl_info = yaml.safe_load(owl_info_handle)

    owl_files = []
    for owl_file in os.listdir(working_dir):
        if owl_file.endswith(".owl"):
            owl_files.append(os.path.join(working_dir, owl_file))

    # Parse owl files
    onto_dict = {}
    for owl_file in owl_files:
        world = owlready2.World()
        onto = world.get_ontology(owl_file)
        onto.load()
        onto_dict[onto.name] = {}

        print(f"Processing {onto.name}")

        for onto_class in onto.classes():
            term_id = onto_class.name.replace("_", ":")

            # Skip terms that are not direct children from this ontology
            if onto.name != term_id.split(":")[0]:
                continue

            # If there are specified target terms then only work with them
            if onto.name in owl_info and "only" in owl_info[onto.name] and term_id not in owl_info[onto.name]["only"]:
                continue

            # Gets label
            onto_dict[onto.name][term_id] = dict()
            try:
                onto_dict[onto.name][term_id]["label"] = onto_class.label[0]
            except IndexError:
                onto_dict[onto.name][term_id]["label"] = ""

            # Add the "deprecated" status
            onto_dict[onto.name][term_id]["deprecated"] = False
            if onto_class.deprecated and onto_class.deprecated.first():
                # if deprecated, include information to determine replacement term(s)
                onto_dict[onto.name][term_id]["deprecated"] = True
                if onto_class.comment:
                    onto_dict[onto.name][term_id]["comments"] = [str(c) for c in onto_class.comment]
                # stores term tracking URL, such as a github issue discussing deprecated term
                if hasattr(onto_class, "IAO_0000233") and onto_class.IAO_0000233:
                    onto_dict[onto.name][term_id]["term_tracker"] = str(onto_class.IAO_0000233[0])

                # only need to record replaced_by OR considers
                if onto_class.IAO_0100001 and onto_class.IAO_0100001.first():
                    # url --> term
                    ontology_term = re.findall(r"[^\W_]+", str(onto_class.IAO_0100001[0]))
                    onto_dict[onto.name][term_id]["replaced_by"] = f"{ontology_term[-2]}:{ontology_term[-1]}"
                else:
                    if hasattr(onto_class, "consider") and onto_class.consider:
                        onto_dict[onto.name][term_id]["consider"] = [str(c) for c in onto_class.consider]
            # Gets ancestors
            ancestors = _get_ancestors(onto_class, onto.name)

            # If "children_of" specified in owl info then skip the current term if it is
            # not a children of those indicated.
            if (onto.name in owl_info and "children_of" in owl_info[onto.name]) and (
                not list(set(ancestors) & set(owl_info[onto.name]["children_of"]))
            ):
                onto_dict[onto.name].pop(term_id)
                continue

            # only add the ancestors if it's not NCBITaxon, as this saves a lot of disk space
            if onto.name == "NCBITaxon":
                onto_dict[onto.name][term_id]["ancestors"] = []
            else:
                onto_dict[onto.name][term_id]["ancestors"] = ancestors

    with gzip.open(output_json_file, "wt") as output_json:
        json.dump(onto_dict, output_json, indent=2)


def _get_ancestors(onto_class: owlready2.entity.ThingClass, ontololgy_name: str) -> List[str]:
    """
    Returns a list of ancestors ids of the given onto class, only returns those belonging to ontology_name,
    it will format the id from the form CL_xxxx to CL:xxxx

    :param owlready2.entity.ThingClass onto_class: the class for which ancestors will be retrieved
    :param str ontololgy_name: only ancestors from this ontology will be kept

    :rtype List[str]
    :return list of ancestors (term ids), it could be empty
    """

    ancestors = []

    for ancestor in onto_class.ancestors():
        if onto_class.name == ancestor.name:
            continue
        if ancestor.name.split("_")[0] == ontololgy_name:
            ancestors.append(ancestor.name.replace("_", ":"))

    return ancestors


# Download and parse owls upon execution
if __name__ == "__main__":
    _download_owls()
    _parse_owls()
