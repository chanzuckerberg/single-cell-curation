import gzip
import json
import os
import sys
import urllib.request
from collections import defaultdict
from threading import Thread
from urllib.error import HTTPError, URLError

import yaml
import warnings
import pronto
warnings.filterwarnings("ignore", category=pronto.warnings.ProntoWarning)
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

        output_file = os.path.join(output_dir, _ontology + ".obo") if _ontology == "NCBITaxon" else os.path.join(output_dir, _ontology + ".owl")
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

    owl_files = dict()
    for owl_file in os.listdir(working_dir):
        if owl_file.endswith(".owl") or owl_file.endswith(".obo"):
            ontology_name = owl_file.split(".")[0]
            owl_files[ontology_name] = os.path.join(working_dir, owl_file)

    # Parse owl files
    onto_dict = {}
    for ontology_name, owl_file in owl_files.items():
        onto = pronto.Ontology(owl_file)
        print(f"Processing {ontology_name}")
        onto_dict[ontology_name] = defaultdict(dict)

        for onto_term in onto.terms():
            # Skip terms that are not direct children from this ontology
            if ontology_name != onto_term.id.split(":")[0]:
                continue
            # Gets label
            onto_dict[ontology_name][onto_term.id]["label"] = "" if onto_term.name is None else onto_term.name

            # Add the "deprecated" status
            onto_dict[ontology_name][onto_term.id]["deprecated"] = False
            if onto_term.obsolete:
                # if deprecated, include information to determine replacement term(s)
                onto_dict[ontology_name][onto_term.id]["deprecated"] = True
                if onto_term.comment:
                    onto_dict[ontology_name][onto_term.id]["comments"] = [onto_term.comment]
                # stores term tracking URL, such as a github issue discussing deprecated term
                for annotation in onto_term.annotations:
                    if annotation.property == "http://purl.obolibrary.org/obo/IAO_0000233":
                        onto_dict[ontology_name][onto_term.id]["term_tracker"] = annotation.literal
                # only need to record replaced_by OR considers
                if onto_term.replaced_by:
                    onto_dict[ontology_name][onto_term.id]["replaced_by"] = next(iter(onto_term.replaced_by.ids))
                else:
                    if onto_term.consider:
                        onto_dict[ontology_name][onto_term.id]["consider"] = list(onto_term.consider.ids)

            # Gets ancestors
            ancestors = _get_ancestors(onto_term, ontology_name)

            # If "children_of" specified in owl info then skip the current term if it is
            # not a children of those indicated.
            if (ontology_name in owl_info and "children_of" in owl_info[ontology_name]) and (
                not list(set(ancestors) & set(owl_info[ontology_name]["children_of"]))
            ):
                onto_dict[ontology_name].pop(onto_term.id)
                continue

            # only add the ancestors if it's not NCBITaxon, as this saves a lot of disk space
            if ontology_name == "NCBITaxon":
                onto_dict[ontology_name][onto_term.id]["ancestors"] = []
            else:
                onto_dict[ontology_name][onto_term.id]["ancestors"] = ancestors

    with gzip.open(output_json_file, "wt") as output_json:
        json.dump(onto_dict, output_json, indent=2)


def _get_ancestors(onto_term: pronto.Term, ontology_name: str) -> List[str]:
    """
    Returns a list of ancestors ids of the given onto class, only returns those belonging to ontology_name

    :param pronto.Term onto_term: the class for which ancestors will be retrieved
    :param str ontology_name: only ancestors from this ontology will be kept

    :rtype List[str]
    :return list of ancestors (term ids), it could be empty
    """
    return [
        ancestor.id for ancestor in onto_term.superclasses(with_self=False)
        if ancestor.id.split(":")[0] == ontology_name
    ]


# Download and parse owls upon execution
if __name__ == "__main__":
    _download_owls()
    _parse_owls(output_json_file=os.path.join(env.ONTOLOGY_DIR, "temp_all_ontology.json.gz"))
