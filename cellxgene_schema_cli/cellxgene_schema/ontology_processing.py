import owlready2
import yaml
import urllib.request
import os
import gzip
import json


def _dowload_owls(
    owl_info_yml=os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "ontology_files/owl_info.yml"
    ),
    output_dir=os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "ontology_files/"
    ),
):

    """Downloads the ontology owl files specified in 'owl_info_yml' into 'output_dir'"""

    with open(owl_info_yml, "r") as owl_info_handle:
        owl_info = yaml.safe_load(owl_info_handle)

    for ontology, info in owl_info.items():

        print(f"Downloading {ontology}")

        # Get owl info
        latest_version = owl_info[ontology]["latest"]
        url = owl_info[ontology]["urls"][latest_version]

        # Format of owl (handles cases where they are compressed)
        download_format = url.split(".")[-1]

        output_file = os.path.join(output_dir, ontology + ".owl")
        if download_format == "gz":
            urllib.request.urlretrieve(url, output_file + ".gz")
            _decompress(output_file + ".gz", output_file)
            os.remove(output_file + ".gz")
        else:
            urllib.request.urlretrieve(url, output_file)


def _decompress(infile, tofile):
    """
    Decompresses a gziped file
    """
    with open(infile, "rb") as inf, open(tofile, "w", encoding="utf8") as tof:
        decom_str = gzip.decompress(inf.read()).decode("utf-8")
        tof.write(decom_str)


def _parse_owls(
    working_dir=os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "ontology_files"
    ),
    owl_info_yml=os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "ontology_files/owl_info.yml"
    ),
    output_json_file=os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "ontology_files/all_ontology.json.gz",
    ),
):

    """
    Parser all owl files in working_dir. Extracts information from all classes in the owl file.
    The extracted information is written into a gzipped a json file with the following structure:
    {
        "ontology_name":
            {
            "term_id": {
                "label": "..."
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
    """

    with open(owl_info_yml, "r") as owl_info_handle:
        owl_info = yaml.safe_load(owl_info_handle)

    owl_files = [
        os.path.join(working_dir, i)
        for i in os.listdir(working_dir)
        if i.endswith(".owl")
    ]

    # Parse owl files
    onto_dict = {}
    for i in owl_files:
        onto = owlready2.get_ontology(i).load()
        onto_dict[onto.name] = {}

        print(f"Processing {onto.name}")
        for onto_class in onto.classes():

            term_id = onto_class.name.replace("_", ":")

            # Skip terms that are not direct children from this ontology
            if not onto.name == term_id.split(":")[0]:
                continue

            # If there are specified target terms then only work with them
            if onto.name in owl_info:
                if "only" in owl_info[onto.name]:
                    if term_id not in owl_info[onto.name]["only"]:
                        continue

            onto_dict[onto.name][term_id] = dict()
            try:
                onto_dict[onto.name][term_id]["label"] = onto_class.label[0]
            except:
                onto_dict[onto.name][term_id]["label"] = ""

            onto_dict[onto.name][term_id]["ancestors"] = [
                i.name.replace("_", ":")
                for i in onto_class.ancestors()
                if i.name.split("_")[0] == onto.name
            ]

    with gzip.open(output_json_file, "wt") as output_json:
        json.dump(onto_dict, output_json, indent=2)


def _parse_gtf(gtf_path, output_file):
    """
    Parses a gziped GTF file to get gene and transcript info into a gziped comma-separated file with the following structure,
    with three columns and no header: 1) gene/transcript id, 2) gene/transcript name, 3) gene/transcript version

    :param str gtf_path: path to gzipped gtf file
    :param str output_json_file: path to output json

    :return: none
    :rtype: none
    """

    output_to_print = ""
    with gzip.open(gtf_path, "rb") as gtf:
        for line in gtf:
            line = line.decode("utf-8")

            if line[0] == "#":
                continue

            line = line.rstrip().split("\t")

            # Desired features based on whether is gene or transcript
            if line[2] == "gene":
                features = ["gene_id", "gene_name", "gene_version"]
            elif line[2] == "transcript":
                features = ["transcript_id", "transcript_name", "transcript_version"]
            else:
                continue

            # Extract features (column 9 of GTF)
            current_features = {
                i.strip().split(" ")[0]: i.strip().split(" ")[1]
                for i in line[8].split(";")
                if len(i) > 1
            }

            # Select  features of interest, raise error feature of interest not found
            target_features = [""] * len(features)
            for i in range(len(features)):
                feature = features[i]
                if feature in current_features:
                    current_features[feature] = current_features[feature].replace(
                        '"', ""
                    )
                    target_features[i] = current_features[feature]

            output_to_print += ",".join(target_features) + "\n"

    with gzip.open(output_file, "wt") as output:
        output.write(output_to_print)


# Download and parse owls upon execution
if __name__ == "__main__":
    _dowload_owls()
    _parse_owls()
