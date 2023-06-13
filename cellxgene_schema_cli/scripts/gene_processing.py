import csv
import gzip
import os
import sys
import urllib.request
from multiprocessing import Process
from typing import Dict

import gtf_tools
import yaml

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../cellxgene_schema"))
import env


def _parse_gtf(gtf_path: str, output_file: str):
    """
    Parses a gziped GTF file to get gene and transcript info into a gziped comma-separated file with the following
    structure, with three columns and no header: 1) gene id, 2) gene,
    3) gene

    :param str gtf_path: path to gzipped gtf file
    :param str output_json_file: path to output json

    :rtype: None
    """

    output_to_print = ""
    gene_lengths = _get_gene_lengths_from_gtf(gtf_path)
    with gzip.open(gtf_path, "rb") as gtf:
        for line in gtf:
            line = line.decode("utf-8")

            if line[0] == "#":
                continue

            line = line.rstrip().split("\t")

            # Desired features based on whether is gene or transcript
            if line[2] == "gene":
                features = ["gene_id", "gene_name", "gene_version"]
            else:
                continue

            # Extract features (column 9 of GTF)
            current_features = gtf_tools._get_features(line)
            # Filter genes suffixed with "PAR_Y"
            if current_features["gene_id"].endswith("PAR_Y"):
                continue

            # get gene length
            current_length = gene_lengths[current_features["gene_id"]]

            # Select  features of interest, raise error if feature of interest not found
            target_features = [""] * len(features)
            for i in range(len(features)):
                feature = features[i]
                if feature in current_features:
                    target_features[i] = current_features[feature]

                # if the symbol starts with ENSG and it does not match the Ensembl ID, then the symbol used should be
                # the Ensembl ID
                if (
                    feature in ["gene_name"]
                    and current_features[feature].startswith("ENSG")
                    and current_features[feature] != current_features["gene_id"]
                ):
                    target_features[i] = current_features["gene_id"]

                # Add gene version if available from gene id
                if feature in ["gene_id"]:
                    if "." in target_features[i]:
                        (feature_id, feature_version) = target_features[i].split(".")
                    else:
                        feature_id = target_features[i]
                        feature_version = ""

                    target_features[i] = feature_id
                    current_features[feature.replace("id", "version")] = feature_version

            output_to_print += ",".join(target_features + [str(current_length)]) + "\n"

    with gzip.open(output_file, "wt") as output:
        output.write(output_to_print)


def _get_gene_lengths_from_gtf(gtf_path: str) -> Dict[str, int]:
    """
    Parses a GTF file and calculates gene lengths, which are calculated as follows for each gene:
       1. Get all different isoforms
       2. Merge exons to create a set of non-overlapping "meta" exons
       3. Sum the lengths of these "meta" exons

    Code inspired from http://www.genemine.org/gtftools.php

    :param str gtf_path: path to gzipped gtf file

    :rtype  Dict[str]
    :return A dictionary with keys being gene ids and values the corresponding length in base pairs
    """

    with gzip.open(gtf_path, "rb") as gtf:
        # Dictionary of list of tuples that will store exon in bed format-like. Elements of the tuples will be the
        # equivalent  fields from the bed format: chromosome, start, end, strand). Each list of tuples will correspond
        # to one gene.

        exons_in_bed = {}

        for line in gtf:
            line = line.decode("utf-8")

            if line[0] == "#":
                continue

            line = line.rstrip().split("\t")

            if line[2] != "exon":
                continue

            # Convert line to bed-like format: (chromosome, start, end, strand)
            exon_bed = (line[0], int(line[3]) - 1, int(line[4]), line[6])

            current_features = gtf_tools._get_features(line)
            gene_id = current_features["gene_id"]
            if gene_id in exons_in_bed:
                exons_in_bed[gene_id].append(exon_bed)
            else:
                exons_in_bed[gene_id] = [exon_bed]

    # Merge exons from the same gene to create non-overlapping "meta" genes
    # Then calculate gene length
    gene_lengths = {}
    for gene in exons_in_bed:
        meta_exons = gtf_tools.merge_bed_ranges(exons_in_bed[gene])

        # get length for this gene, i.e. sum of lengths of "meta" exons
        gene_lengths[gene] = 0
        for exon in meta_exons:
            gene_lengths[gene] += exon[2] - exon[1]

    return gene_lengths


def _process_ercc(ercc_path: str, output_file: str):
    """
    process the ERCC download, keeps only first column with no header

    :param str gtf_path: path to ercch download from thermo fisher
    :param str output_json_file: path to output file
    """

    output_to_print = ""
    with open(ercc_path, "r") as ercc:
        lines = ercc.readlines()[1:]
        for line in lines:
            line = line.rstrip().split("\t")
            ercc_id = line[0]
            errc_length = str(len(line[4]))
            errc_version = "1"
            output_to_print += ",".join([ercc_id, ercc_id + " (spike-in control)", errc_version, errc_length]) + "\n"

    with gzip.open(output_file, "wt") as output:
        output.write(output_to_print)


def process_gene_info(gene_info: dict) -> None:
    """
    Download the gene_info and convert it into a csv
    :param gene_info: context to download and process the gene information.
    """
    print("download", gene_info["description"])
    temp_file_path, _ = urllib.request.urlretrieve(gene_info["url"])
    previous_ref_filepath = os.path.join(env.ONTOLOGY_DIR, f"previous_{gene_info['description']}.csv.gz")
    # temporarily backup previous processed csv
    os.rename(gene_info["new_file"], previous_ref_filepath)
    try:
        print("process", gene_info["description"])
        gene_info["processor"](temp_file_path, gene_info["new_file"])
        print("generating gene reference diff for", gene_info["description"])
        generate_gene_ref_diff(gene_info["description"], gene_info["new_file"], previous_ref_filepath)
    except Exception as e:
        print("processing failed, reverting", gene_info["description"])
        os.replace(previous_ref_filepath, gene_info["new_file"])
        raise e
    # remove previous processed csv once diff is complete
    os.remove(previous_ref_filepath)
    print("finish", gene_info["description"])


def generate_gene_ref_diff(output_filename: str, current_ref_filepath: str, previous_ref_filepath: str) -> None:
    """
    Compare the previous gene reference CSV to the newest generated CSV. Report a text file with every
    Ensembl ID available in the previous gene reference CSV that is no longer found in the newest.
    :param output_filename: Filename for output diff textfile
    :param current_ref_filepath: Full filepath for the new processed gene csv
    :param previous_ref_filepath: Full filepath for the previous processed gene csv
    """
    new_ref_gene_ids = set()
    with gzip.open(current_ref_filepath, "rt") as f:
        for row in csv.reader(f):
            gene_id = row[0]
            new_ref_gene_ids.add(gene_id)

    removed_gene_ids = []
    with gzip.open(previous_ref_filepath, "rt") as f:
        for row in csv.reader(f):
            gene_id = row[0]
            if gene_id not in new_ref_gene_ids:
                removed_gene_ids.append(gene_id)

    diff_filepath = os.path.join(env.ONTOLOGY_DIR, f"{output_filename}_diff.txt")
    if removed_gene_ids:
        with open(diff_filepath, "w") as f:
            for gene_id in removed_gene_ids:
                f.write(f"{gene_id}\n")
    else:
        print(f"No genes removed from {output_filename} reference. No diff created.")


def main():
    with open(env.GENE_INFO_YAML, "r") as gene_info_handle:
        gene_infos: dict = yaml.safe_load(gene_info_handle)

    jobs = []
    for gene_info in gene_infos.values():
        gene_info["new_file"] = os.path.join(env.ONTOLOGY_DIR, f"genes_{gene_info['description']}.csv.gz")
        # determine how to process based on file type
        if gene_info["url"].endswith("gtf.gz"):
            gene_info["processor"] = _parse_gtf
        elif gene_info["url"].endswith("txt"):
            gene_info["processor"] = _process_ercc
        else:
            raise TypeError(f"unknown file type: {gene_info['file_type']}")

        # add version to URL if needed
        if gene_info.get("version"):
            gene_info["url"] = gene_info["url"].format(version=gene_info["version"])
        job = Process(target=process_gene_info, args=(gene_info,))
        job.start()
        jobs.append(job)

    for job in jobs:
        job.join()


if __name__ == "__main__":
    main()
