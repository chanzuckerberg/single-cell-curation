import gzip
import os
import sys
from typing import Dict

import gtf_tools

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

                # if the symbol starts with ENSG and it does not match the Ensembl ID, then the symbol used should be the Ensembl ID
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


if __name__ == "__main__":
    _parse_gtf("./temp/mus_musculus.gtf.gz", os.path.join(env.ONTOLOGY_DIR, "genes_mus_musculus.csv.gz"))
    _parse_gtf("./temp/homo_sapiens.gtf.gz", os.path.join(env.ONTOLOGY_DIR, "genes_homo_sapiens.csv.gz"))
    _parse_gtf("./temp/sars_cov_2.gtf.gz", os.path.join(env.ONTOLOGY_DIR, "genes_sars_cov_2.csv.gz"))
    _process_ercc("./temp/ercc.txt", os.path.join(env.ONTOLOGY_DIR, "genes_ercc.csv.gz"))
