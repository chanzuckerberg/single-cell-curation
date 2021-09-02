import os
import gzip
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../cellxgene_schema"))
import env


def _parse_gtf(gtf_path: str, output_file: str):
    """
    Parses a gziped GTF file to get gene and transcript info into a gziped comma-separated file with the following
    structure, with three columns and no header: 1) gene/transcript id, 2) gene/transcript name,
    3) gene/transcript version

    :param str gtf_path: path to gzipped gtf file
    :param str output_json_file: path to output json

    :rtype: None
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
            current_features = _get_features(line)

            # Select  features of interest, raise error feature of interest not found
            target_features = [""] * len(features)
            for i in range(len(features)):
                feature = features[i]
                if feature in current_features:
                    current_features[feature] = current_features[feature].replace(
                        '"', ""
                    )
                    target_features[i] = current_features[feature]

                # Add gene version if available from gene id
                if feature in ["gene_id", "transcript_id"]:
                    if "." in target_features[i]:
                        (feature_id, feature_version) = target_features[i].split(".")
                        target_features[i] = feature_id
                        current_features[feature.replace("id", "version")] = feature_version

            output_to_print += ",".join(target_features) + "\n"

    with gzip.open(output_file, "wt") as output:
        output.write(output_to_print)


def _get_features(gtf_line: str) -> dict:

    """
    Parses the features found in column 8 of GTF, returns a dict with keys as feature names and values as the feature
    values

    :param str gtf_line: a line from a GTF file
    """

    return_features = {}

    for feature in gtf_line[8].split(";"):
        if len(feature) > 1:
            feature = feature.strip().split(" ")
            feature_name = feature[0]
            feature_value = feature[1]
            return_features[feature_name] = feature_value

    return return_features


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
            ercc_id = line.rstrip().split("\t")[0]
            output_to_print += ",".join([ercc_id, ercc_id + " (spike-in control)", "1"]) + "\n"

    with gzip.open(output_file, "wt") as output:
        output.write(output_to_print)


if __name__ == "__main__":
    _parse_gtf("./temp/mus_musculus.gtf.gz", os.path.join(env.ONTOLOGY_DIR, "genes_mus_musculus.csv.gz"))
    _parse_gtf("./temp/homo_sapiens.gtf.gz", os.path.join(env.ONTOLOGY_DIR, "genes_homo_sapiens.csv.gz"))
    _parse_gtf("./temp/sars_cov_2.gtf.gz", os.path.join(env.ONTOLOGY_DIR, "genes_sars_cov_2.csv.gz"))
    _process_ercc("./temp/ercc.txt", os.path.join(env.ONTOLOGY_DIR, "genes_ercc.csv.gz"))
