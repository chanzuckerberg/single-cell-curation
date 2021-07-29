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


if __name__ == "__main__":
    _parse_gtf("./mus_musculus.gtf", os.path.join(env.ONTOLOGY_DIR, "genes_mus_musculus.csv.gz"))
    _parse_gtf("./homo_sapiens.gtf", os.path.join(env.ONTOLOGY_DIR, "genes_homo_sapiens.csv.gz"))
