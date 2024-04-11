import csv
import gzip
import hashlib
import os
import sys
import urllib.request
from typing import Dict

import gtf_tools
import yaml

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../cellxgene_schema"))
import env


class GeneProcessingResult:
    def __init__(self, gene_id: str, gene_name: str, gene_version: str, gene_length: str):
        self.gene_id = gene_id  # ex: ENSG00000141510, should be unique
        self.gene_name = gene_name  # ex: TP53, not necessarily unique
        self.gene_version = gene_version
        self.gene_length = gene_length


class GeneProcessor:
    def __init__(self):
        # Global mapping of gene id to all the metadata we need to generate output files
        self.gene_metadata: dict[str, GeneProcessingResult] = {}
        # Mapping from description (ex: "mus musculus") to a list of gene ids
        self.gene_ids_by_description: dict[str, list[str]] = {}

    def write_gzip(self, data: str, output_filename: str):
        """
        Writes data to a gziped file. The date modified is not written to the gzip file. This allows
        comparing the hashed contents of two gziped files to determine if they are the same.

        :param str data: data to write
        :param str output_file: path to output file

        :rtype: None
        """
        with open(output_filename, "wb") as fileobj, gzip.GzipFile(mode="wb", fileobj=fileobj, mtime=0) as myzip:
            myzip.write(data.encode("utf-8"))

    def digest(self, file_name: str) -> str:
        with open(file_name, "rb") as f:
            # Read the contents of the file in chunks
            chunk_size = 1024
            hasher = hashlib.sha256()
            while chunk := f.read(chunk_size):
                hasher.update(chunk)
        return hasher.hexdigest()

    def _parse_gtf(self, gtf_path: str, gene_info_description: str):
        """
        Parses a gziped GTF file to get gene and transcript info into a gziped comma-separated file with the following
        structure, with three columns and no header: 1) gene id, 2) gene,
        3) gene

        :param str gtf_path: path to gzipped gtf file
        :param str output_json_file: path to output json

        :rtype: None
        """

        gene_lengths = self._get_gene_lengths_from_gtf(gtf_path)
        with gzip.open(gtf_path, "rb") as gtf:
            for byte_line in gtf:
                line = byte_line.decode("utf-8")

                if line[0] == "#":
                    continue

                line = line.rstrip().split("\t")  # type: ignore

                # Desired features based on whether is gene or transcript
                if line[2] == "gene":
                    features = ["gene_id", "gene_name", "gene_version"]
                else:
                    continue

                # Extract features (column 9 of GTF)
                current_features = gtf_tools._get_features(line)  # type: ignore
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

                gene_id = target_features[0]
                self.gene_metadata[gene_id] = GeneProcessingResult(
                    gene_id=target_features[0],
                    gene_name=target_features[1],
                    gene_version=target_features[2],
                    gene_length=str(current_length),
                )
                if gene_info_description in self.gene_ids_by_description:
                    self.gene_ids_by_description[gene_info_description].append(gene_id)
                else:
                    self.gene_ids_by_description[gene_info_description] = [gene_id]

    def _get_gene_lengths_from_gtf(self, gtf_path: str) -> Dict[str, int]:
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

            exons_in_bed = {}  # type: ignore

            for byte_line in gtf:
                line = byte_line.decode("utf-8")

                if line[0] == "#":
                    continue

                line = line.rstrip().split("\t")  # type: ignore

                if line[2] != "exon":
                    continue

                # Convert line to bed-like format: (chromosome, start, end, strand)
                exon_bed = (line[0], int(line[3]) - 1, int(line[4]), line[6])

                current_features = gtf_tools._get_features(line)  # type: ignore
                gene_id = current_features["gene_id"]
                if gene_id in exons_in_bed:
                    exons_in_bed[gene_id].append(exon_bed)
                else:
                    exons_in_bed[gene_id] = [exon_bed]

        # Merge exons from the same gene to create non-overlapping "meta" genes
        # Then calculate gene length
        gene_lengths = {}
        for gene in exons_in_bed:
            meta_exons = gtf_tools.merge_bed_ranges(exons_in_bed[gene])  # type: ignore

            # get length for this gene, i.e. sum of lengths of "meta" exons
            gene_lengths[gene] = 0
            for exon in meta_exons:
                gene_lengths[gene] += exon[2] - exon[1]

        return gene_lengths

    def _process_ercc(self, ercc_path: str, gene_info_description: str):
        """
        process the ERCC download, keeps only first column with no header

        :param str gtf_path: path to ercch download from thermo fisher
        :param str output_json_file: path to output file
        """
        with open(ercc_path, "r") as ercc:
            lines = ercc.readlines()[1:]
            for line in lines:
                line = line.rstrip().split("\t")  # type: ignore
                ercc_id = line[0]
                errc_length = str(len(line[4]))
                errc_version = "1"

                self.gene_metadata[ercc_id] = GeneProcessingResult(
                    gene_id=ercc_id,
                    gene_name=ercc_id + " (spike-in control)",
                    gene_version=errc_version,
                    gene_length=errc_length,
                )
                if gene_info_description in self.gene_ids_by_description:
                    self.gene_ids_by_description[gene_info_description].append(ercc_id)
                else:
                    self.gene_ids_by_description[gene_info_description] = [ercc_id]

    def process_gene_infos(self, gene_infos: dict) -> None:
        for gene_info_key in gene_infos:
            # Add to self.gene_labels and self.gene_ids_by_description
            self.process_individual_gene_info(gene_infos[gene_info_key])

        # Deduplicate gene names across all gene_metadata
        print("Deduplicating gene names...")
        gene_names = set()
        duplicated_gene_names = set()
        for gene_metadata in self.gene_metadata.values():
            gene_name = gene_metadata.gene_name
            if gene_name in gene_names:
                duplicated_gene_names.add(gene_name)
            else:
                gene_names.add(gene_name)
        for gene_id, gene_metadata in self.gene_metadata.items():
            gene_name = gene_metadata.gene_name
            if gene_name in duplicated_gene_names:
                self.gene_metadata[gene_id].gene_name = gene_name + "_" + gene_id

        # Write output for each file, and process file diffs
        for gene_info_key in gene_infos:
            gene_info_description = gene_infos[gene_info_key]["description"]
            print("Writing output for ", gene_info_description)
            new_file = os.path.join(env.GENCODE_DIR, f"new_genes_{gene_info_description}.csv.gz")
            previous_ref_filepath = os.path.join(env.GENCODE_DIR, f"genes_{gene_info_description}.csv.gz")
            gene_ids = self.gene_ids_by_description[gene_info_description]
            output_to_print = ""
            try:
                for gene_id in gene_ids:
                    gene_metadata = self.gene_metadata[gene_id]
                    output_to_print += (
                        ",".join(
                            [
                                gene_metadata.gene_id,
                                gene_metadata.gene_name,
                                gene_metadata.gene_version,
                                gene_metadata.gene_length,
                            ]
                        )
                        + "\n"
                    )
                self.write_gzip(output_to_print, new_file)

                if self.digest(new_file) == self.digest(previous_ref_filepath):
                    print("New gene reference is identical to previous gene reference", gene_info_description)
                    os.remove(new_file)
                else:
                    print("generating gene reference diff for", gene_info_description)
                    self.generate_gene_ref_diff(gene_info_description, new_file, previous_ref_filepath)
                    os.replace(new_file, previous_ref_filepath)
            except Exception as e:
                print("Writing to new file failed. Using previous version.", gene_info_description)
                raise e

    def process_individual_gene_info(self, gene_info: dict) -> None:
        """
        Download the gene_info and convert it into a csv
        :param gene_info: context to download and process the gene information.
        """
        gene_info_description = gene_info["description"]
        print("Start", gene_info_description)
        # determine how to process based on file type
        if gene_info["url"].endswith("gtf.gz"):
            processer = self._parse_gtf
        elif gene_info["url"].endswith("txt"):
            processer = self._process_ercc  # type: ignore
        else:
            raise TypeError(f"unknown file type: {gene_info['file_type']}")

        # add version to URL if needed
        url = gene_info["url"].format(version=gene_info["version"]) if gene_info.get("version") else gene_info["url"]

        print("download", gene_info_description)
        temp_file_path, _ = urllib.request.urlretrieve(url)
        try:
            print("process", gene_info_description)
            processer(temp_file_path, gene_info_description)
        except Exception as e:
            print("processing failed", gene_info_description)
            raise e
        print("finish", gene_info_description)

    def generate_gene_ref_diff(
        self, output_filename: str, current_ref_filepath: str, previous_ref_filepath: str
    ) -> None:
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

        diff_filepath = os.path.join(env.GENCODE_DIR, f"{output_filename}_diff.txt")
        if removed_gene_ids:
            with open(diff_filepath, "w") as f:
                for gene_id in removed_gene_ids:
                    f.write(f"{gene_id}\n")
        else:
            print(f"No genes removed from {output_filename} reference. No diff created.")


def main():
    with open(env.GENE_INFO_YAML, "r") as gene_info_handle:
        gene_infos: dict = yaml.safe_load(gene_info_handle)

    gene_processor = GeneProcessor()
    gene_processor.process_gene_infos(gene_infos)


if __name__ == "__main__":
    main()
