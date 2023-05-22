#!/usr/bin/env python
import ftplib
import os
import sys

import yaml

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../cellxgene_schema"))
import env


def get_latest_release_from_gencode(species: str) -> str:
    """
    return a URL to the latest release of annotations from https://www.gencodegenes.org/ for mouse or human.
    :param species:mouse or human
    :return:
    """

    server = "ftp.ebi.ac.uk"
    ftp = ftplib.FTP(server)
    ftp.login()
    ftp.cwd(f"/pub/databases/gencode/Gencode_{species}/latest_release/")
    file_name = [*filter(lambda x: x.endswith(".primary_assembly.annotation.gtf.gz"), ftp.nlst())][0]
    release_version = file_name.split(".")[1][1:]
    return release_version


def update_gene_info():
    """
    Get the current version and compare it to the latest. Update the version in the yaml files
    :return:
    """
    with open(env.GENE_INFO_YAML, "r") as gene_info_handle:
        gene_infos: dict = yaml.safe_load(gene_info_handle)

    gene_info_updated = False
    for key, gene_info in gene_infos.items():
        if gene_info.get("version") is None:
            print(f"{key} has no version")
            continue
        latest_version = get_latest_release_from_gencode(key)
        if gene_info["version"] == latest_version:
            print(f"{key} is up to date with latest: {gene_info['version']}")
        else:
            print(f"updating {key}: {gene_info['version']} -> {latest_version}")
            gene_info["version"] = latest_version
            gene_info_updated = True

    if gene_info_updated:
        with open(env.GENE_INFO_YAML, "w") as gene_info_handle:
            yaml.dump(gene_infos, gene_info_handle)


if __name__ == "__main__":
    update_gene_info()
