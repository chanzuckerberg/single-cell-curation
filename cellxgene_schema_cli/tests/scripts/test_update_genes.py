from unittest.mock import patch

import yaml
from scripts.update_genes import get_latest_release_from_gencode, update_gene_info


def test_no_version(tmp_path):
    """
    When gene_info.yml has no version,
    then it is not updated.
    """
    test_data = {"species": {"description": "test"}}
    expected = {"species": {"description": "test"}}
    tempfile = f"{tmp_path}/gene_info.yml"
    with patch("scripts.update_genes.env") as env:
        env.GENE_INFO_YAML = tempfile
        with open(tempfile, "w") as fp:
            yaml.dump(test_data, fp)
        update_gene_info()
        with open(tempfile, "r") as fp:
            actual = yaml.safe_load(fp)
        assert actual == expected


def test_at_latest_version(tmp_path):
    """
    When the version in gene_info.yml match the latest version from Gencode,
    then no changes are made to gene_info.yml
    """
    test_data = {"species": {"description": "test", "version": "lastest"}}
    expected = {"species": {"description": "test", "version": "lastest"}}
    tempfile = f"{tmp_path}/gene_info.yml"
    with patch("scripts.update_genes.env") as env, patch(
        "scripts.update_genes.get_latest_release_from_gencode"
    ) as get_latest_release_from_gencode:
        get_latest_release_from_gencode.return_value = "lastest"
        env.GENE_INFO_YAML = tempfile
        with open(tempfile, "w") as fp:
            yaml.dump(test_data, fp)
        update_gene_info()
        with open(tempfile, "r") as fp:
            actual = yaml.safe_load(fp)
        assert actual == expected


def test_update_to_latest_version(tmp_path):
    """
    When the version in gene_info.yml is older than the lastest version from Gencode,
    then the version in gene_info.yml is updated to lastest.
    """
    test_data = {"species": {"description": "test", "version": "old"}}
    expected = {"species": {"description": "test", "version": "lastest"}}
    tempfile = f"{tmp_path}/gene_info.yml"
    with patch("scripts.update_genes.env") as env, patch(
        "scripts.update_genes.get_latest_release_from_gencode"
    ) as get_latest_release_from_gencode:
        get_latest_release_from_gencode.return_value = "lastest"
        env.GENE_INFO_YAML = tempfile
        with open(tempfile, "w") as fp:
            yaml.dump(test_data, fp)
        update_gene_info()
        with open(tempfile, "r") as fp:
            actual = yaml.safe_load(fp)
        assert actual == expected


def test_get_latest_release_from_gencode(tmp_path):
    """
    When the version
    """
    with patch("scripts.update_genes.ftplib") as ftplib:
        ftplib.FTP.return_value.nlst.return_value = ["Gencode_mouse.vM25.primary_assembly.annotation.gtf.gz"]
        actual = get_latest_release_from_gencode("mouse")
        assert actual == "M25"
