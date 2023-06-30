import os
from typing import List

import semver
import yaml

from . import __version__, env


def get_schema_file_path(version: str) -> str:
    """
    Given a schema version, returns the potential path for the corresponding yaml file of its definition
    :param str version: Schema version
    :return Path to yaml files
    :rtype str
    """

    return os.path.join(env.SCHEMA_DEFINITIONS_DIR, version.replace(".", "_") + ".yaml")


def get_schema_versions_supported() -> List[str]:
    """
    Retrieves a list of the schema versions supported by this version of the validator
    :param str version: Schema version
    :return list of supported schema versions
    :rtype list[str]
    """

    versions = []
    for file in os.listdir(env.SCHEMA_DEFINITIONS_DIR):
        if "yaml" in file:
            version = file.replace("_", ".")
            version = version.replace(".yaml", "")
            versions.append(version)
    return versions


def get_schema_definition(version: str) -> dict:
    """
    Look up and read a schema definition based on a version number like "2.0.0".
    :param str version: Schema version
    :return The schema definition
    :rtype dict
    """

    path = get_schema_file_path(version)

    if not os.path.isfile(path):
        raise ValueError(f"No definition for version '{version}' found.")
    with open(path) as fp:
        return yaml.load(fp, Loader=yaml.FullLoader)


def get_current_schema_version() -> str:
    current_version: semver.Version = semver.Version.parse(__version__)
    return f"{str(current_version.major)}.{str(current_version.minor)}.0"
