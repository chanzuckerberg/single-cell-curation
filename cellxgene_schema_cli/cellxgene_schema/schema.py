import semver
import yaml

from . import __version__, env


def get_schema_definition() -> dict:
    """
    Look up and read the current schema definition
    :return The schema definition
    :rtype dict
    """

    path = env.SCHEMA_DEFINITION_FILE
    with open(path) as fp:
        return yaml.load(fp, Loader=yaml.FullLoader)


def get_current_schema_version() -> str:
    current_version: semver.Version = semver.Version.parse(__version__)
    return f"{str(current_version.major)}.{str(current_version.minor)}.0"
