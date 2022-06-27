import json
import os
import requests

from src.utils.logger import get_custom_logger


logger = get_custom_logger()


def url_builder(path_segment):
    access_token = os.getenv("access_token")
    if not access_token:
        logger.warning("Access token is not set!")
    api_url_base = os.getenv("api_url_base")
    if not api_url_base:
        raise Exception("The required 'api_url_base' env var is not set. Call set_api_urls() from src.utils.config")
    route_path = f"/curation/v1{path_segment}"
    logger.debug(f"route path: {route_path}")
    url = f"{api_url_base}{route_path}"
    logger.debug(f"url: {url}")
    headers = {"Authorization": f"Bearer {access_token}", "Content-Type": "application/json"}
    return url, headers


def create_collection(collection_form_metadata: str) -> str:
    """
    Create a new private Collection
    :param collection_form_metadata: the Collection metadata to use to instantiate a Collection
    :return: the Collection uuid
    """
    url, headers = url_builder("/collections")
    try:
        res = requests.post(url, data=json.dumps(collection_form_metadata), headers=headers)
        res.raise_for_status()
        data = res.json()
    except Exception as e:
        logger.error("\n\033[1m\033[38;5;9mFAILED\033[0m")  # 'FAILED' in bold red
        raise e
    else:
        collection_uuid = data.get("collection_uuid")
        logger.info("\n\033[1m\033[38;5;10mSUCCESS\033[0m\n")  # 'SUCCESS' in bold green
        logger.info(f"New private Collection uuid:\n{collection_uuid}\n")
        logger.info(f"New private Collection url:\n{os.getenv('site_url')}/collections/{collection_uuid}")
        return collection_uuid


def create_revision(collection_uuid: str) -> str:

    url, headers = url_builder(f"/collections/{collection_uuid}/revision")
    try:
        res = requests.post(url, headers=headers)
        res.raise_for_status()
        data = res.json()
    except Exception as e:
        logger.error("\n\033[1m\033[38;5;9mFAILED\033[0m")  # 'FAILED' in bold red
        raise e
    else:
        revision_uuid = data.get("revision_id")
        logger.info("\n\033[1m\033[38;5;10mSUCCESS\033[0m\n")  # 'SUCCESS' in bold green
        logger.info(f"Revision uuid:\n{revision_uuid}\n")
        logger.info(f"Revision url:\n{os.getenv('site_url')}/collections/{revision_uuid}")
        return revision_uuid


def get_collection_uuid(collection_uuid: str) -> dict:

    url, headers = url_builder(f"/collections/{collection_uuid}")
    res = requests.get(url, headers=headers)
    res.raise_for_status()
    return res.json()


def get_collections() -> dict:

    url, headers = url_builder("/collections")
    res = requests.get(url, headers=headers)
    res.raise_for_status()
    return res.json()
