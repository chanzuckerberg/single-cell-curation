import json
import os
import requests

from src.utils.logger import get_custom_logger, failure, success
from src.utils.http import url_builder, get_headers


logger = get_custom_logger()


def create_collection(collection_form_metadata: str) -> str:
    """
    Create a new private Collection
    :param collection_form_metadata: the Collection metadata to use to instantiate a Collection
    :return: the Collection uuid
    """
    url = url_builder("/collections")
    headers = get_headers()
    try:
        res = requests.post(url, data=json.dumps(collection_form_metadata), headers=headers)
        res.raise_for_status()
        data = res.json()
    except Exception as e:
        failure(logger, e)
    else:
        collection_uuid = data.get("collection_uuid")
        logger.info("\n\033[1m\033[38;5;10mSUCCESS\033[0m\n")  # 'SUCCESS' in bold green
        logger.info(f"New private Collection uuid:\n{collection_uuid}\n")
        logger.info(f"New private Collection url:\n{os.getenv('site_url')}/collections/{collection_uuid}")
        return collection_uuid


def create_revision(collection_uuid: str) -> str:

    url = url_builder(f"/collections/{collection_uuid}/revision")
    headers = get_headers()
    try:
        res = requests.post(url, headers=headers)
        res.raise_for_status()
        data = res.json()
    except Exception as e:
        failure(logger, e)
    else:
        revision_uuid = data.get("revision_id")
        logger.info("\n\033[1m\033[38;5;10mSUCCESS\033[0m\n")  # 'SUCCESS' in bold green
        logger.info(f"Revision uuid:\n{revision_uuid}\n")
        logger.info(f"Revision url:\n{os.getenv('site_url')}/collections/{revision_uuid}")
        return revision_uuid


def delete_collection(collection_uuid: str) -> None:

    url = url_builder(f"/collections/{collection_uuid}")
    headers = get_headers()
    try:
        res = requests.delete(url, headers=headers)
        res.raise_for_status()
    except Exception as e:
        failure(logger, e)
    else:
        success(logger, f"Deleted the Collection at url:\n{os.getenv('site_url')}/collections/{collection_uuid}")


def get_collection(collection_uuid: str) -> dict:

    url = url_builder(f"/collections/{collection_uuid}")
    headers = get_headers()
    res = requests.get(url, headers=headers)
    res.raise_for_status()
    return res.json()


def get_collections(visibility: str = "PUBLIC") -> list:

    url = url_builder("/collections")
    headers = get_headers()
    res = requests.get(url, headers=headers, params={"visibility": visibility})
    res.raise_for_status()
    return res.json().get("collections")


def update_collection(collection_uuid: str, collection_form_metadata: dict) -> None:

    url = url_builder(f"/collections/{collection_uuid}")
    headers = get_headers()
    try:
        res = requests.put(url, headers=headers, data=json.dumps(collection_form_metadata))
        res.raise_for_status()
    except Exception as e:
        failure(logger, e)
    else:
        success(logger, f"Updated the Collection at url:\n{os.getenv('site_url')}/collections/{collection_uuid}")
