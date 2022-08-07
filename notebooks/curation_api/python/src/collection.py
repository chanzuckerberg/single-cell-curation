import json
import os
import requests

from src.utils.logger import get_custom_logger, failure, success
from src.utils.http import url_builder, get_headers


logger = get_custom_logger()


def create_collection(collection_form_metadata: dict) -> str:
    """
    Create a new private Collection
    :param collection_form_metadata: the Collection metadata to use to instantiate a Collection
    :return: the Collection id
    """
    url = url_builder("/collections")
    headers = get_headers()
    try:
        res = requests.post(url, data=json.dumps(collection_form_metadata), headers=headers)
        res.raise_for_status()
        data = res.json()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    else:
        collection_id = data.get("collection_id")
        success(logger, f"New private Collection id:\n{collection_id}\n",
                f"New private Collection url:\n{os.getenv('site_url')}/collections/{collection_id}")

        return collection_id


def create_revision(collection_id: str) -> str:

    url = url_builder(f"/collections/{collection_id}/revision")
    headers = get_headers()
    try:
        res = requests.post(url, headers=headers)
        res.raise_for_status()
        data = res.json()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    else:
        revision_id = data.get("revision_id")
        success(logger, f"Revision id:\n{revision_id}\n",
                f"Revision url:\n{os.getenv('site_url')}/collections/{revision_id}")

        return revision_id


def delete_collection(collection_id: str) -> None:

    url = url_builder(f"/collections/{collection_id}")
    headers = get_headers()
    try:
        res = requests.delete(url, headers=headers)
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    else:
        success(logger, f"Deleted the Collection at url:\n{os.getenv('site_url')}/collections/{collection_id}")


def get_collection(collection_id: str) -> dict:

    url = url_builder(f"/collections/{collection_id}")
    headers = get_headers()
    res = requests.get(url, headers=headers)
    res.raise_for_status()
    return res.json()


def get_collections(visibility: str = None) -> list:

    url = url_builder("/collections")
    headers = get_headers()
    params = {"visibility": visibility} if visibility else None
    res = requests.get(url, headers=headers, params=params)
    res.raise_for_status()
    return res.json().get("collections")


def update_collection(collection_id: str, collection_form_metadata: dict) -> None:

    url = url_builder(f"/collections/{collection_id}")
    headers = get_headers()
    res = requests.patch(url, headers=headers, data=json.dumps(collection_form_metadata))
    res.raise_for_status()
    success(logger, f"Updated the Collection at url:\n{os.getenv('site_url')}/collections/{collection_id}")
    return res.json()
