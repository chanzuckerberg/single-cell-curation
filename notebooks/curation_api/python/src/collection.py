import requests
from src.utils.config import format_c_url
from src.utils.http import get_headers_and_cookies, url_builder
from src.utils.logger import failure, get_custom_logger, success

logger = get_custom_logger()


def create_collection(collection_form_metadata: dict) -> str:
    """
    Create a new private Collection
    :param collection_form_metadata: the Collection metadata to use to instantiate a Collection
    :return: the Collection id
    """
    url = url_builder("/collections")
    try:
        res = requests.post(url, json=collection_form_metadata, **get_headers_and_cookies())
        res.raise_for_status()
        data = res.json()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    collection_id = data.get("collection_id")
    success(
        logger,
        f"New private Collection id:\n{collection_id}\n",
        f"New private Collection url:\n{format_c_url(collection_id)}",
    )
    return collection_id


def create_revision(collection_id: str) -> str:
    url = url_builder(f"/collections/{collection_id}/revision")
    try:
        res = requests.post(url, **get_headers_and_cookies())
        res.raise_for_status()
        data = res.json()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    revision_id = data.get("collection_id")
    success(logger, f"Revision id:\n{revision_id}\n", f"Revision url:\n{format_c_url(revision_id)}")
    return revision_id


def delete_collection(collection_id: str) -> None:
    url = url_builder(f"/collections/{collection_id}")
    try:
        res = requests.delete(url, **get_headers_and_cookies())
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    success(logger, f"Deleted the Collection at url:\n{format_c_url(collection_id)}")


def get_collection(collection_id: str) -> dict:
    url = url_builder(f"/collections/{collection_id}")
    try:
        res = requests.get(url, **get_headers_and_cookies())
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    return res.json()


def get_collection_version(collection_version_id: str) -> dict:
    url = url_builder(f"/collection_versions/{collection_version_id}")
    try:
        res = requests.get(url, **get_headers_and_cookies())
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    return res.json()


def get_collection_versions(collection_id: str) -> list:
    url = url_builder(f"/collections/{collection_id}/versions")
    try:
        res = requests.get(url, **get_headers_and_cookies())
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    return res.json()


def get_collections(visibility: str = None, curator: str = None) -> list:
    params = {}
    if visibility:
        params["visibility"] = visibility
    if curator:
        params["curator"] = curator
    url = url_builder("/collections")
    try:
        res = requests.get(url, params=params, **get_headers_and_cookies())
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    return res.json()


def update_collection(collection_id: str, collection_form_metadata: dict) -> None:
    url = url_builder(f"/collections/{collection_id}")
    try:
        res = requests.patch(url, json=collection_form_metadata, **get_headers_and_cookies())
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    success(logger, f"Updated the Collection at url:\n{format_c_url(collection_id)}")
    return res.json()
