import json
import os
import requests
import typing

from src.utils.logger import get_custom_logger


def create_collection(collection_form_metadata: str, log_level: typing.Union[str, int] = "INFO") -> str:

    logger = get_custom_logger(log_level)

    collections_path = "/curation/v1/collections"
    collections_url = f"{os.getenv('api_url_base')}{collections_path}"
    headers = {"Authorization": f"Bearer {os.getenv('access_token')}", "Content-Type": "application/json"}
    try:
        res = requests.post(collections_url, data=json.dumps(collection_form_metadata), headers=headers).json()
    except Exception as e:
        logger.error("\n\033[1m\033[38;5;9mFAILED\033[0m")  # 'FAILED' in bold red
        raise e
    else:
        collection_uuid = res.get("collection_uuid")
        logger.info("\n\033[1m\033[38;5;10mSUCCESS\033[0m\n")  # 'SUCCESS' in bold green
        logger.info(f"New private Collection uuid:\n{collection_uuid}\n")
        logger.info(f"New private Collection url:\n{os.getenv('site_url')}/collections/{collection_uuid}")
        return collection_uuid


def create_revision(collection_uuid: str, log_level: typing.Union[str, int] = "INFO") -> str:

    logger = get_custom_logger(log_level)

    revision_path = f"/curation/v1/collections/{collection_uuid}/revision"
    revision_url = f"{os.getenv('api_url_base')}{revision_path}"
    headers = {"Authorization": f"Bearer {os.getenv('access_token')}", "Content-Type": "application/json"}
    try:
        res = requests.post(revision_url, headers=headers).json()
    except Exception as e:
        logger.error("\n\033[1m\033[38;5;9mFAILED\033[0m")  # 'FAILED' in bold red
        raise e
    else:
        revision_uuid = res.get("revision_id")
        logger.info("\n\033[1m\033[38;5;10mSUCCESS\033[0m\n")  # 'SUCCESS' in bold green
        logger.info(f"Revision uuid:\n{revision_uuid}\n")
        logger.info(f"Revision url:\n{os.getenv('site_url')}/collections/{revision_uuid}")
        return revision_uuid
