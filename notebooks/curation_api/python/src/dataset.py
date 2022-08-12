import json

import boto3
import logging
import os
import re
import requests
import threading
from botocore.credentials import RefreshableCredentials
from botocore.session import get_session
from typing import Tuple

from src.utils.config import format_c_url
from src.utils.logger import get_custom_logger, failure, success
from src.utils.http import url_builder, get_headers


logger = get_custom_logger()

UUID_REGEX = r"[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}"
ID_REGEX = f"(?P<id>{UUID_REGEX})"
CURATOR_TAG_PREFIX_REGEX = r"(?P<tag_prefix>.*)"
EXTENSION_REGEX = r"(?P<extension>h5ad)"


def get_identifier_type_and_value(identifier: str) -> Tuple[str, str]:
    identifier_type = None
    identifier_value = identifier
    if re.match(f"^{UUID_REGEX}$", identifier):
        # identifier is a uuid
        identifier_type = "dataset_id"
    else:
        # CURATOR_TAG_PREFIX_REGEX is superfluous; leaving in to match lambda handler code; may use later
        matched = re.match(f"({UUID_REGEX}|{CURATOR_TAG_PREFIX_REGEX})\\.{EXTENSION_REGEX}$", identifier)
        if matched:
            matches = matched.groupdict()
            if _id := matches.get("id"):
                identifier_type = "dataset_id"
                identifier_value = _id
            else:
                identifier_type = "curator_tag"

    if not identifier_type:
        raise Exception(f"The identifier '{identifier}' must be either 1) a curator tag that includes a '.h5ad' suffix "
                        f"OR 2) a uuid")

    return identifier_value, identifier_type


def delete_dataset(collection_id: str, identifier: str):
    """
    Delete a private Dataset
    :param collection_id: the id of the Collection to which the Dataset belongs
    :param identifier: the curator tag or Dataset id
    :return: True if deletion is successful otherwise False
    """
    url = url_builder(f"/collections/{collection_id}/datasets")
    headers = get_headers()

    identifier_value, identifier_type = get_identifier_type_and_value(identifier)

    params_dict = dict()
    params_dict[identifier_type] = identifier_value

    success_message = f"Deleted the Dataset with {identifier_type} '{identifier_value}' from its Collection: " \
                      f"\n{format_c_url(collection_id)}"
    try:
        res = requests.delete(url, headers=headers, params=params_dict)
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    else:
        success(logger, success_message)


def get_assets(collection_id: str, identifier: str):
    """
    Fetch download links for assets for a Dataset
    :param collection_id: the id of the Collection to which the Dataset belongs
    :param identifier: the curator tag or Dataset id
    :return: download links
    """
    url = url_builder(f"/collections/{collection_id}/datasets/assets")
    headers = get_headers()

    identifier_value, identifier_type = get_identifier_type_and_value(identifier)
    params_dict = dict()
    params_dict[identifier_type] = identifier_value

    try:
        res = requests.get(url, headers=headers, params=params_dict)
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    return res.json()


def update_curator_tag(collection_id: str, identifier: str, new_tag: str):
    """
    Update a private Dataset's curator tag
    :param collection_id: the id of the Collection to which the Dataset belongs
    :param identifier: the curator tag or Dataset id
    :param new_tag: the new curator tag to assign to the Dataset
    """
    url = url_builder(f"/collections/{collection_id}/datasets")
    headers = get_headers()

    identifier_value, identifier_type = get_identifier_type_and_value(identifier)

    params_dict = dict()
    params_dict[identifier_type] = identifier_value

    new_curator_tag_dict = dict(curator_tag=new_tag)

    success_message = f"Dataset with {identifier_type} '{identifier_value}' updated to have curator_tag '{new_tag}'"
    try:
        res = requests.patch(url, headers=headers, params=params_dict, data=json.dumps(new_curator_tag_dict))
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    else:
        success(logger, success_message)


def upload_datafile_from_link(link: str, collection_id: str, identifier: str = None):
    """
    Create/update a Dataset from the datafile found at the source link.
    :param link: the source datafile link to upload to the Data Portal to become a Dataset
    :param collection_id: the id of the Collection to which the resultant Dataset will belong
    :param identifier: the curator tag or Dataset id. Must be suffixed with '.h5ad'. See heading
    of create_dataset_from_local_file.ipynb for details about how to use the identifier to 'create new' vs 'replace existing'
    """
    url = url_builder(f"/collections/{collection_id}/datasets/upload-link")
    headers = get_headers()

    data_dict = dict(link=link)
    if identifier:
        identifier_value, identifier_type = get_identifier_type_and_value(identifier)

        identifier_param_name = "curator_tag" if identifier_type == "curator_tag" else "id"
        data_dict[identifier_param_name] = identifier_value

        success_message = f"Uploading Dataset with {identifier_type} '{identifier_value}' to Collection " \
                          f"{os.getenv('site_url')}/collections/{collection_id} sourcing from datafile at {link}"
    else:
        success_message = f"Uploading Dataset to Collection {os.getenv('site_url')}/collections/{collection_id} " \
                          f"sourcing from datafile at {link}"

    try:
        res = requests.put(url, headers=headers, data=json.dumps(data_dict))
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    else:
        success(logger, success_message)


def upload_local_datafile(datafile_path: str, collection_id: str, identifier: str):
    """
    :param datafile_path: the fully qualified path of the datafile to be uploaded
    :param collection_id: the id of the Collection to which the resultant Dataset will belong
    :param identifier: the curator tag or Dataset id. Must be suffixed with '.h5ad'. See heading
    of upload_local_datafile.ipynb for details about how to use the identifier to 'create new' vs 'replace existing'
    :param log_level: the logging level
    Datasets.
    :return: None
    """
    url = url_builder(f"/collections/{collection_id}/datasets/s3-upload-credentials")
    headers = get_headers()

    def retrieve_s3_credentials_and_upload_key_prefix():
        try:
            res = requests.get(url, headers=headers)
            res.raise_for_status()
        except requests.HTTPError as e:
            failure(logger, e)
            raise e
        return res.json()

    log_level = os.getenv("log_level", "INFO")  # hack to determine log level in callback (separate process)

    def s3_refreshable_credentials_cb():
        res_data = retrieve_s3_credentials_and_upload_key_prefix()
        s3_credentials = res_data.get("Credentials")
        s3_credentials_formatted = {
            "access_key": s3_credentials.get("AccessKeyId"),
            "secret_key": s3_credentials.get("SecretAccessKey"),
            "token": s3_credentials.get("SessionToken"),
            "expiry_time": s3_credentials.get("Expiration"),
        }
        if getattr(logging, log_level) < 20:  # if log level NOTSET or DEBUG
            print("Retrieved/refreshed s3 credentials")
        return s3_credentials_formatted

    filesize = os.path.getsize(datafile_path)

    def get_progress_cb(collection_id: str, identifier: str):
        lock = threading.Lock()
        uploaded_bytes = 0
        prev_percent = 0

        def progress_cb(num_bytes):
            nonlocal uploaded_bytes
            nonlocal prev_percent
            should_update_progress_printout = False

            lock.acquire()
            uploaded_bytes += num_bytes
            percent_of_total_upload = float("{:.1f}".format(uploaded_bytes / filesize * 100))
            if percent_of_total_upload > prev_percent:
                should_update_progress_printout = True
            prev_percent = percent_of_total_upload
            lock.release()

            if should_update_progress_printout:
                color = "\033[38;5;10m" if percent_of_total_upload == 100 else ""
                if getattr(logging, log_level) < 40:
                    print(f"{collection_id}/{identifier}: "
                          f"\033[1m{color}{percent_of_total_upload}% uploaded\033[0m\r", end="")

        return progress_cb

    credentials_and_path = retrieve_s3_credentials_and_upload_key_prefix()
    bucket, key_prefix = credentials_and_path["Bucket"], credentials_and_path["UploadKeyPrefix"]
    upload_key = key_prefix + identifier
    logger.debug(f"Full S3 write path is s3://{bucket}/{upload_key}\n")

    session_creds = RefreshableCredentials.create_from_metadata(
        metadata=s3_refreshable_credentials_cb(),
        refresh_using=s3_refreshable_credentials_cb,
        method="sts-assume-role-with-web-identity",
    )
    session = get_session()
    session._credentials = session_creds
    boto3_session = boto3.Session(botocore_session=session)
    s3 = boto3_session.client("s3")

    try:
        logger.info(f"\nUploading {datafile_path} to Collection {collection_id} with identifier '{identifier}'...\n")
        s3.upload_file(
            Filename=datafile_path,
            Bucket=bucket,
            Key=upload_key,
            Callback=get_progress_cb(collection_id, identifier),
        )
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    else:
        success(logger, "UPLOAD COMPLETE -- Dataset is queued for processing and will surface in the system shortly.")
