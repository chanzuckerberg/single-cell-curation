import json

import boto3
import logging
import os
import re
import requests
import threading
from botocore.credentials import RefreshableCredentials
from botocore.session import get_session
from datetime import datetime, timezone
from typing import Tuple

from src.utils.logger import get_custom_logger
from src.utils.http import url_builder, get_headers


logger = get_custom_logger()

UUID_REGEX = r"[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}"
DATASET_ID_REGEX = f"(?P<dataset_uuid>{UUID_REGEX})"
CURATOR_TAG_PREFIX_REGEX = r"(?P<tag_prefix>.*)"
EXTENSION_REGEX = r"(?P<extension>h5ad)"


def get_identifier_type_and_value(identifier: str) -> Tuple[str, str]:
    identifier_type = None
    identifier_value = identifier
    if re.match(f"^{DATASET_ID_REGEX}$", identifier):
        # identifier is a uuid
        identifier_type = "dataset_uuid"
    else:
        # CURATOR_TAG_PREFIX_REGEX is superfluous; leaving in to match lambda handler code; may use later
        matched = re.match(f"({DATASET_ID_REGEX}|{CURATOR_TAG_PREFIX_REGEX})\\.{EXTENSION_REGEX}$", identifier)
        if matched:
            matches = matched.groupdict()
            if dataset_uuid := matches.get("dataset_uuid"):
                identifier_type = "dataset_uuid"
                identifier_value = dataset_uuid
            else:
                identifier_type = "curator_tag"

    if not identifier_type:
        raise Exception(f"The identifier '{identifier}' must either 1) include a '.h5ad' suffix OR 2) be a uuid")

    return identifier_value, identifier_type


def delete_dataset(collection_uuid: str, identifier: str):
    """
    Delete a private Dataset
    :param collection_uuid: the uuid of the Collection to which the Dataset belongs
    :param identifier: the curator tag or cellxgene Dataset uuid
    :return: True if deletion is successful otherwise False
    """
    url = url_builder(f"/collections/{collection_uuid}/datasets")
    headers = get_headers()

    identifier_value, identifier_type = get_identifier_type_and_value(identifier)

    params_dict = dict()
    params_dict[identifier_type] = identifier_value
    try:
        res = requests.delete(url, headers=headers, params=params_dict)
        res.raise_for_status()
    except Exception as e:
        logger.error("\n\033[1m\033[38;5;9mFAILED\033[0m")  # 'FAILED' in bold red
        raise e
    else:
        logger.info("\n\033[1m\033[38;5;10mSUCCESS\033[0m\n")  # 'SUCCESS' in bold green
        logger.info(f"Deleted the Dataset with {identifier_type} '{identifier_value}' from its Collection: "
                    f"\n{os.getenv('site_url')}/collections/{collection_uuid}")


def update_curator_tag(collection_uuid: str, identifier: str, new_tag: str):
    """
    Update a private Dataset's curator tag
    :param collection_uuid: the uuid of the Collection to which the Dataset belongs
    :param identifier: the curator tag or cellxgene Dataset uuid
    :param new_tag: the new curator tag to assign to the Dataset
    """
    url = url_builder(f"/collections/{collection_uuid}/datasets")
    headers = get_headers()

    identifier_value, identifier_type = get_identifier_type_and_value(identifier)

    params_dict = dict()
    params_dict[identifier_type] = identifier_value

    new_curator_tag_dict = dict(curator_tag=new_tag)
    try:
        res = requests.patch(url, headers=headers, params=params_dict, data=json.dumps(new_curator_tag_dict))
        res.raise_for_status()
    except Exception as e:
        logger.error("\n\033[1m\033[38;5;9mFAILED\033[0m")  # 'FAILED' in bold red
        raise e
    else:
        logger.info("\n\033[1m\033[38;5;10mSUCCESS\033[0m\n")  # 'SUCCESS' in bold green
        logger.info(f"Dataset with {identifier_type} '{identifier_value}' updated to have curator_tag '{new_tag}'")


def upload_local_datafile(datafile_path: str, collection_uuid: str, identifier: str):
    """
    :param datafile_path: the fully qualified path of the datafile to be uploaded
    :param collection_uuid: the uuid of the Collection to which the resultant Dataset will belong
    :param identifier: the curator tag or cellxgene Dataset uuid. Must be suffixed with '.h5ad'. See heading
    of upload_local_datafile.ipynb for details about how to use the identifier to 'create new' vs 'replace existing'
    :param log_level: the logging level
    Datasets.
    :return: None
    """
    url = url_builder(f"/collections/{collection_uuid}/datasets/s3-upload-credentials")
    headers = get_headers()

    def retrieve_s3_credentials_and_upload_key_prefix():
        return requests.post(url, headers=headers).json()

    log_level = os.getenv("log_level", "INFO")  # hack to determine log level in callback (separate process)
    time_zone_info = datetime.now(timezone.utc).astimezone().tzinfo

    def s3_refreshable_credentials_cb():
        res_data = retrieve_s3_credentials_and_upload_key_prefix()
        s3_credentials = res_data.get("Credentials")
        s3_credentials_formatted = {
            "access_key": s3_credentials.get("AccessKeyId"),
            "secret_key": s3_credentials.get("SecretAccessKey"),
            "token": s3_credentials.get("SessionToken"),
            "expiry_time": datetime.fromtimestamp(s3_credentials.get("Expiration")).replace(
                tzinfo=time_zone_info).isoformat(),
        }
        if getattr(logging, log_level) < 20:  # if log level NOTSET or DEBUG
            print("Retrieved/refreshed s3 credentials")
        return s3_credentials_formatted

    filesize = os.path.getsize(datafile_path)

    def get_progress_cb(collection_uuid: str, identifier: str):
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
                    print(f"{collection_uuid}/{identifier}: "
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
        logger.info(f"\nUploading {datafile_path} to Collection {collection_uuid} with identifier '{identifier}'...\n")
        s3.upload_file(
            Filename=datafile_path,
            Bucket=bucket,
            Key=upload_key,
            Callback=get_progress_cb(collection_uuid, identifier),
        )
    except Exception as e:
        logger.error(f"\n\033[1m\033[38;5;9mFAILED uploading:\033[0m {collection_uuid}/{identifier}")
        raise e
    else:
        logger.info("\n\033[1m\033[38;5;10mSUCCESS\033[0m\n")  # 'SUCCESS' in bold green
