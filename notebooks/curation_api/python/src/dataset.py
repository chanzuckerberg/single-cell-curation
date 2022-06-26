import boto3
import logging
import os
import requests
import threading
from typing import Union
from botocore.credentials import RefreshableCredentials
from botocore.session import get_session
from datetime import datetime, timezone

from src.utils.logger import get_custom_logger


def upload_local_datafile(
        datafile_path: str,
        collection_uuid: str,
        identifier: str,
        log_level: Union[int, str] = "INFO",
) -> bool:
    """
    :param datafile_path: the fully qualified path of the datafile to be uploaded
    :param collection_uuid: the uuid of the Collection to which the resultant Dataset will belong
    :param identifier: the curator tag or cellxgene Dataset uuid. Must be suffixed with '.h5ad'. See heading
    of upload_local_datafile.ipynb for details about how to use the identifier to 'create new' vs 'replace existing'
    :param log_level: the logging level
    Datasets.
    :return: True if upload succeeds otherwise False
    """
    logger = get_custom_logger(log_level)

    s3_credentials_path = f"/curation/v1/collections/{collection_uuid}/datasets/s3-upload-credentials"
    s3_credentials_url = f"{os.getenv('api_url_base')}{s3_credentials_path}"
    s3_cred_headers = {"Authorization": f"Bearer {os.getenv('access_token')}"}

    def retrieve_s3_credentials_and_upload_key_prefix():
        return requests.post(s3_credentials_url, headers=s3_cred_headers).json()

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
        if getattr(logging, log_level) < 20:
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
        logger.error(f"\n{e}")
        return False
    return True
