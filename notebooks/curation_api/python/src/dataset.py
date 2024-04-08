import logging
import os
import threading

import boto3
import requests
from botocore.credentials import RefreshableCredentials
from botocore.session import get_session
from src.utils.config import format_c_url
from src.utils.http import get_headers_and_cookies, url_builder
from src.utils.logger import failure, get_custom_logger, success

logger = get_custom_logger()


def create_dataset(collection_id: str):
    """
    Create a new, empty Dataset
    :param collection_id: the id of the Collection to which the Dataset belongs
    :return: None
    """
    url = url_builder(f"/collections/{collection_id}/datasets")

    try:
        res = requests.post(url, **get_headers_and_cookies())
        res.raise_for_status()
        data = res.json()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    dataset_id = data["dataset_id"]
    success(logger, f"Created new Dataset {dataset_id} in the Collection at {format_c_url(collection_id)}")
    return dataset_id


def delete_dataset(collection_id: str, dataset_id: str):
    """
    Delete a private Dataset
    :param collection_id: the id of the Collection to which the Dataset belongs
    :param dataset_id: Dataset id
    :return: True if deletion is successful otherwise False
    """
    url = url_builder(f"/collections/{collection_id}/datasets/{dataset_id}")

    success_message = (
        f"Deleted the Dataset with id '{dataset_id}' from its Collection: " f"\n{format_c_url(collection_id)}"
    )
    try:
        res = requests.delete(url, **get_headers_and_cookies())
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    else:
        success(logger, success_message)


def download_assets(datasets: list, log_version_id: bool = False):
    """
    Download assets locally for a list of Datasets
    :param log_version_id: bool to log dataset version ID rather than dataset ID
    :param datasets: list of full metadata Dataset json response objects.
    :return: None
    """
    try:
        for dataset in datasets:
            dataset_id = dataset["dataset_version_id"] if log_version_id else dataset["dataset_id"]
            assets = dataset["assets"]
            for asset in assets:
                download_filename = f"{dataset_id}.{asset['filetype']}"
                print(f"\nDownloading {download_filename}... ")
                with requests.get(asset["url"], stream=True) as res:
                    res.raise_for_status()
                    filesize = int(res.headers["Content-Length"])
                    with open(download_filename, "wb") as df:
                        total_bytes_received = 0
                        for chunk in res.iter_content(chunk_size=1024 * 1024):
                            df.write(chunk)
                            total_bytes_received += len(chunk)
                            percent_of_total_upload = float("{:.1f}".format(total_bytes_received / filesize * 100))
                            color = "\033[38;5;10m" if percent_of_total_upload == 100 else ""
                            print(f"\033[1m{color}{percent_of_total_upload}% downloaded\033[0m\r", end="")
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    success(logger, "Finished downloading assets")


def get_dataset(collection_id: str, dataset_id: str):
    """
    Get full metadata for a Dataset
    :param collection_id: the id of the Collection to which the Dataset belongs
    :param dataset_id: Dataset id
    :return: the full Dataset metadata
    """
    url = url_builder(f"/collections/{collection_id}/datasets/{dataset_id}")

    try:
        res = requests.get(url, **get_headers_and_cookies())
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    return res.json()


def get_dataset_version(dataset_version_id: str):
    """
    Get full metadata for a Dataset Version
    :param dataset_version_id: Dataset Version id
    :return: the full Dataset Version metadata
    """
    url = url_builder(f"/dataset_versions/{dataset_version_id}")

    try:
        res = requests.get(url, **get_headers_and_cookies())
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    return res.json()


def get_datasets(visibility: str = None, schema_version: str = None):
    """
    Get full metadata for all public Datasets
    """
    params = {}
    if visibility:
        params["visibility"] = visibility
    if schema_version:
        params["schema_version"] = schema_version
    url = url_builder("/datasets")

    try:
        res = requests.get(url, **get_headers_and_cookies(), params=params)
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    return res.json()


def get_dataset_versions(dataset_id: str):
    """
    Get list of metadata for all Published Versions of a Dataset
    :param dataset_id: Dataset id
    :return: the full Dataset metadata
    """
    url = url_builder(f"/datasets/{dataset_id}/versions")

    try:
        res = requests.get(url, **get_headers_and_cookies())
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    return res.json()


def upload_datafile_from_link(link: str, collection_id: str, dataset_id: str):
    """
    Create/update a Dataset from the datafile found at the source link.
    :param link: the source datafile link to upload to the Data Portal to become a Dataset
    :param collection_id: the id of the Collection to which the resultant Dataset will belong
    :param dataset_id: Dataset id.
    :return: None
    """
    url = url_builder(f"/collections/{collection_id}/datasets/{dataset_id}")

    data_dict = dict(link=link)

    success_message = (
        f"Uploading Dataset with id '{dataset_id}' to Collection "
        f"{os.getenv('SITE_URL')}/collections/{collection_id} sourcing from datafile at {link}"
    )

    try:
        res = requests.put(url, json=data_dict, **get_headers_and_cookies())
        res.raise_for_status()
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    else:
        success(logger, success_message)


def upload_local_datafile(datafile_path: str, collection_id: str, dataset_id: str):
    """
    :param datafile_path: the fully qualified path of the datafile to be uploaded
    :param collection_id: the id of the Collection to which the resultant Dataset will belong
    :param dataset_id: Dataset id.
    :param log_level: the logging level
    :return: None
    """
    url = url_builder(f"/collections/{collection_id}/s3-upload-credentials")

    def retrieve_s3_credentials_and_upload_key_prefix():
        try:
            res = requests.get(url, **get_headers_and_cookies())
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

    def get_progress_cb(collection_id: str, dataset_id: str):
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
                    print(
                        f"{collection_id}/{dataset_id}: " f"\033[1m{color}{percent_of_total_upload}% uploaded\033[0m\r",
                        end="",
                    )

        return progress_cb

    credentials_and_path = retrieve_s3_credentials_and_upload_key_prefix()
    bucket, key_prefix = credentials_and_path["Bucket"], credentials_and_path["UploadKeyPrefix"]
    upload_key = key_prefix + dataset_id
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
        logger.info(f"\nUploading {datafile_path} to Collection {collection_id} with dataset_id '{dataset_id}'...\n")
        s3.upload_file(
            Filename=datafile_path,
            Bucket=bucket,
            Key=upload_key,
            Callback=get_progress_cb(collection_id, dataset_id),
        )
    except requests.HTTPError as e:
        failure(logger, e)
        raise e
    else:
        success(logger, "UPLOAD COMPLETE -- Dataset is queued for processing and will surface in the system shortly.")
