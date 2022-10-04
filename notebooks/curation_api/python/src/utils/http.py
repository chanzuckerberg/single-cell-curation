import os

from src.utils.logger import get_custom_logger


logger = get_custom_logger()


def url_builder(path_segment):
    api_url_base = os.getenv("api_url_base")
    if not api_url_base:
        raise Exception("The required 'api_url_base' env var is not set. Call set_api_urls() from src.utils.config")
    route_path = f"/curation/v1{path_segment}"
    logger.debug(f"route path: {route_path}")
    url = f"{api_url_base}{route_path}"
    logger.debug(f"url: {url}")
    return url


def get_headers():
    access_token = os.getenv("access_token")
    headers = {"Content-Type": "application/json"}
    if not access_token:
        logger.warning("No access token included in request")
    else:
        headers["Authorization"] = f"Bearer {access_token}"
    return headers
