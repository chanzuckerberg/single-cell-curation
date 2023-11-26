import os

from src.utils.logger import get_custom_logger

logger = get_custom_logger()


def url_builder(path_segment):
    api_url_base = os.getenv("API_URL_BASE")
    if not api_url_base:
        raise Exception("The required 'API_URL_BASE' env var is not set. Call set_api_urls() from src.utils.config")
    route_path = f"/curation/v1{path_segment}"
    logger.debug(f"route path: {route_path}")
    url = f"{api_url_base}{route_path}"
    logger.debug(f"url: {url}")
    return url


def get_headers():
    access_token = os.getenv("ACCESS_TOKEN")
    oauth_proxy_token = os.getenv("OAUTH_PROXY_ACCESS_TOKEN")
    headers = {"Content-Type": "application/json"}
    if oauth_proxy_token:
        # rdev Auth0 uses Authorization header for proxy token
        headers["Authorization"] = f"Bearer {oauth_proxy_token}"
        if access_token:
            # Use alternative header for Discover API access token
            headers["X-Curation-Authorization"] = access_token
        else:
            logger.warning("No access token included in request")
    elif access_token:
        headers["Authorization"] = f"Bearer {access_token}"
    else:
        logger.warning("No access token included in request")
    return headers
