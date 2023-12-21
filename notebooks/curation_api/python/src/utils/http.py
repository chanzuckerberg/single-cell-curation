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


def get_headers_and_cookies() -> dict:
    headers = {"Content-Type": "application/json"}
    cookies = {}
    if access_token := os.getenv("ACCESS_TOKEN"):
        # Discover API access token will also allow request to get through rdev proxy
        headers["Authorization"] = f"Bearer {access_token}"
    elif oauth_cookie := os.getenv("OAUTH_COOKIE"):
        # Only use oauth cookie if access_token not set -- does not confer any Discover API permissions
        cookies["_oauth2_proxy"] = oauth_cookie
    elif oauth_proxy_token := os.getenv("OAUTH_PROXY_ACCESS_TOKEN"):
        # Or use generic 'test app' token -- does not confer any Discover API permissions
        headers["Authorization"] = f"Bearer {oauth_proxy_token}"
    else:
        logger.warning("No access token included in request")
    return {"headers": headers, "cookies": cookies}
