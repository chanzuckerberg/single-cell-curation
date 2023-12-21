import json
import os

import boto3
import requests


def set_api_urls(env: str, stack: str = "") -> None:
    """
    This function sets url environment variables that other Curator API notebook modules use when calling
    Curator API endpoints.
    :param env: the deployment environment
    :param stack: the stack name (if rdev)
    :return: None
    """
    if stack:
        root_domain = "rdev.single-cell.czi.technology"
        os.environ["SITE_URL"] = f"https://{stack}-frontend.{root_domain}"
        os.environ["API_URL_BASE"] = f"https://{stack}-backend.{root_domain}"
        print(f"Stack: {stack}")
    else:
        if env == "prod":  # For all official use of the cellxgene product
            domain_name = "cellxgene.cziscience.com"
        elif env == "dev" or env == "staging":  # For testing purposes only
            domain_name = f"cellxgene.{env}.single-cell.czi.technology"
        else:
            raise Exception("Must provide env arg: 'dev', 'staging', or 'prod'.")
        os.environ["SITE_URL"] = f"https://{domain_name}"
        os.environ["API_URL_BASE"] = f"https://api.{domain_name}"

    print(f"Set 'SITE_URL' env var to {os.getenv('SITE_URL')}")
    print(f"Set 'API_URL_BASE' env var to {os.getenv('API_URL_BASE')}")


def set_api_access_config(
    api_key_file_path: str = None, env: str = "prod", stack: str = "", oauth_cookie: str = ""
) -> None:
    """
    This function uses the API key to retrieve a temporary access token from the Curator API. It then sets
    the 'ACCESS_TOKEN' environment variable, which other Curator API notebook modules use when calling
    Curator API endpoints. If no api_key_file_path arg is provided, then no access token will be set, and the user
    will be limited to only publicly-accessible data.
    :param api_key_file_path: the relative path to the file containing the API key
    :param env: the deployment environment
    :param stack: the stack name (if rdev)
    :param oauth_cookie: the user's _oauth2_proxy cookie associated with .rdev.single-cell.czi.technology
    :return: None
    """
    set_api_urls(env, stack)

    access_token_headers = {}
    access_token_cookies = {}
    if stack:
        if oauth_cookie:
            os.environ["OAUTH_COOKIE"] = oauth_cookie
            access_token_cookies["_oauth2_proxy"] = oauth_cookie
        else:
            # First generate a 'test app' oauth proxy token
            set_rdev_oauth_proxy_access_token()
            access_token_headers["Authorization"] = f"Bearer {os.getenv('OAUTH_PROXY_ACCESS_TOKEN')}"

    if not api_key_file_path:
        print(
            "No API key file provided. Without an access token, which requires an API key to retrieve, no private "
            "data or actions are accessible (except for fetching an individual Collection or Dataset, which only "
            "requires the url)."
        )
        return

    with open(api_key_file_path) as fp:
        api_key = fp.read().strip()
        access_token_headers["x-api-key"] = api_key

    access_token_path = "/curation/v1/auth/token"
    api_url_base = os.getenv("API_URL_BASE")
    access_token_url = f"{api_url_base}{access_token_path}"

    res = requests.post(access_token_url, headers=access_token_headers, cookies=access_token_cookies)
    res.raise_for_status()
    access_token = res.json().get("access_token")
    os.environ["ACCESS_TOKEN"] = access_token
    print("Successfully set 'ACCESS_TOKEN' env var!")


def set_rdev_oauth_proxy_access_token() -> None:
    """
    Get an access token for the outer Auth0 layer that gates all access to single-cell rdev
    """
    secrets_client = boto3.client("secretsmanager")
    secret = json.loads(secrets_client.get_secret_value(SecretId="corpora/backend/rdev/auth0-secret")["SecretString"])
    payload = {
        "client_id": secret["test_app_id"],
        "client_secret": secret["test_app_secret"],
        "grant_type": "client_credentials",
        "audience": "https://api.cellxgene.dev.single-cell.czi.technology/dp/v1/curator",
    }
    headers = {"content-type": "application/json"}
    res = requests.post("https://czi-cellxgene-dev.us.auth0.com/oauth/token", json=payload, headers=headers)
    res.raise_for_status()
    os.environ["OAUTH_PROXY_ACCESS_TOKEN"] = res.json()["access_token"]
    print("Successfully set 'OAUTH_PROXY_ACCESS_TOKEN' env var!")


def format_c_url(collection_id: str) -> str:
    return f"{os.getenv('SITE_URL')}/collections/{collection_id}"
