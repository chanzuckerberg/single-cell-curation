import os
import requests


def set_api_urls(env: str = "prod"):
    """
    This function sets url environment variables that other Curator API notebook modules use when calling
    Curator API endpoints.
    :param env: the deployment environment
    :return: None
    """
    if env == "prod":  # For all official use of the cellxgene product
        domain_name = "cellxgene.cziscience.com"
    elif env == "dev" or env == "staging":  # For testing purposes only
        domain_name = f"cellxgene.{env}.single-cell.czi.technology"
    else:
        raise Exception("Must provide env arg: 'dev', 'staging', or 'prod'.")

    os.environ["site_url"] = f"https://{domain_name}"
    print(f"Set 'site_url' env var to {os.getenv('site_url')}")
    os.environ["api_url_base"] = f"https://api.{domain_name}"
    print(f"Set 'api_url_base' env var to {os.getenv('api_url_base')}")


def set_access_token(api_key_file_path: str):
    """
    This function uses the API key to retrieve a temporary access token from the Curator API. It then sets
    the 'access_token' environment variable, which other Curator API notebook modules use when calling
    Curator API endpoints. The usage of this function relies on a *prior* set_api_urls() call (from api_urls.py).
    :param api_key_file_path: the relative path to the file containing the API key
    :param env: the deployment environment
    :return: None
    """
    api_key = open(api_key_file_path).read().strip()
    access_token_headers = {"x-api-key": api_key}
    access_token_path = "/curation/v1/auth/token"

    api_url_base = os.getenv("api_url_base")
    if not api_url_base:
        raise Exception("You must first set the API url environment variables with set_api_urls()")

    access_token_url = f"{api_url_base}{access_token_path}"
    res = requests.post(access_token_url, headers=access_token_headers)
    res.raise_for_status()
    access_token = res.json().get("access_token")
    os.environ["access_token"] = access_token
    print("Successfully set 'access_token' env var!")
