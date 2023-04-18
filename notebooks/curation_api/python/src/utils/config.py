import os

import requests


def set_api_urls(env: str) -> None:
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

    os.environ["SITE_URL"] = f"https://{domain_name}"
    print(f"Set 'SITE_URL' env var to {os.getenv('SITE_URL')}")
    os.environ["API_URL_BASE"] = f"https://api.{domain_name}"
    print(f"Set 'API_URL_BASE' env var to {os.getenv('API_URL_BASE')}")


def set_api_access_config(api_key_file_path: str = None, env: str = "prod") -> None:
    """
    This function uses the API key to retrieve a temporary access token from the Curator API. It then sets
    the 'ACCESS_TOKEN' environment variable, which other Curator API notebook modules use when calling
    Curator API endpoints. If no api_key_file_path arg is provided, then no access token will be set, and the user
    will be limited to only publicly-accessible data.
    :param api_key_file_path: the relative path to the file containing the API key
    :param env: the deployment environment
    :return: None
    """
    set_api_urls(env)
    if not api_key_file_path:
        print(
            "No API key file provided. Without an access token, which requires an API key to retrieve, no private "
            "data or actions are accessible (except for fetching an individual Collection or Dataset, which only "
            "requires the url)."
        )
        return

    with open(api_key_file_path) as fp:
        api_key = fp.read().strip()
    access_token_headers = {"x-api-key": api_key}
    access_token_path = "/curation/v1/auth/token"
    api_url_base = os.getenv("API_URL_BASE")
    access_token_url = f"{api_url_base}{access_token_path}"

    res = requests.post(access_token_url, headers=access_token_headers)
    res.raise_for_status()
    access_token = res.json().get("access_token")
    os.environ["ACCESS_TOKEN"] = access_token
    print("Successfully set 'ACCESS_TOKEN' env var!")


def format_c_url(collection_id: str) -> str:
    return f"{os.getenv('SITE_URL')}/collections/{collection_id}"
