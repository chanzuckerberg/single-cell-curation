import json
import os

import requests


def fetch_public_datasets(base_url):  # type: ignore
    return requests.get(f"https://{base_url}/curation/v1/datasets").json()


def fetch_private_collections(base_url, headers):  # type: ignore
    return requests.get(f"https://{base_url}/curation/v1/collections?visibility=PRIVATE", headers=headers).json()


def fetch_private_dataset(base_url, headers, collection_id, dataset_id):  # type: ignore
    return requests.get(
        f"https://{base_url}/curation/v1/collections/{collection_id}/datasets/{dataset_id}", headers=headers
    ).json()


def get_headers(base_url):  # type: ignore
    auth0_secrets = json.loads(os.getenv("AUTH0_SECRETS"))  # type: ignore
    response = requests.post(
        f"https://{base_url}/curation/v1/auth/token", headers={"x-api-key": f"{auth0_secrets['super_curator_api_key']}"}
    )
    access_token = response.json()["access_token"]
    return {"Authorization": f"Bearer {access_token}"}


BASE_API = {
    "prod": "api.cellxgene.cziscience.com",
    "staging": "api.cellxgene.staging.single-cell.czi.technology",
    "dev": "api.cellxgene.dev.single-cell.czi.technology",
}
