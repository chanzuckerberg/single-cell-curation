import json
import os

import requests


def fetch_public_datasets(base_url):
    return requests.get(f"https://{base_url}/curation/v1/datasets").json()


def fetch_private_collections(base_url, headers):
    return requests.get(f"https://{base_url}/curation/v1/collections?visibility=PRIVATE", headers=headers).json()


def fetch_private_dataset(base_url, headers, collection_id, dataset_id):
    return requests.get(
        f"https://{base_url}/curation/v1/collections/{collection_id}/datasets/{dataset_id}", headers=headers
    ).json()


def get_headers(base_url):
    auth0_secrets = json.loads(os.getenv("AUTH0_SECRETS"))
    response = requests.post(
        f"https://{base_url}/curation/v1/auth/token", headers={"x-api-key": f"{auth0_secrets['super_curator_api_key']}"}
    )
    access_token = response.json()["access_token"]
    return {"Authorization": f"Bearer {access_token}"}


BASE_API = {
    "prod": "https://api.cellxgene.cziscience.com",
    "staging": "https://api.cellxgene.staging.single-cell.czi.technology",
    "dev": "https://api.cellxgene.dev.single-cell.czi.technology",
}
