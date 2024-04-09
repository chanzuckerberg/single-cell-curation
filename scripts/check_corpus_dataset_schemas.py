import json
import os
from collections import defaultdict
from time import sleep
from typing import Any, Dict, List

import requests

"""
This script leverages the Discover API to give a breakdown of the Datasets of each schema version that is 
represented by

1) Public Datasets, or
2) Private and Revision Datasets


Requirements:

- The 'visibility' variable must be set to "PUBLIC" for 1) and to "PRIVATE" for 2)

1) produces a public.json output
2) produces a private.json and revision.json output

- The 'api_url' variable for the Discover API must be set appropriately for the chosen environment


Invocation:

[ACCESS_TOKEN=<access_token>] python check_corpus_dataset_schemas.py


Results:

In both cases, json structure:

{
  <schema_version>: [
    [<collection_id>, <collection_version_id>, <dataset_id>, <dataset_version_id>],
    ...
  ],
  <schema_version>: [
    ...
  ],
  ...
}

Example:
{
  "5.0.0": [
    [00000000-0000-0000-0000-000000000000, 11111111-1111-1111-1111-111111111111, 22222222-2222-2222-2222-222222222222, 33333333-3333-3333-3333-333333333333]
  ]
}
"""


api_url = "https://api.cellxgene.dev.single-cell.czi.technology/curation/v1"  # dev
# api_url = "https://api.cellxgene.cziscience.com/curation/v1"  # prod

visibility = "PRIVATE"

if visibility == "PUBLIC":
    collections = requests.get(f"{api_url}/collections?visibility={visibility}").json()  # For public Collections
    public_datasets: Dict[str, List[Any]] = defaultdict(list)
    for c in collections:
        sleep(1)
        resp = requests.get(f"{api_url}/collections/{c['collection_id']}")
        if resp.status_code != 200:
            print(resp)
            print(f"Unable to retrieve collection {c['collection_id']}")
            continue
        datasets = resp.json()["datasets"]
        for d in datasets:
            schema = d["schema_version"]
            print(
                f"adding dataset version {d['dataset_version_id']} from collection version {c['collection_version_id']} to public. Count {schema}: {len(public_datasets[schema])}"
            )
            public_datasets[schema].append(
                (c["collection_id"], c["collection_version_id"], d["dataset_id"], d["dataset_version_id"])
            )
    with open("public.json", "w") as fp:
        json.dump(public_datasets, fp)

elif visibility == "PRIVATE":
    access_token = os.getenv("ACCESS_TOKEN")
    assert access_token and len(access_token) > 0

    collections = requests.get(
        f"{api_url}/collections?visibility={visibility}", headers={"Authorization": f"Bearer {access_token}"}
    ).json()  # For private Collections

    revision_datasets: Dict[str, List[Any]] = defaultdict(list)
    private_datasets: Dict[str, List[Any]] = defaultdict(list)

    no_schema_version = []
    no_schema_private = 0
    no_schema_revisions = 0
    for c in collections:
        sleep(1)
        if c["revision_of"]:
            group = revision_datasets
            name = "revisions"
        else:
            group = private_datasets
            name = "private"
        resp = requests.get(f"{api_url}/collections/{c['collection_id']}")
        if resp.status_code != 200:
            print(f"ERROR: {resp.status_code}")
            print(c["collection_id"])
        datasets = resp.json()["datasets"]
        for d in datasets:
            if "schema_version" not in d:
                if c["revision_of"]:
                    print(
                        f"no schema version for dataset {d['dataset_id']} {d['dataset_version_id']} in REVISION collection {c['collection_id']}"
                    )
                    no_schema_version.append(
                        (c["collection_id"], c["collection_version_id"], d["dataset_id"], d["dataset_version_id"])
                    )
                    no_schema_revisions += 1
                else:
                    print(
                        f"no schema version for dataset {d['dataset_id']} {d['dataset_version_id']} in PRIVATE collection {c['collection_id']}"
                    )
                    no_schema_private += 1
                continue
            schema = d["schema_version"]
            print(
                f"adding dataset version {d['dataset_version_id']} from collection version {c['collection_version_id']} to {name}. Count {schema}: {len(group[schema])}"
            )
            group[schema].append(
                (c["collection_id"], c["collection_version_id"], d["dataset_id"], d["dataset_version_id"])
            )

    print(f"Count of no_schema_version: {len(no_schema_version)}")
    print(f"Count of no_schema_revisions: {no_schema_revisions}")
    print(f"Count of no_schema_private: {no_schema_private}")

    # Write results to output files
    with open("no_schema_version.json", "w") as fp:
        json.dump(
            no_schema_version, fp
        )  # All datasets that do not have a schema version. Indicates processing incomplete / error. Array of tuples where each tuple is a dataset: (c_id, c_v_id, d_id, d_v_id).
    with open("private.json", "w") as fp:
        json.dump(private_datasets, fp)
    with open("revisions.json", "w") as fp:
        json.dump(revision_datasets, fp)
