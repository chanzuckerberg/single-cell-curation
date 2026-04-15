#!/usr/bin/env python3
"""Submit guidescan2 index Argo jobs via kubectl for all species in gene_assembly.yml.

Reads gene_assembly.yml, submits a Workflow manifest referencing the
guidescan2-index-v1 ClusterWorkflowTemplate for each species, and polls
kubectl until every workflow reaches a terminal phase (Succeeded/Failed/Error).

This script is invoked by the guidescan-index.yml GHA workflow whenever
gene_assembly.yml changes. It always processes every species in the file;
the GHA path trigger is the change-detection mechanism.

kubectl must be on PATH and configured to reach the staging cluster
(e.g. via `aws eks update-kubeconfig`). No Argo CLI required.

Exit codes:
  0 — all workflows Succeeded
  1 — one or more workflows Failed/Errored, timed out, or config missing
"""

import json
import subprocess
import sys
import time
from pathlib import Path

import yaml

GENE_ASSEMBLY_PATH = Path(__file__).parents[1] / "cellxgene_schema" / "gencode_files" / "gene_assembly.yml"

ARGO_NAMESPACE = "argo-workflows"
WORKFLOW_TEMPLATE = "guidescan2-index-v1"
ARGO_SERVICE_ACCOUNT = "sci-data-staging-guidescan2-sa"
S3_OUTPUT_PATH = "s3://czi-biohub-references/guidescan-indexes-test/"

TERMINAL_PHASES = {"Succeeded", "Failed", "Error"}
POLL_INTERVAL_SECONDS = 30
POLL_TIMEOUT_SECONDS = 14400  # 4 hours — FASTA download + indexing is slow


def build_workflow_manifest(species_key: str, fasta_url: str) -> dict:
    return {
        "apiVersion": "argoproj.io/v1alpha1",
        "kind": "Workflow",
        "metadata": {
            "generateName": f"guidescan2-index-{species_key}-",
            "namespace": ARGO_NAMESPACE,
        },
        "spec": {
            "workflowTemplateRef": {
                "name": WORKFLOW_TEMPLATE,
                "clusterScope": True,
            },
            "arguments": {
                "parameters": [
                    {"name": "fasta_url", "value": fasta_url},
                    {"name": "output_path", "value": S3_OUTPUT_PATH},
                ]
            },
            "serviceAccountName": ARGO_SERVICE_ACCOUNT,
        },
    }


def submit_workflow(species_key: str, fasta_url: str) -> str:
    """Submit the workflow manifest via kubectl and return the created workflow name."""
    manifest_yaml = yaml.dump(
        build_workflow_manifest(species_key, fasta_url),
        default_flow_style=False,
    )

    result = subprocess.run(
        ["kubectl", "create", "-f", "-", "-n", ARGO_NAMESPACE, "-o", "json"],
        input=manifest_yaml,
        capture_output=True,
        text=True,
        check=True,
    )

    workflow_name: str = json.loads(result.stdout)["metadata"]["name"]
    print(f"[{species_key}] Submitted workflow: {workflow_name}")
    return workflow_name


def wait_for_workflow(species_key: str, workflow_name: str) -> str:
    """Poll kubectl until the workflow reaches a terminal phase; return that phase."""
    elapsed = 0
    while elapsed < POLL_TIMEOUT_SECONDS:
        result = subprocess.run(
            [
                "kubectl",
                "get",
                "workflow",
                workflow_name,
                "-n",
                ARGO_NAMESPACE,
                "-o",
                "jsonpath={.status.phase}",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        phase = result.stdout.strip()
        print(f"[{species_key}] {workflow_name}: {phase or 'Pending'} ({elapsed}s elapsed)")

        if phase in TERMINAL_PHASES:
            return phase

        time.sleep(POLL_INTERVAL_SECONDS)
        elapsed += POLL_INTERVAL_SECONDS

    raise TimeoutError(f"[{species_key}] Timed out after {POLL_TIMEOUT_SECONDS}s waiting for {workflow_name}")


def submit_and_wait(species_key: str, fasta_url: str) -> None:
    print(f"\n[{species_key}] fasta_url   : {fasta_url}")
    print(f"[{species_key}] output_path : {S3_OUTPUT_PATH}")
    print(f"[{species_key}] template    : ClusterWorkflowTemplate/{WORKFLOW_TEMPLATE}")

    workflow_name = submit_workflow(species_key, fasta_url)
    phase = wait_for_workflow(species_key, workflow_name)

    if phase != "Succeeded":
        raise RuntimeError(f"[{species_key}] Workflow {workflow_name} ended with phase: {phase}")
    print(f"[{species_key}] Workflow {workflow_name} Succeeded")


def main() -> None:
    with open(GENE_ASSEMBLY_PATH) as f:
        gene_assembly = yaml.safe_load(f) or {}

    errors: list[str] = []

    for species_key, info in gene_assembly.items():
        fasta_url = info.get("fasta_url")
        if not fasta_url:
            print(f"ERROR: no fasta_url configured for '{species_key}'", file=sys.stderr)
            errors.append(species_key)
            continue

        try:
            submit_and_wait(species_key, fasta_url)
        except (subprocess.CalledProcessError, RuntimeError, TimeoutError) as e:
            print(f"ERROR: {e}", file=sys.stderr)
            errors.append(species_key)

    if errors:
        print(f"\nFailed: {', '.join(errors)}", file=sys.stderr)
        sys.exit(1)

    print("\nAll jobs completed successfully.")


if __name__ == "__main__":
    main()
