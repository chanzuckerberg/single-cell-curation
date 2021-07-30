## Updating ontologies

Ontology information for cellxgene-shcema is  stored in `cellxgene_schema_cli/cellxgene_schema/ontology_files/all_ontology.json.gz` 

To update these ontologies, a developer has to:

### 1. Modify owl references

The yaml file `cellxgene_schema_cli/cellxgene_schema/ontology_files/owl_info` contains the current ontologies and versions used.

Foe example, for CL:

```yaml
CL:
  latest: 2021-06-21
  urls:
    2021-06-21: https://github.com/obophenotype/cell-ontology/raw/v2021-06-21/cl.owl
```

`latest` cotains the current version used in cellxgene-schema. When updating an ontology, a new version has to be added to `urls` and 
`latest` has to be modified accordingly. It would like look this:

```yaml
CL:
  latest: new_date
  urls:
    2021-06-21: https://github.com/obophenotype/cell-ontology/raw/v2021-06-21/cl.owl
    new_date: https://foo.owl
```

**Note:** The processing scripts currently support `.owl` and `.owl.gz` files

### 2. Download and parse ontologies

From this directory (`cellxgene_schema_cli/scripts`) the following script has to be called:

```bash
python3 ./ontology_processing.py
```

This will:

1. Download ontology owl files to `cellxgene_schema_cli/cellxgene_schema/ontology_files/`
2. Parse each owl, extract important information and store all of it into **one** file `cellxgene_schema_cli/cellxgene_schema/ontology_files/all_ontology.json`

## Updating gene references

Gene references are stored under  `cellxgene_schema_cli/cellxgene_schema/ontology_files/`

The following gene references are currently stored (to see specific versions refer to schema definition):

1. Humnan `genes_homo_sapines.csv.gz`
2. Mouse `genes_mus_musculus_sapines.csv.gz`
3. sars‑CoV‑2 `genes_sars_cov_2.csv.gz`
4. ERCC Spike-Ins `genes_ercc.csv.gz`

To update these gene references, a developer has to:

### 1. Download gene reference files

The reference GTF files for mouse and human are obtained from the 10X reference data. The GTF file for sars‑CoV‑2 is obtained from ENSEMBL. The ERCC IDs are obtained from ThermoFisher.
To download all these files, from this directory (`cellxgene_schema_cli/scripts`) the following script has to be called:

```bash
bash gene_download.sh
```

This script will:

1. Download the reference gene files.
2. For 10X data, uncompress it, get the GTF files and remove the rest.
3. At the end the following files will be in this folder:

    - `mus_musculus.gtf.gz`
    - `homo_sapiens.gtf.gz`
    - `sars_cov_2.gtf.gz`
    - `ercc.txt`

If reference gene versions need to be updated then `gtf_download.sh` has to be changed accordingly.

### 2. Parse reference files

From this directory (`cellxgene_schema_cli/scripts`) the following script has to be called:

```bash
python3 ./gene_processing.py
```

This will:

1. Parse the downloaded reference files 
2. Write the parsed versions on `cellxgene_schema_cli/cellxgene_schema/ontology_files/` into:

    - `genes_homo_sapiens.csv.gz`
    - `genes_mus_musculus.csv.gz`
    - `genes_sars_cov_2.csv.gz`
    - `genes_ercc.csv.gz`
