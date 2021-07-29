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

To update these gene references, a developer has to:

### 1. Download gene reference files

#### 10X GTFs

The reference GTF files are obtained from the 10X reference data. To get the GTF files, from this directory (`cellxgene_schema_cli/scripts`) 
the following script has to be called:

```bash
bash gtf_download.sh
```

This script will:

1. Download the reference 10X data for both mouse and human.
2. Uncompress the data.
3. Copy the gtf file to to this directory.
4. Remove all the other data.

If GTF versions need to be updated then `gtf_download.sh` has to be changed accordingly.

### 2. Parse reference files

From this directory (`cellxgene_schema_cli/scripts`) the following script has to be called:

```bash
python3 ./gtf_processing.py
```

This will:

1. Parse the downloaded GTF file (in this directory).
2. Write the parsed versions on `cellxgene_schema_cli/cellxgene_schema/ontology_files/` into:

    - `gene_homo_sapiens.csv.gz`
    - `gene_mus_musculus.csv.gz`
