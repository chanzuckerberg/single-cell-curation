# Shendure data curation

## Project details

Last updated: Wed Nov 18 19:08:34 2020

Project title: Survey of human embryonic development

DOI: [http://dx.doi.org/10.1126/science.aba7721](http://dx.doi.org/10.1126/science.aba7721)

Download original (mod by Marcus). S3 buckets are under the `single-cell-dev` role at CZI's AWS account:

- Full data 

`s3://pablo-tmp/sc_datasets/shendure/Survey_of_human_embryonic_development-processed.h5ad`

- Down-sampled 1M cells 

`s3://pablo-tmp/sc_datasets/shendure/subsampled.h5ad`

Download after remixing

`s3://pablo-tmp/sc_datasets/shendure/Survey_of_human_embryonic_development-processed_remixed.h5ad`
`s3://pablo-tmp/sc_datasets/shendure/subsampled_remixed.h5ad`

Reqs:
- scanpy
- cellxgene with schema functions
- R with yaml package


Notes:

- Original data has unique identifiers, no need to worry about that anymore
- Down-sampled data is for cellxgene pruposes -- it gets overwhelmed by full dataset
- From Ambrose "We should remix the full dataset and view the down sampling as a modification for cellxgene visualization. The computational community, and therefore the data portal will benefit from the full dataset"

## Curation

### Quick-glance pipeline

- Download original h5ad files from the sources shown above. They should be stored in the same path as is this repo

- Create a table of suggested ontology terms for cell types and tissues

```bash
curl https://raw.githubusercontent.com/chanzuckerberg/cellxgene/78176f971159804da93a01aa60012a34d670bd62/server/converters/schema/ontology.py > scripts/ontology.py 

python3 scripts/create_ontology_lookup_cell_type.py > scripts/ontology_lookup_cell_type.tsv
python3 scripts/create_ontology_lookup_tissue.py > scripts/ontology_lookup_tissue.tsv
python3 scripts/create_ontology_lookup_dev_stage.py > scripts/ontology_lookup_dev_stage_curated.tsv

```

- Manually curate and select appropriate ontology terms for each item in the previous tables, do so by changing `True` to `False` in column 5.
Currently the tables are [here](https://docs.google.com/spreadsheets/d/14NfyiUWGOzgcRg7nUDiLNO9rMHUef6-gM1f7rlR5ymE/edit?usp=sharing). 

Save curated tables to:

```bash
scripts/ontology_lookup_cell_type_curated.tsv
scripts/ontology_lookup_tissue_curated.tsv
scripts/ontology_lookup_dev_stage_curated.tsv
```

- Append disease status (there's one sample with trisomy 18):

```bash
python3 scripts/append_dissease_state.py subsampled.h5ad subsampled_w_disease.h5ad
python3 scripts/append_dissease_state.py Survey_of_human_embryonic_development-processed.h5ad Survey_of_human_embryonic_development-processed_w_disease.h5ad
```


- Create config yaml for remixing with 

```bash
Rscript scripts/create_yaml.R scripts/ontology_lookup_cell_type_curated.tsv scripts/ontology_lookup_tissue_curated.tsv scripts/ontology_lookup_dev_stage_curated.tsv schema-shendure.yml
```

- Remix and verify datasets

```bash
cellxgene schema apply --source-h5ad subsampled_w_disease.h5ad --remix-config schema-shendure.yml --output-filename subsampled_remixed.h5ad
cellxgene schema apply --source-h5ad Survey_of_human_embryonic_development-processed_w_disease.h5ad --remix-config schema-shendure.yml --output-filename Survey_of_human_embryonic_development-processed_remixed.h5ad

cellxgene schema validate subsampled_remixed.h5ad
cellxgene schema validate Survey_of_human_embryonic_development-processed_remixed.h5ad
```

### Curation process

#### Cell type ontologies

Luckily the field `sub_cluster_name` in the `obs` data frame from the anndata object has very clear cell type definitions, and this seems to be the most granular cell types among all other fields that contain some sort of cell type information.

These arscriion of cell types

**1. Create a table of suggested ontology terms from the ontology look-up service**

This uses the look-up function that Marcus created [here]("https://github.com/chanzuckerberg/cellxgene/blob/78176f971159804da93a01aa60012a34d670bd62/server/converters/schema/ontology.py#L65")

Fist download the script containing that function
```bash
wget -O scripts/ontology.py https://raw.githubusercontent.com/chanzuckerberg/cellxgene/78176f971159804da93a01aa60012a34d670bd62/server/converters/schema/ontology.py
```

Then run custom script to create the table
```
python3 scripts/create_ontology_lookup_cell_type.py > scripts/ontology_lookup_cell_type.tsv
```

This creates a table with five columns:
```
original_cell_type	stripped_cell_type	ontology_term_id	ontology_term_name	final
Eye-Retinal progenitors and Muller glia-1	Retinal progenitors and Muller glia	CL:0000511	androgen binding protein secreting cell	False
Eye-Retinal progenitors and Muller glia-1	Retinal progenitors and Muller glia	CL:0000593	androgen secreting cell	False
Eye-Retinal progenitors and Muller glia-1	Retinal progenitors and Muller glia	CL:0000519	phagocyte (sensu Nematoda and Protostomia)	False
Eye-Retinal progenitors and Muller glia-1	Retinal progenitors and Muller glia	CL:0000338	neuroblast (sensu Nematoda and Protostomia)	False
Eye-Retinal progenitors and Muller glia-1	Retinal progenitors and Muller glia	CL:0000340	glioblast (sensu Nematoda and Protostomia)	False
Eye-Retinal progenitors and Muller glia-1	Retinal progenitors and Muller glia	CL:0000385	prohemocyte (sensu Nematoda and Protostomia)	False
Eye-Retinal progenitors and Muller glia-1	Retinal progenitors and Muller glia	CL:0000468	neuroglioblast (sensu Nematoda and Protostomia)	False
Eye-Retinal progenitors and Muller glia-1	Retinal progenitors and Muller glia	CL:0011107	Muller cell	False
Eye-Retinal progenitors and Muller glia-1	Retinal progenitors and Muller glia	CL:0000568	amine precursor uptake and decarboxProtostom
```

Each unique cell type in the original data is assigned to 0 or more ontology results (rows). The columns are:

- The first two columns are from the original dataset (`original_cell_type` is the cell type as found originally, `stripped_cell_type` is the version use in the search
- `ontology_term_id` and `ontology_term_name` are the results from the search)
- `final` is a boolean, if True this is the ontology mapping used for remixing

**2. Manually select best ontology match for each cell type**

I uploaded the table to a [google sheet](https://docs.google.com/spreadsheets/d/14NfyiUWGOzgcRg7nUDiLNO9rMHUef6-gM1f7rlR5ymE/edit?usp=sharing) for collaboration. The last column `final` should be manually changed to `True` for the best match of each cell type. If there was only one result for a given cell type, this should already be set as 'True'. If no results were found then we should do a manual search in the ontology look-up service and modify the table accordingly.

There are some cell types that:

- Didn't exist in the ontology look up service
- Were difficult to assign to an ontology term because the level of granularity of the definition didn't fully match the original cell type

The final version of the table needs to be stored here

```bash
scripts/ontology_lookup_cell_type_curated.tsv
```

#### Tissue ontologies

The field `Organ` in the `obs` data frame from the anndata object has very clear tissue definitions.

These are the steps to follow for the curation of tissues.

**1. Create a table of suggested ontology terms from the ontology look-up service**

Run custom script to create the table

```
python3 scripts/create_ontology_lookup_tissue.py > scripts/ontology_lookup_tissue.tsv
```

This cretes a table with five colums:
```
original_tisue	stripped_tissue	ontology_term_id	ontology_term_name	final
Eye	eye	UBERON:0000970	eye	True
Eye	eye	UBERON:0001702	eyelash	False
Eye	eye	UBERON:0010230	eyeball of camera-type eye	False
Eye	eye	UBERON:0004859	eye gland	False
Eye	eye	UBERON:0001711	eyelid	False
```

Each tissue in the original data is assigned >1 ontology results (rows). The columns are:

- The first two columns are from the original ddefinitionataset (`original_tissue` is the cell type as found originally, `stripped_tissue` is the version use in the search
- `ontology_term_id` and `ontology_term_name` are the results from the search
- `final` is a boolean, if True this is the ontology mapping used for remixing

**2. Manually select best ontology match for tissue**

I uploaded the table to a [google sheet](https://docs.google.com/spreadsheets/d/14NfyiUWGOzgcRg7nUDiLNO9rMHUef6-gM1f7rlR5ymE/edit?usp=sharing) for collaboration. The last column `final` should be manually changed to `True` for the best match of each cell type.

The final version of the table needs to be stored here

```bash
scripts/ontology_lookup_cell_type_curated.tsv
```

#### Developmental stage ontologies

The field `Development_day` in the `obs` data frame from the anndata object has the number of days post-fertilization. 

Run custom script to create the table. This create exactly one mapping per original entry

```
python3 scripts/create_ontology_lookup_dev_stage.py > ontology_lookup_dev_stage_curated.tsv
```

I added these mappings to the same [google sheet](https://docs.google.com/spreadsheets/d/14NfyiUWGOzgcRg7nUDiLNO9rMHUef6-gM1f7rlR5ymE/edit?usp=sharing) mentioned in the previous sections.

#### Disease state

There's one sample of cerebrum that comes from a donor with trisomy 18, this information needs to be added to the datasets. I made a custom python script to do so:

```bash
python3 scripts/append_dissease_state.py subsampled.h5ad subsampled_w_disease.h5ad
python3 scripts/append_dissease_state.py Survey_of_human_embryonic_development-processed.h5ad Survey_of_human_embryonic_development-processed_w_disease.h5ad
```

This will create two files with the `_w_disease` suffixes that are then use for remixing.

#### Remixing datasets

These are the steps to remix the data

**1. Create config yaml file with ontology mappings**

I made a custom R script that takes the curated ontology tables from the previous section and creates a yaml config file with the mappings ready for `cellxgene schema apply`

```bash
Rscript scripts/create_yaml.R scripts/ontology_lookup_cell_type.tsv scripts/ontology_lookup_tissue.tsv scripts/ontology_lookup_dev_stage_curated.tsv schema-shendure.yml
```

This will create the `schema-shendure.yml` with the mappings and along these fields:

Obs:

- assay_ontology_term_id: "EFO:0010550"
- disease_ontology_term_id: "PATO:0000461"
- ethnicity_ontology_term_id: "unknown"

Uns: 

- organism: "Homo sapiens"
- organism_ontology_term_id: "NCBITaxon:9606"
- layer_descriptions 
- publication_doi: "http://dx.doi.org/10.1126/science.aba7721"
- title: "Survey of human embryonic development"


**2. Remix datasets**

Next I applied the mapping and injected the metadata in the datasets

```bash
cellxgene schema apply --source-h5ad subsampled_w_disease.h5ad --remix-config schema-shendure.yml --output-filename subsampled_remixed.h5ad
cellxgene schema apply --source-h5ad Survey_of_human_embryonic_development-processed_w_disease.h5ad --remix-config schema-shendure.yml --output-filename Survey_of_human_embryonic_development-processed_remixed.h5ad
```

**3. Validate remixing**

Once the remixing is done the final step I did was to validate the schema was properly applied 

```bash
cellxgene schema validate subsampled_remixed.h5ad
cellxgene schema validate Survey_of_human_embryonic_development-processed_remixed.h5ad
```



### Curation notes and exceptions

- Samples were obtained by the University of Washington Birth Defects Research Laboratory (BDRL) under a protocol approved by the University of Washington Institutional Review Board
- Cell types were generally mapped to the most general ontology term when possible. For instance, a stromal cell from the lung was mapped to the ontology "stromal cell" and not to "stromal cell of lung"
- Cell types that didn't exist in the ontology look-up service were not mapped to an ontology term, and were assigned the same label as the original. For instance, "CSH1_CSH2 positive" cells.
organ
