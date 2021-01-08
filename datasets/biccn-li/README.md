# Li data curation (BICCN)

## Project details

Last updated: Mon Jan 8th 09:54:32 2021

Author: Maximilian Lombardo

Project title: An Atlas of Gene Regulatory Elements in Adult Mouse Cerebrum

DOI: [https://doi.org/10.1101/2020.05.10.087585](https://doi.org/10.1101/2020.05.10.087585)

Reqs:

- scanpy
- pandas
- numpy
- cellxgene with schema functions

Note:
- snucATACseq profiles were generated and split into seperate data objects representing Gabaergic Neurons (GabaN), Glutmatergic Neurons (GlutaN), and Non-Neuronal Cells (NonN)

## Curation steps

Clone this repo and download data from http://catlas.org/mousebrain/#!/

- Data can be manually downloaded by navigating to the CellBrowser Tab and selecting the dataset of choice (i.e. `NonN` - there are two entries for each dataset,
but the relevant selection is the first of the two). You can then navigate to the datadownload tab and download the relevant files (exprMatrix.tsv.gz, meta.tsv, and umap.coords.tsv.gz).
- Alternatively, one can download these data into seperate directories (pertaining to each major cell type) by running the following commands

### NonN

From your working directory (i.e. `./LiCuration2020/Data/`), download files pertaining to normalized "gene activity scores", nuclei metadata, and UMAP embeddings

```bash
mkdir NonN
cd NonN
wget http://catlas.org/mousebrain/cellbrowser/NonN_all/exprMatrix.tsv.gz
wget http://catlas.org/mousebrain/cellbrowser/NonN_all/meta.tsv
wget http://catlas.org/mousebrain/cellbrowser/NonN_all/umap.coords.tsv.gz
```

### GabaN
Do the same as above for the Gabaergic Neurons

```bash
mkdir GabaN
cd GabaN
wget http://catlas.org/mousebrain/cellbrowser/GABA_all/exprMatrix.tsv.gz
wget http://catlas.org/mousebrain/cellbrowser/GABA_all/meta.tsv
wget http://catlas.org/mousebrain/cellbrowser/GABA_all/umap.coords.tsv.gz
```

### GlutaN

And Finally for the Glutamatergic Neurons

```bash
mkdir GlutaN
cd GlutaN
wget http://catlas.org/mousebrain/cellbrowser/Glutamate_all/exprMatrix.tsv.gz
wget http://catlas.org/mousebrain/cellbrowser/Glutamate_all/meta.tsv
wget http://catlas.org/mousebrain/cellbrowser/Glutamate_all/umap.coords.tsv.gz
```


## Object Creation and annotation

Objects can then be created using the `LiCuration.py` script included in this repo (paths are not relative and need to be changed to the appropriate locations on your machine).
Codes pertaining to cell type and brain region can be translated using the Supplementary tables provided in this repo (already done in the provided script). These Supplementary Tables are directly provided by the authors.
A third conversion table (`NonNCellTypes_CLOntologies_byMax.tsv`) was manually created to map CL ontologies to the Non Neuronal Cell types and to store that information as cell/nuclei metadata.

## Apply cellxgene schema to newly created objects

Cellxgene schema can be applied by using the provided config files (for each object, found in `.datasets/biccn-li/objects/`) and by running the following commands in their respective directories

### NonN
```bash
cellxgene schema apply --source-h5ad non_neuronal_annotated.h5ad --remix-config config.yaml --output-filename non_neuronal_remixed.h5ad
cellxgene schema validate non_neuronal_remixed.h5ad
```

### GabaN
```bash
cellxgene schema apply --source-h5ad gaba_neurons_annotated.h5ad --remix-config config.yaml --output-filename gaba_neurons_remixed.h5ad
cellxgene schema validate gaba_neurons_remixed.h5ad
```

### GlutaN
```bash
cellxgene schema apply --source-h5ad gluta_neurons_annotated.h5ad --remix-config config.yaml --output-filename gluta_neurons_remixed.h5ad
cellxgene schema validate gluta_neurons_remixed.h5ad
```
