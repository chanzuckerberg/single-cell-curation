This file describes the expected format of *[gene_info.yml](./gene_info.yml)*. It is used to track the
current version of gene data to be download from various sources.

### FORMAT

Bellow is a description of the expected fields.

| Field       | type   | required | default | description                                                                                                                        |
|-------------|--------|----------|---------|------------------------------------------------------------------------------------------------------------------------------------|
| {id}        | object | True     |         | Used to identify the dataset.                                                                                                      |
| description | str    | True     |         | Can only contain characters valid in a file name.                                                                                  |
| version     | str    | False    | None    | Use to track the current version of the data used. A version must be present to support auto updating.                             |
| url         | str    | True     |         | The URL pointing at the data to download or a python formatting string that can be used to build the URL by adding in the version. |

### example

```
human:
  description: homo_sapiens
  version: '38'
  url: http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{version}/gencode.v{version}.primary_assembly.annotation.gtf.gz      
```

### Updating

There are two way to update *[gene_info.yml](./gene_info.yml)*.
The first approach is to run the recipe `make genes_update`. This is the automated approach and will pull the latest
version from Gencode and update `gene_info.yml` to reflect those changes. The second approach is to update the fields
`version` in *[gene_info.yml](./gene_info.yml)* to the desired version.
