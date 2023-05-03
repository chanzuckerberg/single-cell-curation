This file describes the expected format of *gene_info.yml*. It is used to track is used to track the current version of 
gene data to be download from various sources.

### FORMAT
Bellow is a description of the expected fields.

```
{id}:                 (object, required)                  used to identify the dataset
    file_description: (str, required)                     can only container characters valid in a file name.
    file_type:        (str, required)                     the file type to save the file as.
    version:          (str, optional, default is None)    use to track the current version of the data used
    url:              (str, required)                     the URL pointing at the data to download.
```

```
human:
  file_description: homo_sapiens
  file_type: gtf.gz
  version: 38
  url: http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz
mouse:
  file_description: mus_musculus
  file_type: gtf.gz
  version: M27
  url: http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.primary_assembly.annotation.gtf.gz
sars_cov_2:
  file_description: sars_cov_2
  file_type: gtf.gz
  url: ftp://ftp.ensemblgenomes.org/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf.gz
ercc:
  file_description: ercc
  file_type: txt
  url: https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095047.txt

```
