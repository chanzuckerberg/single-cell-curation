This file describes the expected format of *gene_info.yml*. It is used to track is used to track the current version of 
gene data to be download from various sources.

### FORMAT
Bellow is a description of the expected fields.

```
{id}:                 (object, required)                  used to identify the dataset
    description: (str, required)                     can only container characters valid in a file name.
    file_type:        (str, required)                     the file type to save the file as.
    version:          (str, optional, default is None)    use to track the current version of the data used
    url:              (str, required)                     the URL pointing at the data to download.
```
