import tempfile

from cellxgene_schema.validate import validate

H5AD_FILE = "/data/dataset_34gb_unlabeled.h5ad"


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmp_file = "/".join([tmpdirname, "/data/labeled.h5ad"])
        validate(H5AD_FILE, tmp_file)
