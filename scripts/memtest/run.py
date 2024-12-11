import tempfile

from cellxgene_schema.validate import validate

H5AD_FILE = "/memtest/data/slow_cxg.h5ad"


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmp_file = "/".join([tmpdirname, "/data/labeled.h5ad"])
        validate(H5AD_FILE, tmp_file)
