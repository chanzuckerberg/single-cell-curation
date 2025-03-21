import tempfile

from cellxgene_schema.atac_seq import process_fragment

H5AD_FILE = "./data/revised.h5ad"
FRAGMENT_FILE = "./data/fragments.tsv.gz"

if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmp_file = "/".join([tmpdirname, "fragment.tsv.bgz"])
        process_fragment(FRAGMENT_FILE, H5AD_FILE, True, output_file=tmp_file)
