from cellxgene_schema import atac_seq


def test_process_fragment():
    atac_seq.process_fragment(
        "fragment.tsv.gz",
        "anndata.h5ad",
        generate_index=True,
    )
