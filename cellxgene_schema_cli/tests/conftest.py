import anndata as ad
import pandas as pd
import pytest


@pytest.fixture
def atac_anndata():
    obs = pd.DataFrame(index=["A", "B", "C"])
    uns = {"organism_ontology_term_id": "NCBITaxon:9606"}
    obs["assay_ontology_term_id"] = ["EFO:0030059"] * 3
    var = pd.DataFrame()
    var["var_names"] = ["gene1", "gene2", "gene3"]
    X = pd.DataFrame(index=["A", "B", "C"], data=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    return ad.AnnData(obs=obs, var=var, X=X, uns=uns)
