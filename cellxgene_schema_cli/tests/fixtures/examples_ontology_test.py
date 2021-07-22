from cellxgene_schema import ontology

# For GeneChecker
invalid_species = ["Caenorhabditis elegans"]

valid_genes = {
    ontology.SupportedOrganisms.HOMO_SAPIENS: {"ENSG00000141510": "TP53"},
    ontology.SupportedOrganisms.MUS_MUSCULUS: {"ENSMUSG00000059552": "Trp53"},
}

invalid_genes = {
    ontology.SupportedOrganisms.HOMO_SAPIENS: ["ENSMUSG00000059552", "GENE"],
    ontology.SupportedOrganisms.MUS_MUSCULUS: ["ENSG00000141510", "GENE"],
}

# For ontology checker
valid_ontologies = [
    "CL",
    "EFO",
    "MONDO",
    "MmusDv",
    "NCBITaxon",
    "UBERON",
    "PATO",
]

invalid_ontologies = ["NOT_ONTOLOGY", "OOO"]

valid_terms = {
    "CL": {"CL:0000066": "epithelial cell", "CL:0000192": "smooth muscle cell"},
    "EFO": {
        "EFO:0009899": "10x 3' v2",
        "EFO:0009922": "10x 3' v3",
        "EFO:0011025": "10x 5' v1",
        "EFO:0008930": "Smart-seq",
        "EFO:0008931": "Smart-seq2",
    },
    "MONDO": {"MONDO:0100096": "COVID-19"},
    "MmusDv": {
        "MmusDv:0000062": "2 month-old stage",
        "MmusDv:0000003": "Theiler stage 01",
    },
    "NCBITaxon": {
        "NCBITaxon:9606": "Homo sapiens",
        "NCBITaxon:10090": "Mus musculus",
    },
    "UBERON": {"UBERON:0002048": "lung"},
    "PATO": {"PATO:0000461": "normal"},
}

invalid_terms = {
    "CL": ["EFO:0009899", "NO_TERM"],
    "EFO": ["UBERON:0002048", "NO_TERM"],
    "MONDO": ["EFO:0009899", "NO_TERM"],
    "MmusDv": ["EFO:0009899", "NO_TERM"],
    "NCBITaxon": ["EFO:0009899", "NO_TERM"],
    "UBERON": ["EFO:0009899", "NO_TERM"],
    "PATO": ["EFO:0009899", "NO_TERM"],
}
