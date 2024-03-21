from cellxgene_schema import gencode

# For GeneChecker
invalid_species = ["Caenorhabditis elegans"]

valid_genes = {
    gencode.SupportedOrganisms.HOMO_SAPIENS: {"ENSG00000141510": ("TP53", 6836)},
    gencode.SupportedOrganisms.MUS_MUSCULUS: {"ENSMUSG00000059552": ("Trp53", 4045)},
}

valid_genes_same_name_diff_species = {
    gencode.SupportedOrganisms.HOMO_SAPIENS: {"ENSG00000166278": ("C2_ENSG00000166278", 7151)},
    gencode.SupportedOrganisms.MUS_MUSCULUS: {
        "ENSMUSG00000024371": ("C2_ENSMUSG00000024371", 3872),
    },
}

valid_genes_same_name_and_species = {
    gencode.SupportedOrganisms.MUS_MUSCULUS: {
        "ENSMUSG00000091071": ("1700030C10Rik_ENSMUSG00000091071", 126),
        "ENSMUSG00000099759": ("1700030C10Rik_ENSMUSG00000099759", 1735),
    },
}

invalid_genes = {
    gencode.SupportedOrganisms.HOMO_SAPIENS: ["ENSMUSG00000059552", ("GENE", 1000)],
    gencode.SupportedOrganisms.MUS_MUSCULUS: ["ENSG00000141510", ("GENE", 200)],
}

# For ontology checker
valid_ontologies = [
    "CL",
    "EFO",
    "HANCESTRO",
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
    "HANCESTRO": {
        "HANCESTRO:0005": "European",
        "HANCESTRO:0014": "Hispanic or Latin American",
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
