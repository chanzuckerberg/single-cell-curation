from cellxgene_schema import gencode

# For GeneChecker
invalid_species = ["Caenorhabditis elegans"]

valid_genes = {
    gencode.SupportedOrganisms.HOMO_SAPIENS: {"ENSG00000141510": ("TP53", 2426)},
    gencode.SupportedOrganisms.MUS_MUSCULUS: {"ENSMUSG00000059552": ("Trp53", 1797)},
    gencode.SupportedOrganisms.DROSOPHILA_MELANOGASTER: {"RR45003_transposable_element": ("S{}RR4500", 1234)},
    gencode.SupportedOrganisms.DROSOPHILA_MELANOGASTER: {"FBgn0037293": ("RabGGTa", 1853)},
    gencode.SupportedOrganisms.CAENORHABDITIS_ELEGANS: {"WBGene00000003": ("aat-2", 1738)},
    gencode.SupportedOrganisms.CALLITHRIX_JACCHUS: {"ENSCJAG00000071296": ("U4", 141)},
    gencode.SupportedOrganisms.DANIO_RERIO: {"ENSDARG00000009657": ("FGFR1OP2", 1088)},
    gencode.SupportedOrganisms.GORILLA_GORILLA: {"ENSGGOG00000010861": ("CAMSAP2", 7438)},
    gencode.SupportedOrganisms.MACACA_FASCICULARIS: {"ENSMFAG00000001539": ("DFFB", 1174)},
    gencode.SupportedOrganisms.MACACA_MULATTA: {"ENSMMUG00000000634": ("ZNF692", 1944)},
    gencode.SupportedOrganisms.MICROCEBUS_MURINUS: {"ENSMICG00000026886": ("CIR1", 1807)},
    gencode.SupportedOrganisms.ORYCTOLAGUS_CUNICULUS: {"ENSOCUG00000025472": ("SNORD42", 67)},
    gencode.SupportedOrganisms.PAN_TROGLODYTES: {"ENSPTRG00000000799": ("HOOK1", 5839)},
    gencode.SupportedOrganisms.RATTUS_NORVEGICUS: {"ENSRNOG00000070901": ("Irgq", 6200)},
    gencode.SupportedOrganisms.SUS_SCROFA: {"ENSSSCG00000031382": ("C9orf40", 3815)},
}

valid_genes_same_name_diff_species = {
    gencode.SupportedOrganisms.HOMO_SAPIENS: {"ENSG00000166278": ("C2", 1876)},
    gencode.SupportedOrganisms.MUS_MUSCULUS: {
        "ENSMUSG00000024371": ("C2", 694),
    },
}

valid_genes_same_name_and_species = {
    gencode.SupportedOrganisms.MUS_MUSCULUS: {
        "ENSMUSG00000091071": ("1700030C10Rik", 126),
        "ENSMUSG00000099759": ("1700030C10Rik", 991),
    },
}

invalid_genes = {
    gencode.SupportedOrganisms.HOMO_SAPIENS: ["ENSMUSG00000059552", ("GENE", 1000)],
    gencode.SupportedOrganisms.MUS_MUSCULUS: ["ENSG00000141510", ("GENE", 200)],
    gencode.SupportedOrganisms.CAENORHABDITIS_ELEGANS: {"WBGene_00000003": ("aat-2", 1738)},
    gencode.SupportedOrganisms.CALLITHRIX_JACCHUS: {"ENSCJAG_00000071296": ("U4", 141)},
    gencode.SupportedOrganisms.DANIO_RERIO: {"ENSDARG_00000009657": ("fgfr1op2", 1088)},
    gencode.SupportedOrganisms.GORILLA_GORILLA: {"ENSGGOG_00000010861": ("CAMSAP2", 7438)},
    gencode.SupportedOrganisms.MACACA_FASCICULARIS: {"ENSMFAG_00000001539": ("DFFB", 1174)},
    gencode.SupportedOrganisms.MACACA_MULATTA: {"ENSMMUG_00000000634": ("ZNF692", 1944)},
    gencode.SupportedOrganisms.MICROCEBUS_MURINUS: {"ENSMICG_00000026886": ("CIR1", 1807)},
    gencode.SupportedOrganisms.ORYCTOLAGUS_CUNICULUS: {"ENSOCUG_00000025472": ("SNORD42", 67)},
    gencode.SupportedOrganisms.PAN_TROGLODYTES: {"ENSPTRG_00000000799": ("HOOK1", 5839)},
    gencode.SupportedOrganisms.RATTUS_NORVEGICUS: {"ENSRNOG_00000070901": ("Irgq1", 6116)},
    gencode.SupportedOrganisms.SUS_SCROFA: {"ENSSSCG_00000031382": ("C9orf40", 3815)},
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
        "HANCESTRO:0019": "Japanese",
        "HANCESTRO:0014": "Afghan",
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
