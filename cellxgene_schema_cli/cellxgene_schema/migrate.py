import anndata as ad

from . import utils

# fmt: off
# ONTOLOGY TERMS TO UPDATE ACROSS ALL DATASETS IN CORPUS
# Initialization is AUTOMATED for newly deprecated terms that have 'Replaced By' terms in their ontology files

# Curators should review the monthly 'Curator Report' and add deprecated term replacements to corresponding map if
# 'Replaced By' is not available for a deprecated term.

# If Curators have non-deprecated term changes to apply to all datasets in the corpus where applicable,
# add them here.
ONTOLOGY_TERM_MAPS = {
    "assay": {
    },
    "cell_type": {
        "CL:0000003": "unknown",
        "CL:0002371": "unknown",
    },
    "development_stage": {
    },
    "disease": {
    },
    "organism": {
    },
    "self_reported_ethnicity": {
    },
    "sex": {
    },
    "tissue": {
    },
}

DEPRECATED_FEATURE_IDS = [
    "ENSG00000269933",
    "ENSG00000261737",
    "ENSG00000259834",
    "ENSG00000256374",
    "ENSG00000263464",
    "ENSG00000203812",
    "ENSG00000272196",
    "ENSG00000272880",
    "ENSG00000284299",
    "ENSG00000270188",
    "ENSG00000287116",
    "ENSG00000237133",
    "ENSG00000224739",
    "ENSG00000227902",
    "ENSG00000239467",
    "ENSG00000272551",
    "ENSG00000280374",
    "ENSG00000284741",
    "ENSG00000236886",
    "ENSG00000229352",
    "ENSG00000286601",
    "ENSG00000227021",
    "ENSG00000259855",
    "ENSG00000228206",
    "ENSG00000273301",
    "ENSG00000271870",
    "ENSG00000237838",
    "ENSG00000286996",
    "ENSG00000269028",
    "ENSG00000286699",
    "ENSG00000273370",
    "ENSG00000261490",
    "ENSG00000272567",
    "ENSG00000270394",
    "ENSG00000280250",
    "ENSG00000272370",
    "ENSG00000272354",
    "ENSG00000251044",
    "ENSG00000272040",
    "ENSG00000271043",
    "ENSG00000288639",
    "ENSG00000182230",
    "ENSG00000285476",
    "ENSG00000204092",
    "ENSG00000261068",
    "ENSG00000236740",
    "ENSG00000236996",
    "ENSG00000255633",
    "ENSG00000232295",
    "ENSG00000269966",
    "ENSG00000271734",
    "ENSG00000236673",
    "ENSG00000280058",
    "ENSG00000227220",
    "ENSG00000236166",
    "ENSG00000112096",
    "ENSG00000285162",
    "ENSG00000286228",
    "ENSG00000205485",
    "ENSG00000237513",
    "ENSG00000285106",
    "ENSG00000226380",
    "ENSG00000270672",
    "ENSG00000225932",
    "ENSG00000244693",
    "ENSG00000268955",
    "ENSG00000272267",
    "ENSG00000253878",
    "ENSG00000259820",
    "ENSG00000226403",
    "ENSG00000288541",
    "ENSG00000233776",
    "ENSG00000269900",
    "ENSG00000283486",
    "ENSG00000272934",
    "ENSG00000272904",
    "ENSG00000212951",
    "ENSG00000261534",
    "ENSG00000237548",
    "ENSG00000282246",
    "ENSG00000239665",
    "ENSG00000256892",
    "ENSG00000249860",
    "ENSG00000271409",
    "ENSG00000224745",
    "ENSG00000261438",
    "ENSG00000231575",
    "ENSG00000260461",
    "ENSG00000255823",
    "ENSG00000256863",
    "ENSG00000254740",
    "ENSG00000254561",
    "ENSG00000282080",
    "ENSG00000256427",
    "ENSG00000287388",
    "ENSG00000276814",
    "ENSG00000277077",
    "ENSG00000280710",
    "ENSG00000215271",
    "ENSG00000258414",
    "ENSG00000258808",
    "ENSG00000277050",
    "ENSG00000273888",
    "ENSG00000258861",
    "ENSG00000259444",
    "ENSG00000244952",
    "ENSG00000137808",
    "ENSG00000279765",
    "ENSG00000273923",
    "ENSG00000262668",
    "ENSG00000232196",
    "ENSG00000288630",
    "ENSG00000261963",
    "ENSG00000256618",
    "ENSG00000221995",
    "ENSG00000226377",
    "ENSG00000286065",
    "ENSG00000273576",
    "ENSG00000267637",
    "ENSG00000282965",
    "ENSG00000225178",
    "ENSG00000279948",
    "ENSG00000273837",
    "ENSG00000286949",
    "ENSG00000256222",
    "ENSG00000276612",
    "ENSG00000279769",
    "ENSG00000280095",
    "ENSG00000280346",
    "ENSG00000279226",
    "ENSG00000278927",
    "ENSG00000278955",
    "ENSG00000273614",
    "ENSG00000277352",
    "ENSG00000161149",
    "ENSG00000285762",
    "ENSG00000239446",
    "ENSG00000288546",
    "ENSG00000256045",
    "ENSG00000228906",
    "ENSG00000228139",
    "ENSG00000261773",
    "ENSG00000278198",
    "ENSG00000273496",
    "ENSG00000277666",
    "ENSG00000278782",
    "ENSG00000277761",
    "ENSMUSG00000022591",
    "ENSMUSG00000094127",
    "ENSMUSG00000066936",
    "ENSMUSG00000116275",
    "ENSMUSG00000091312",
    "ENSMUSG00000098794",
    "ENSMUSG00000079353",
    "ENSMUSG00000096240",
    "ENSMUSG00000079286",
    "ENSMUSG00000085431",
    "ENSMUSG00000075015",
    "ENSMUSG00000075014",
    "ENSMUSG00000078091",
    "ENSMUSG00000075006",
    "ENSMUSG00000079175",
    "ENSMUSG00000079171",
    "ENSMUSG00000079170",
    "ENSMUSG00000079169",
    "ENSMUSG00000090353",
    "ENSMUSG00000100963",
    "ENSMUSG00000079039",
    "ENSMUSG00000078197",
    "ENSMUSG00000074735",
    "ENSMUSG00000092123",
    "ENSMUSG00000078912",
    "ENSMUSG00000090625",
    "ENSMUSG00000094487",
    "ENSMUSG00000094962",
    "ENSMUSG00000053706",
    "ENSMUSG00000090408",
    "ENSMUSG00000105204",
    "ENSMUSG00000090441",
    "ENSMUSG00000078620",
    "ENSMUSG00000078590",
    "ENSMUSG00000094856",
    "ENSMUSG00000096930",
    "ENSMUSG00000092345",
    "ENSMUSG00000097078",
    "ENSMUSG00000085147",
    "ENSMUSG00000095386",
    "ENSMUSG00000094958",
    "ENSMUSG00000073682",
    "ENSMUSG00000093574",
    "ENSMUSG00000095346",
    "ENSMUSG00000079439",
    "ENSMUSG00000079438",
    "ENSMUSG00000079416",
    "ENSMUSG00000091096",
    "ENSMUSG00000072693",
    "ENSMUSG00000105875",
    "ENSMUSG00000106198",
    "ENSMUSG00000079511",
    "ENSMUSG00000068181",
    "ENSMUSG00000092004",
    "ENSMUSG00000092047",
    "ENSMUSG00000079264",
    "ENSMUSG00000067292",
    "ENSMUSG00000116184",
    "ENSMUSG00000074210",
    "ENSMUSG00000094462",
    "ENSMUSG00000066378",
    "ENSMUSG00000091441",
    "ENSMUSG00000067627",
    "ENSMUSG00000074473",
    "ENSMUSG00000074302",
    "ENSMUSG00000079016",
    "ENSMUSG00000091041",
    "ENSMUSG00000092329",
    "ENSMUSG00000102141",
    "ENSMUSG00000078840",
    "ENSMUSG00000091028",
    "ENSMUSG00000096385",
    "ENSMUSG00000096519",
    "ENSMUSG00000074564",
    "ENSMUSG00000095547",
    "ENSMUSG00000095186",
    "ENSMUSG00000095891",
    "ENSMUSG00000096736",
    "ENSMUSG00000096201",
    "ENSMUSG00000091159",
    "ENSMUSG00000066810",
    "ENSMUSG00000094204",
    "ENSMUSG00000079433",
    "ENSMUSG00000096734",
    "ENSMUSG00000078488",
    "ENSMUSG00000078481",
    "ENSMUSG00000097704",
    "ENSMUSG00000091731",
    "ENSMUSG00000097854",
    "ENSMUSG00000069518",
    "ENSMUSG00000079010",
    "ENSMUSG00000096923",
    "ENSMUSG00000091195",
    "ENSMUSG00000098292",
    "ENSMUSG00000096083",
    "ENSMUSG00000092157",
    "ENSMUSG00000114046",
    "ENSMUSG00000079061",
    "ENSMUSG00000113103",
    "ENSMUSG00000078984",
    "ENSMUSG00000090863",
    "ENSMUSG00000090853",
    "ENSMUSG00000091768",
    "ENSMUSG00000114923",
    "ENSMUSG00000091109",
    "ENSMUSG00000067122",
    "ENSMUSG00000115235",
    "ENSMUSG00000091306",
    "ENSMUSG00000090889",
    "ENSMUSG00000096463",
    "ENSMUSG00000094296",
    "ENSMUSG00000079024",
    "ENSMUSG00000091604",
    "ENSMUSG00000079740",
    "ENSMUSG00000094472",
    "ENSMUSG00000095464",
    "ENSMUSG00000079546",
    "ENSMUSG00000094030",
    "ENSMUSG00000067929",
    "ENSMUSG00000095330",
    "ENSMUSG00000045506",
    "ENSMUSG00000079333",
    "ENSMUSG00000053861",
    "ENSMUSG00000117732",
    "ENSMUSG00000096597",
    "ENSMUSG00000116166",
    "ENSMUSG00000067085",
    "ENSMUSG00000073291",
    "ENSMUSG00000073290",
    "ENSMUSG00000095316",
    "ENSMUSG00000092463",
    "ENSMUSG00000079600",
    "ENSMUSG00000095693",
    "ENSMUSG00000101725",
    "ENSMUSG00000096850",
]

# Manually curated migration of v38 terms to v44
# https://github.com/chanzuckerberg/single-cell-curation/issues/572
GENCODE_MAPPER = {
    'ENSG00000112096': 'ENSG00000291237',
    'ENSG00000277077': 'ENSG00000289133',
    'ENSG00000215271': 'ENSG00000290292',
    'ENSG00000225932': 'ENSG00000288784',
    'ENSG00000244693': 'ENSG00000289604',
    'ENSG00000203812': 'ENSG00000288825',
    'ENSG00000272196': 'ENSG00000288859',
    'ENSG00000263464': 'ENSG00000288867',
    'ENSG00000256374': 'ENSG00000289549',
    'ENSG00000288546': 'ENSG00000289084',
    'ENSG00000205485': 'ENSG00000290873'
}
# fmt: on


def migrate(input_file, output_file, collection_id, dataset_id):
    print(f"Converting {input_file} into {output_file}")

    dataset = ad.read_h5ad(input_file, backed="r")
    if dataset.raw is not None and DEPRECATED_FEATURE_IDS:
        dataset = dataset.to_memory()

    # AUTOMATED, DO NOT CHANGE
    for ontology_name, deprecated_term_map in ONTOLOGY_TERM_MAPS.items():
        utils.replace_ontology_term(dataset.obs, ontology_name, deprecated_term_map)

    # CURATOR-DEFINED, DATASET-SPECIFIC UPDATES
    # Use the template below to define dataset and collection specific ontology changes. Will only apply to dataset
    # if it matches a condition.
    # If no such changes are needed, leave blank
    # Examples:
    # if dataset_id == "<dataset_1_id>":
    #   <no further logic necessary>
    #   utils.replace_ontology_term(df, <ontology_name>, {"term_to_replace": "replacement_term", ...})
    # elif dataset_id == "<dataset_2_id>":
    #   <custom transformation logic beyond scope of util functions>
    # elif collection_id == "<collection_1_id>":
    #   <no further logic necessary>
    #   utils.replace_ontology_term(df, <ontology_name>, {"term_to_replace": "replacement_term", ...})
    # elif collection_id == "<collection_2_id>":
    #   <custom transformation logic beyond scope of replace_ontology_term>
    # ...

    # Delete any uns keys with an empty value, logic taken from:
    # https://github.com/chanzuckerberg/single-cell-curation/blob/43f891005fb9439dbbb747fa0df8f0435ebf3f7c/cellxgene_schema_cli/cellxgene_schema/validate.py#L761-L762
    for key, value in list(dataset.uns.items()):
        if value is not None and type(value) is not bool and len(value) == 0:
            del dataset.uns[key]

    if "X_.umap_MinDist_0.2_N_Neighbors_15" in dataset.obsm:
        # Applies to dataset ids 63bb6359-3945-4658-92eb-3072419953e4 and 9b188f26-c8e1-4a78-af15-622a35a371fc
        dataset.obsm["X_umap_MinDist_0.2_N_Neighbors_15"] = dataset.obsm["X_.umap_MinDist_0.2_N_Neighbors_15"]
        del dataset.obsm["X_.umap_MinDist_0.2_N_Neighbors_15"]

    if collection_id == "91c8e321-566f-4f9d-b89e-3a164be654d5":
        utils.map_ontology_term(
            dataset.obs,
            "cell_type",
            "author_cell_type",
            {
                "L4-5IT_RORB_TSHZ2": "CL:4030062",
                "L2-4IT_CUX2": "CL:4030059",
                "L4-5IT_RORB_LRRK1": "CL:4030062",
                "L4-5IT_RORB_ARHGAP15": "CL:4030062",
                "L6IT_THEMIS_LINC00343": "CL:4030065",
                "L6IT_THEMIS_CUX1": "CL:4030065",
                "L3-5IT_RORB_PLCH1": "CL:4030061",
            },
        )

    if collection_id == "e1fa9900-3fc9-4b57-9dce-c95724c88716":
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:4023040": "CL:4030059"})

    if collection_id == "4a9fd4d7-d870-4265-89a5-ad51ab811d89":
        utils.replace_ontology_term(dataset.obs, "assay", {"EFO:0008913": "EFO:0022490"})

    if collection_id == "03608e22-227a-4492-910b-3cb3f16f952e":
        utils.replace_ontology_term(dataset.obs, "disease", {"MONDO:1011336": "MONDO:1010239"})
        utils.replace_ontology_term(dataset.obs, "assay", {"EFO:0008930": "EFO:0022488"})

    if collection_id == "59c9ecfe-c47d-4a6a-bab0-895cc0c1942b":
        utils.map_ontology_term(
            dataset.obs,
            "cell_type",
            "cellID2",
            {
                "L5/6 excitatory neuron": "CL:4030067",
                "L5/6 CC excitatory neuron": "CL:4030067",
                "L2/3 excitatory neuron": "CL:4030059",
                "L4 excitatory neuron": "CL:4030063",
            },
        )

    if dataset.uns["title"] == "Major cell cluster: CNS neurons":
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0003001": "CL:4033053"})
    if dataset.uns["title"] == "Major cell cluster: Epithelial cells":
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0000068": "CL:4030066"})
    if dataset.uns["title"] == "Whole dataset: Raw counts only":
        utils.replace_ontology_term(dataset.obs, "cell_type", {"CL:0003001": "CL:4033053", "CL:0000068": "CL:4030066"})

    # Manually remap v38 terms to v44
    dataset = utils.remap_deprecated_features(adata=dataset, remapped_features=GENCODE_MAPPER)

    # AUTOMATED, DO NOT CHANGE -- IF GENCODE UPDATED, DEPRECATED FEATURE FILTERING ALGORITHM WILL GO HERE.
    if DEPRECATED_FEATURE_IDS:
        dataset = utils.remove_deprecated_features(adata=dataset, deprecated=DEPRECATED_FEATURE_IDS)

    dataset.write(output_file, compression="gzip")
