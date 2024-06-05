import pandas as pd
import pandera as pa
from cellxgene_ontology_guide.ontology_parser import OntologyParser
from pandera import Check, Column, DataFrameSchema, Index

ONTOLOGY_PARSER = OntologyParser(schema_version="v5.1.0")

uns_spatial_is_single = False

ASSAY_TO_SUSPENSION_TYPE = {
    "EFO:0030080": ["cell", "nucleus"],
    "EFO:0007045": ["nucleus"],
    "EFO:0700003": ["cell"],
    "EFO:0700004": ["cell"],
    "EFO:0010010": ["cell", "nucleus"],
    "EFO:0009294": ["cell"],
    "EFO:0008720": ["nucleus"],
    "EFO:0008722": ["cell", "nucleus"],
    "EFO:0700011": ["cell", "nucleus"],
    "EFO:0008780": ["cell", "nucleus"],
    "EFO:0008796": ["cell"],
    "EFO:0030002": ["cell"],
    "EFO:0008853": ["cell"],
    "EFO:0030026": ["nucleus"],
    "EFO:0010550": ["cell", "nucleus"],
    "EFO:0008919": ["cell"],
    "EFO:0010184": ["cell", "nucleus"],
    "EFO:0009918": ["na"],
    "EFO:0008939": ["nucleus"],
    "EFO:0030027": ["nucleus"],
    "EFO:0700000": ["na"],
    "EFO:0008994": ["na"],
    "EFO:0009919": ["cell", "nucleus"],
    "EFO:0008953": ["cell"],
    "EFO:0700010": ["cell", "nucleus"],
}


def check_assay_ontology_term_id(pandas_obj) -> bool:
    return (
        pandas_obj.isin(ONTOLOGY_PARSER.get_term_descendants("EFO:0002772")).all()
        or pandas_obj.isin(ONTOLOGY_PARSER.get_term_descendants("EFO:0010183")).all()
        or pandas_obj.isin(["EFO:0010961"]).all()
        or pandas_obj.isin(["EFO:0030062"]).all()
    )


def check_cell_type_ontology_term_id(df) -> bool:
    valid_cell_types = ONTOLOGY_PARSER.get_term_descendants("CL:0000000") + ["unknown"]

    if uns_spatial_is_single:
        condition1 = (
            (df["assay_ontology_term_id"] == "EFO:0010961")
            & (df["in_tissue"] == 0)
            & (df["cell_type_ontology_term_id"] != "unknown")
        )
    else:
        condition1 = pd.Series([False] * len(df))
    condition2 = ~df["cell_type_ontology_term_id"].isin(valid_cell_types)

    return not (condition1 | condition2).any()


FORBIDDEN_HANCESTRO_TERMS = [
    "HANCESTRO:0002",
    "HANCESTRO:0003",
    "HANCESTRO:0004",
    "HANCESTRO:0018",
    "HANCESTRO:0290",
    "HANCESTRO:0304",
    "HANCESTRO:0323",
    "HANCESTRO:0324",
    "HANCESTRO:0551",
    "HANCESTRO:0554",
    "HANCESTRO:0555",
    "HANCESTRO:0557",
    "HANCESTRO:0558",
    "HANCESTRO:0559",
    "HANCESTRO:0560",
    "HANCESTRO:0561",
    "HANCESTRO:0564",
    "HANCESTRO:0565",
    "HANCESTRO:0566",
    "HANCESTRO:0029",
    "HANCESTRO:0030",
    "HANCESTRO:0031",
    "HANCESTRO:0032",
    "HANCESTRO:0033",
    "HANCESTRO:0034",
]


def check_self_reported_ethnicity_ontology_term_id(row):
    organism = row["organism_ontology_term_id"]
    ethnicity = row["self_reported_ethnicity_ontology_term_id"]
    if organism == "NCBITaxon:9606":
        # For Homo sapiens, the value must be one or more comma-separated HANCESTRO terms in ascending lexical order
        # with no duplication of terms or "unknown" if unavailable.
        terms = ethnicity.split(",")
        return (
            terms == sorted(terms)  # sorted
            and len(terms) == len(set(terms))  # no duplicates
            and not any(term in FORBIDDEN_HANCESTRO_TERMS for term in terms)  # no forbidden terms
            and all(ONTOLOGY_PARSER.is_valid_term_id(term) for term in terms)
        )  # all valide ontology terms
    else:
        # For all other organisms the value must be "na".
        return ethnicity == "na"


def check_suspension_type(row):
    assay = row["assay_ontology_term_id"]
    suspension_type = row["suspension_type"]
    if assay in ASSAY_TO_SUSPENSION_TYPE:
        return suspension_type in ASSAY_TO_SUSPENSION_TYPE[assay]
    return True


def check_tissue_ontology_term_id(row):
    tissue_type = row["tissue_type"]
    tissue_ontology_term_id = row["tissue_ontology_term_id"]
    if tissue_type == "tissue":
        return ONTOLOGY_PARSER.is_valid_term_id(tissue_ontology_term_id, "UBERON")
    return True


# Define the schema for each column
index_schema = Index(pa.String, Check(lambda s: s.is_unique))
array_col_schema = Column(
    pa.Int, Check(lambda x: 0 <= x <= 127), required=False
)  # Must only be set if if assay_ontology_term_id is "EFO:0010961" for Visium Spatial Gene
# Expression and uns['spatial']['is_single'] is True; otherwise, this key MUST NOT be present.
array_row_schema = Column(
    pa.Int, Check(lambda x: 0 <= x <= 77), required=False
)  # Must only be set if if assay_ontology_term_id is "EFO:0010961" for Visium Spatial Gene
# Expression and uns['spatial']['is_single'] is True; otherwise, this key MUST NOT be present.
assay_ontology_term_id_schema = Column(pa.Category, Check(check_assay_ontology_term_id))
cell_type_ontology_term_id_schema = Column(
    pa.Category, checks=[Check.notin(["CL:0000255", "CL:0000257", "CL:0000548"])]
)
development_stage_ontology_term_id_schema = Column(pa.Category)
disease_ontology_term_id_schema = Column(
    pa.Category, Check.isin(["PATO:0000461", "MONDO:0000001", "MONDO:0021178"]), nullable=False
)
donor_id_schema = Column(pa.Category)
in_tissue_schema = Column(
    pa.Int, Check.isin([0, 1]), required=False
)  # Must only be set if if assay_ontology_term_id is "EFO:0010961" for Visium Spatial Gene
# Expression and uns['spatial']['is_single'] is True; otherwise, this key MUST NOT be present.
is_primary_data_schema = Column(
    pa.Bool
)  # This MUST be False if uns['spatial']['is_single'] is False. This MUST be True if this is the
# canonical instance of this cellular observation and False if not. This is commonly False for meta-analyses reusing
# data or for secondary views of data.
organism_ontology_term_id_schema = Column(
    pa.Category, Check.isin(ONTOLOGY_PARSER.get_term_descendants("NCBITaxon:33208"))
)
self_reported_ethnicity_ontology_term_id_schema = Column(pa.Category)
self_reported_ethnicity_ontology_term_id_schema = Column(pa.Category, Check.isin("na"))
sex_ontology_term_id_schema = Column(
    pa.Category, Check.isin(ONTOLOGY_PARSER.get_term_descendants("PATO:0001894") + ["unknown"])
)
suspension_type_schema = Column(pa.Category)
tissue_type_schema = Column(pa.Category, Check.isin(["tissue", "organoid", "cell culture"]))
tissue_ontology_term_id_schema = Column(pa.Category)

# Combine all the column schemas into a DataFrameSchema
obs_schema = DataFrameSchema(
    {
        "array_col": array_col_schema,
        "array_row": array_row_schema,
        "assay_ontology_term_id": assay_ontology_term_id_schema,
        "cell_type_ontology_term_id": cell_type_ontology_term_id_schema,
        "development_stage_ontology_term_id": development_stage_ontology_term_id_schema,
        "disease_ontology_term_id": disease_ontology_term_id_schema,
        "donor_id": donor_id_schema,
        "in_tissue": in_tissue_schema,
        "is_primary_data": is_primary_data_schema,
        "organism_ontology_term_id": organism_ontology_term_id_schema,
        "self_reported_ethnicity_ontology_term_id": self_reported_ethnicity_ontology_term_id_schema,
        "sex_ontology_term_id": sex_ontology_term_id_schema,
        "suspension_type": suspension_type_schema,
        "tissue_type": tissue_type_schema,
        "tissue_ontology_term_id": tissue_ontology_term_id_schema,
    },
    index=index_schema,
    checks=[
        Check(check_cell_type_ontology_term_id),
        Check(check_self_reported_ethnicity_ontology_term_id),
        Check(check_suspension_type),
        Check(check_tissue_ontology_term_id),
    ],
)

import anndata as ad
from tests.fixtures.examples_validate import h5ad_valid


def test_obs_schema():
    adata = ad.read_h5ad(h5ad_valid)
    obs_schema.validate(adata.obs, lazy=True)
