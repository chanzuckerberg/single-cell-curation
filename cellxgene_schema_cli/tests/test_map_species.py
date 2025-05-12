import logging
from tempfile import TemporaryDirectory
from unittest import mock

import numpy
import pandas as pd
import pytest
from cellxgene_schema.map_species import OntologyParser, map_species
from cellxgene_schema.validate import Validator
from fixtures.examples_validate import adata


@pytest.fixture
def zebrafish_adata():
    zebrafish_adata = adata.copy()
    obs = zebrafish_adata.obs
    obs.drop(columns=["cell_type_ontology_term_id", "tissue_ontology_term_id"], inplace=True)
    zebrafish_adata.uns["organism_ontology_term_id"] = "NCBITaxon:7955"
    obs["development_stage_ontology_term_id"] = "ZFS:0000016"
    obs["self_reported_ethnicity_ontology_term_id"] = "na"
    obs["tissue_type"] = "tissue"

    obs["development_stage_ontology_term_id"] = obs["development_stage_ontology_term_id"].astype("category")
    obs["self_reported_ethnicity_ontology_term_id"] = obs["self_reported_ethnicity_ontology_term_id"].astype("category")
    obs["tissue_type"] = obs["tissue_type"].astype("category")

    # single match
    obs["organism_cell_type_ontology_term_id"] = "ZFA:0005769"

    # single match
    obs["organism_tissue_ontology_term_id"] = "ZFA:0000047"

    zebrafish_adata.var = pd.DataFrame(
        [[False], [False], [False], [False], [False], [False], [False]],
        index=[
            "ENSDARG00000103202",
            "ENSDARG00000096156",
            "ENSDARG00000076160",
            "ENSDARG00000117163",
            "ENSDARG00000096187",
            "ENSDARG00000009657",
            "ENSDARG00000098510",
        ],
        columns=["feature_is_filtered"],
    )
    X = numpy.ones([zebrafish_adata.obs.shape[0], zebrafish_adata.var.shape[0]], dtype=numpy.float32)
    for i in range(zebrafish_adata.obs.shape[0]):
        for j in range(zebrafish_adata.var.shape[0]):
            X[i, j] = i + j
    zebrafish_adata.X = X
    zebrafish_adata.raw = zebrafish_adata.copy()
    zebrafish_adata.raw.var.drop("feature_is_filtered", axis=1, inplace=True)
    return zebrafish_adata


@pytest.fixture
def validator_with_adata(zebrafish_adata):
    validator = Validator()
    validator.adata = zebrafish_adata.copy()
    return validator


def test_map_species__valid_output(validator_with_adata):
    with TemporaryDirectory() as tmp:
        input_file = tmp + "input.h5ad"
        validator_with_adata.adata.copy().write_h5ad(input_file, compression="gzip")
        output_file = tmp + "output.h5ad"
        map_species(input_file, output_file)

        assert validator_with_adata.validate_adata(output_file) is True


def test_map_species_log_ouput(validator_with_adata, caplog):
    caplog.set_level(logging.INFO)
    obs = validator_with_adata.adata.obs
    # test single match, error match, multiple options
    obs.loc[obs.index[0], "organism_cell_type_ontology_term_id"] = "ZFA:0009384"
    obs.loc[obs.index[1], "organism_tissue_ontology_term_id"] = "ZFA:0"

    with TemporaryDirectory() as tmp:
        input_file = tmp + "input.h5ad"
        validator_with_adata.adata.copy().write_h5ad(input_file, compression="gzip")
        output_file = tmp + "output.h5ad"
        map_species(input_file, output_file)
        assert (
            "ZFA:0009384 has multiple closest matches for CL in its ancestry tree: ['CL:0000075', 'CL:0002077']"
            in caplog.text
        )
        assert "Setting ZFA:0005769 cell_type_ontology_term_id rows to single closest match: CL:0000018" in caplog.text
        assert (
            "ZFA:0 in organism_tissue_ontology_term_id is not a valid ontology term ID, cannot find UBERON map."
            in caplog.text
        )
        assert "Setting ZFA:0000047 tissue_ontology_term_id rows to single closest match: UBERON:0000004" in caplog.text

        # output should be invalid, per logs
        assert validator_with_adata.validate_adata(output_file) is False


def test_map_species_log_output_no_match(zebrafish_adata, caplog):
    with TemporaryDirectory() as tmp, mock.patch.object(
        OntologyParser, "get_closest_bridge_term_ids", mock.Mock(return_value=[])
    ):
        input_file = tmp + "input.h5ad"
        zebrafish_adata.copy().write_h5ad(input_file, compression="gzip")
        output_file = tmp + "output.h5ad"
        map_species(input_file, output_file)
        assert "ZFA:0005769 has no closest match for CL in its ancestry tree." in caplog.text
        assert "ZFA:0000047 has no closest match for UBERON in its ancestry tree." in caplog.text


def test_map_species__with_existing_cols(validator_with_adata):
    obs = validator_with_adata.adata.obs
    obs["cell_type_ontology_term_id"] = ""
    obs["tissue_ontology_term_id"] = ""
    # set-up organism_cell_type and organism_tissue term IDs with one of multiple matching CL/UBERON options
    # leave single match rows blank, to be set by map_species
    obs.loc[obs.index[0], "organism_cell_type_ontology_term_id"] = "ZFA:0009384"
    obs.loc[obs.index[0], "cell_type_ontology_term_id"] = "CL:0000075"
    obs.loc[obs.index[0], "organism_tissue_ontology_term_id"] = "ZFA:0005914"
    obs.loc[obs.index[0], "tissue_ontology_term_id"] = "UBERON:0001769"

    with TemporaryDirectory() as tmp:
        input_file = tmp + "input.h5ad"
        validator_with_adata.adata.copy().write_h5ad(input_file, compression="gzip")
        output_file = tmp + "output.h5ad"
        map_species(input_file, output_file)

        assert validator_with_adata.validate_adata(output_file) is True
