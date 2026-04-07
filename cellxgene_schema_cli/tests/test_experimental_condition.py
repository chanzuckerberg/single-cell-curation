"""
Tests for schema 7.1.0 experimental_condition fields:
  - experimental_condition_ontology_term_id (curator-annotated, optional)
  - experimental_condition (Discover-annotated label, reserved)
  - perturbation_types (Discover-annotated, cross-column derived, reserved)
"""

import tempfile

from cellxgene_schema.validate import validate
from cellxgene_schema.write_labels import AnnDataLabelAppender
from fixtures.examples_validate import (
    adata_ec_absent,
    adata_ec_invalid_all_na,
    adata_ec_invalid_efo_not_allowed,
    adata_ec_invalid_forbidden_chebi,
    adata_ec_invalid_multi_unsorted,
    adata_ec_invalid_reserved_label_column,
    adata_ec_invalid_reserved_perturbation_types,
    adata_ec_valid_anti_uniprot,
    adata_ec_valid_chebi,
    adata_ec_valid_diet_descendant,
    adata_ec_valid_diet_root,
    adata_ec_valid_multi,
    adata_ec_valid_temperature,
    adata_ec_valid_uniprot,
    adata_ec_valid_with_gp,
    adata_gene_perturbations,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _validate(adata):
    with tempfile.TemporaryDirectory() as tmp:
        p = tmp + "/temp.h5ad"
        adata.copy().write_h5ad(p, compression="gzip")
        success, errors, _ = validate(p)
        return success, errors


def _write_labels(adata):
    import anndata as ad

    with tempfile.TemporaryDirectory() as tmp:
        out = tmp + "/labeled.h5ad"
        appender = AnnDataLabelAppender(adata.copy())
        appender.write_labels(out)
        return ad.read_h5ad(out)


# ---------------------------------------------------------------------------
# Validation — valid cases
# ---------------------------------------------------------------------------


class TestExperimentalConditionValid:
    def test_absent_column_is_valid(self):
        """Field is optional; omitting it entirely must pass validation."""
        success, errors = _validate(adata_ec_absent)
        assert success, errors

    def test_chebi_term(self):
        success, errors = _validate(adata_ec_valid_chebi)
        assert success, errors

    def test_anti_uniprot_term(self):
        """anti-uniprot: prefix should be stripped before ontology validation."""
        success, errors = _validate(adata_ec_valid_anti_uniprot)
        assert success, errors

    def test_uniprot_term(self):
        success, errors = _validate(adata_ec_valid_uniprot)
        assert success, errors

    def test_temperature_term(self):
        success, errors = _validate(adata_ec_valid_temperature)
        assert success, errors

    def test_diet_root_term(self):
        success, errors = _validate(adata_ec_valid_diet_root)
        assert success, errors

    def test_diet_descendant_term(self):
        success, errors = _validate(adata_ec_valid_diet_descendant)
        assert success, errors

    def test_multi_term_sorted(self):
        success, errors = _validate(adata_ec_valid_multi)
        assert success, errors

    def test_ec_with_genetic_perturbation(self):
        """ec + gp together must be valid (perturbation_types will be derived)."""
        success, errors = _validate(adata_ec_valid_with_gp)
        assert success, errors


# ---------------------------------------------------------------------------
# Validation — invalid cases
# ---------------------------------------------------------------------------


class TestExperimentalConditionInvalid:
    def test_all_na_forbidden(self):
        """Column present with every row == 'na' must fail (forbidden_when_all_na)."""
        success, errors = _validate(adata_ec_invalid_all_na)
        assert not success
        assert any("experimental_condition_ontology_term_id" in e for e in errors)

    def test_forbidden_chebi_term(self):
        """CHEBI:23367 (molecular entity) is explicitly forbidden."""
        success, errors = _validate(adata_ec_invalid_forbidden_chebi)
        assert not success
        assert any("CHEBI:23367" in e for e in errors)

    def test_efo_term_not_in_allowed_list(self):
        """An EFO term that is neither EFO:0001702 nor a diet descendant is not allowed."""
        success, errors = _validate(adata_ec_invalid_efo_not_allowed)
        assert not success
        assert any("EFO:0009899" in e for e in errors)

    def test_multi_term_unsorted(self):
        """Multi-term values must be in sorted order."""
        success, errors = _validate(adata_ec_invalid_multi_unsorted)
        assert not success

    def test_reserved_experimental_condition_column(self):
        """Curator must not manually set obs['experimental_condition'] (it is reserved)."""
        success, errors = _validate(adata_ec_invalid_reserved_label_column)
        assert not success
        assert any("experimental_condition" in e for e in errors)

    def test_reserved_perturbation_types_column(self):
        """Curator must not manually set obs['perturbation_types'] (it is reserved)."""
        success, errors = _validate(adata_ec_invalid_reserved_perturbation_types)
        assert not success
        assert any("perturbation_types" in e for e in errors)


# ---------------------------------------------------------------------------
# Label writing — experimental_condition label column
# ---------------------------------------------------------------------------


class TestExperimentalConditionLabels:
    def test_chebi_label_added(self):
        """write_labels should add obs['experimental_condition'] with human-readable label."""
        labeled = _write_labels(adata_ec_valid_chebi)
        assert "experimental_condition" in labeled.obs.columns
        assert labeled.obs["experimental_condition"].iloc[0] == "acetylsalicylic acid"

    def test_anti_uniprot_label_prefixed(self):
        """anti- prefix must be retained in the label (anti-CD180_HUMAN)."""
        labeled = _write_labels(adata_ec_valid_anti_uniprot)
        assert "experimental_condition" in labeled.obs.columns
        assert labeled.obs["experimental_condition"].iloc[0] == "anti-CD180_HUMAN"

    def test_uniprot_label_added(self):
        labeled = _write_labels(adata_ec_valid_uniprot)
        assert "experimental_condition" in labeled.obs.columns
        assert labeled.obs["experimental_condition"].iloc[0] == "CD180_HUMAN"

    def test_absent_column_no_label_added(self):
        """If experimental_condition_ontology_term_id is absent, no label column should be added."""
        labeled = _write_labels(adata_ec_absent)
        assert "experimental_condition" not in labeled.obs.columns

    def test_multi_term_label(self):
        """Multi-term values should produce ' || '-joined labels."""
        labeled = _write_labels(adata_ec_valid_multi)
        label = labeled.obs["experimental_condition"].iloc[0]
        assert " || " in label
        assert "acetylsalicylic acid" in label
        assert "diet" in label


# ---------------------------------------------------------------------------
# Label writing — perturbation_types derivation
# ---------------------------------------------------------------------------


class TestPerturbationTypesLabels:
    def test_no_perturbation_columns_no_perturbation_types(self):
        """Without ec or gp columns, perturbation_types must not be added."""
        labeled = _write_labels(adata_ec_absent)
        assert "perturbation_types" not in labeled.obs.columns

    def test_chebi_ec_gives_chemical_type(self):
        labeled = _write_labels(adata_ec_valid_chebi)
        assert "perturbation_types" in labeled.obs.columns
        assert labeled.obs["perturbation_types"].iloc[0] == "chemical"

    def test_anti_uniprot_ec_gives_protein_type(self):
        labeled = _write_labels(adata_ec_valid_anti_uniprot)
        assert "perturbation_types" in labeled.obs.columns
        assert labeled.obs["perturbation_types"].iloc[0] == "protein"

    def test_uniprot_ec_gives_protein_type(self):
        labeled = _write_labels(adata_ec_valid_uniprot)
        assert "perturbation_types" in labeled.obs.columns
        assert labeled.obs["perturbation_types"].iloc[0] == "protein"

    def test_temperature_ec_gives_temperature_type(self):
        labeled = _write_labels(adata_ec_valid_temperature)
        assert "perturbation_types" in labeled.obs.columns
        assert labeled.obs["perturbation_types"].iloc[0] == "temperature"

    def test_diet_ec_gives_diet_type(self):
        labeled = _write_labels(adata_ec_valid_diet_root)
        assert "perturbation_types" in labeled.obs.columns
        assert labeled.obs["perturbation_types"].iloc[0] == "diet"

    def test_diet_descendant_gives_diet_type(self):
        labeled = _write_labels(adata_ec_valid_diet_descendant)
        assert "perturbation_types" in labeled.obs.columns
        assert labeled.obs["perturbation_types"].iloc[0] == "diet"

    def test_gp_only_gives_genetic_type(self):
        labeled = _write_labels(adata_gene_perturbations)
        assert "perturbation_types" in labeled.obs.columns
        assert all(labeled.obs["perturbation_types"] == "genetic")

    def test_ec_and_gp_gives_combined_type(self):
        """CHEBI ec + genetic gp should produce 'chemical || genetic'."""
        labeled = _write_labels(adata_ec_valid_with_gp)
        assert "perturbation_types" in labeled.obs.columns
        pt = labeled.obs["perturbation_types"].iloc[0]
        assert pt == "chemical || genetic"

    def test_multi_ec_gives_multiple_types(self):
        """CHEBI + diet ec should give 'chemical || diet'."""
        labeled = _write_labels(adata_ec_valid_multi)
        assert "perturbation_types" in labeled.obs.columns
        pt = labeled.obs["perturbation_types"].iloc[0]
        assert pt == "chemical || diet"
