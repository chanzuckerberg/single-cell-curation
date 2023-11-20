from unittest.mock import Mock, patch

import pytest
from owlready2 import Ontology

from scripts.compute_mappings.compute_ancestor_mapping import RO__PART_OF, create_ancestors_mapping, get_ancestors


class FakeThingClass:
    def __init__(self, name, is_a=None, value=None, property=None):
        self.name = name
        self.is_a = is_a or []
        self.value = value
        self.property = property


class TestAncestorMapping:
    def test_get_ancestors(self):
        # Arrange
        A = FakeThingClass("obo.A_")
        B = FakeThingClass("obo.B_")
        C_A = FakeThingClass("C_to_A", value=A, property=FakeThingClass(RO__PART_OF))
        C_B = FakeThingClass("C_to_A", value=B, property=FakeThingClass(RO__PART_OF))
        # Should be skipped if missing 'value' attr
        C_NA = FakeThingClass("should be skipped", value=FakeThingClass("NA"), property=FakeThingClass(RO__PART_OF))
        delattr(C_NA, "value")  # Test obj does not have 'value' attr
        # Should be skipped if property.name is not RO__PART_OF
        C_NA_2 = FakeThingClass(
            "should be skipped",
            value=FakeThingClass("NA_2"),
            property=FakeThingClass("not_ro_part_of"),  # Test wrong type of relation
        )
        C = FakeThingClass("C_", is_a=[C_NA, C_NA_2, C_A, C_B])
        D = FakeThingClass("D_")
        E_C = FakeThingClass("C_to_A", value=C, property=FakeThingClass(RO__PART_OF))
        E_D = FakeThingClass("C_to_A", value=D, property=FakeThingClass(RO__PART_OF))
        E = FakeThingClass("E_", is_a=[E_C, E_D])
        return_values = [D, B, A, C, E]
        attrs = {"search_one.side_effect": lambda **x: return_values.pop()}
        ontology = Mock(spec=Ontology, **attrs)

        #
        # graph:
        #
        #   A   B
        #    \ /
        #     C   D
        #      \ /
        #       E
        #

        # Act
        ancestors = get_ancestors(ontology, "E:")

        # Assert
        assert list(ancestors) == ["E_", "C_", "obo.A_", "obo.B_", "D_"]  # Should not contain "NA" or "NA_2"

    @pytest.mark.parametrize(
        "classes, prefix", [([FakeThingClass("1_"), FakeThingClass("2_"), FakeThingClass("3_")], None)]
    )
    @patch("scripts.compute_mappings.compute_ancestor_mapping.get_ancestors")
    def test_create_ancestors_mapping(self, get_ancestors_mock: Mock, classes, prefix):
        # Arrange
        # print([c.name for c in classes])
        ontology = Mock(spec=Ontology)
        return_values = [["3_", "6_"], [], ["1_", "4_", "5:"]]
        get_ancestors_mock.side_effect = lambda *x: return_values.pop()
        expected_ancestors_dict = {"1:": ["1:", "4:", "5:"], "3:": ["3:", "6:"]}

        # Act
        ancestors_dict = create_ancestors_mapping(ontology, classes, prefix=prefix)  # type: ignore

        # Assert
        assert ancestors_dict == expected_ancestors_dict
