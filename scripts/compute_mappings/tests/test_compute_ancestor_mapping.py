from unittest.mock import Mock

from owlready2 import Ontology

from scripts.compute_mappings.compute_ancestor_mapping import RO__PART_OF, get_ancestors


class TestAncestorMapping:
    class FakeThingClass:
        def __init__(self, name, is_a=None, value=None, property=None):
            self.name = name
            self.is_a = is_a or []
            self.value = value
            self.property = property

    #
    # graph:
    #
    #   A   B
    #    \ /
    #     C   D
    #      \ /
    #       E
    #

    # @patch("scripts.compute_mappings.compute_ancestor_mapping.Ontology.search_one")
    def test_get_ancestors(self):
        # Arrange
        A = self.FakeThingClass("obo.A_")
        B = self.FakeThingClass("obo.B_")
        C_A = self.FakeThingClass("C_to_A", value=A, property=self.FakeThingClass(RO__PART_OF))
        C_B = self.FakeThingClass("C_to_A", value=B, property=self.FakeThingClass(RO__PART_OF))

        # Should be skipped if missing 'value' attr
        C_NA = self.FakeThingClass(
            "should be skipped", value=self.FakeThingClass("NA"), property=self.FakeThingClass(RO__PART_OF)
        )
        delattr(C_NA, "value")

        # Should be skipped if property.name is not RO__PART_OF
        C_NA_2 = self.FakeThingClass(
            "should be skipped", value=self.FakeThingClass("NA_2"), property=self.FakeThingClass("not_ro_part_of")
        )

        C = self.FakeThingClass("C_", is_a=[C_NA, C_NA_2, C_A, C_B])
        D = self.FakeThingClass("D_")
        E_C = self.FakeThingClass("C_to_A", value=C, property=self.FakeThingClass(RO__PART_OF))
        E_D = self.FakeThingClass("C_to_A", value=D, property=self.FakeThingClass(RO__PART_OF))
        E = self.FakeThingClass("E_", is_a=[E_C, E_D])

        return_values = [D, B, A, C, E]
        attrs = {"search_one.side_effect": lambda **x: return_values.pop()}
        ontology = Mock(spec=Ontology, **attrs)

        # Act
        ancestors = get_ancestors(ontology, "E:")

        # Assert
        assert list(ancestors) == ["E_", "C_", "obo.A_", "obo.B_", "D_"]
        print(ontology.call_args_list)

    # def test_create_ancestors_mapping(self):
    #     create_ancestors_mapping()
