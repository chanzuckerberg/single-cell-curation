from unittest.mock import Mock, patch

import pytest

from scripts.compute_mappings.compute_tissue_and_cell_type_mappings import (
    AGraph,
    build_descendants_graph,
)


@pytest.fixture
def fake_thing():
    class FakeThingClass:
        def __init__(self, name):
            self.name = name

    return FakeThingClass


class TestGraphBuild:
    @patch("scripts.compute_mappings.compute_tissue_and_cell_type_mappings.list_direct_descendants")
    def test_build_descendants_graph(self, list_direct_descendants_mock: Mock, fake_thing):
        graph = AGraph()
        # return_values = [["1"], ["2", "3"]]
        return_values = [[fake_thing("1")], [fake_thing("2"), fake_thing("3")]]
        list_direct_descendants_mock.side_effect = lambda: return_values.pop(0)

        build_descendants_graph("0", graph)
