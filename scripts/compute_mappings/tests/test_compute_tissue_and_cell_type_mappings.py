from typing import List
from unittest.mock import Mock, patch

from scripts.compute_mappings.compute_tissue_and_cell_type_mappings import (
    AGraph,
    build_descendants_and_parts_graph,
    build_descendants_graph,
)


class TestGraphBuild:
    class FakeThingClass:
        def __init__(self, name, axiom=False):
            self.name = name
            if axiom:
                self.__setattr__("Classes", None)

    def graph_structure_test(self, graph: AGraph, all_nodes: List[str]):
        assert all([graph.has_node(n) for n in all_nodes])
        assert graph.has_edge("0", "1")
        assert graph.has_edge("1", "2")
        assert graph.has_edge("1", "3")
        assert not graph.has_edge("0", "2")
        assert not graph.has_edge("0", "3")
        assert not graph.has_edge("2", "3")

    @patch("scripts.compute_mappings.compute_tissue_and_cell_type_mappings.list_direct_descendants")
    def test_build_descendants_graph(self, list_direct_descendants_mock: Mock):
        graph = AGraph()
        all_nodes = ["0", "1", "2", "3"]
        return_values = [["1"], ["2", "3"], [], []]  # Return two empty arrays to terminate recursion
        list_direct_descendants_mock.side_effect = lambda x: return_values.pop(0)

        build_descendants_graph("0", graph)

        #
        # graph:
        #
        #    0
        #    |
        #    1
        #   / \
        #  2   3
        #

        # Test graph formation
        assert list_direct_descendants_mock.call_count == 4
        assert [a.args[0] for a in list_direct_descendants_mock.call_args_list] == ["0", "1", "2", "3"]

        # Test graph structure
        self.graph_structure_test(graph, all_nodes)

    @patch("scripts.compute_mappings.compute_tissue_and_cell_type_mappings.list_direct_descendants_and_parts")
    def test_build_descendants_and_parts_graph(self, list_direct_descendants_and_parts_mock: Mock):
        graph = AGraph()
        all_nodes = ["0", "1", "2", "3"]
        return_values = [
            [[self.FakeThingClass("1")]],
            [
                [self.FakeThingClass("2")],
                [self.FakeThingClass("11", axiom=True)],
                [self.FakeThingClass("3")],
            ],  # Axiom should be ignored
            [[self.FakeThingClass("Nothing")]],  # Should be ignored
            [],  # This plus previous return value terminate recursion
        ]
        list_direct_descendants_and_parts_mock.side_effect = lambda x: return_values.pop(0)

        build_descendants_and_parts_graph("0", graph)

        assert list_direct_descendants_and_parts_mock.call_count == 4

        # Test graph structure
        self.graph_structure_test(graph, all_nodes)
