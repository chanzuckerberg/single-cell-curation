from typing import List
from unittest.mock import Mock, patch

import pytest

from scripts.compute_mappings.compute_tissue_and_cell_type_mappings import (
    AGraph,
    build_ancestors_set,
    build_descendants_and_parts_graph,
    build_descendants_graph,
    build_descendants_set,
    write_ancestors_by_entity,
    write_descendants_by_entity,
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

        # Test graph formation
        assert list_direct_descendants_and_parts_mock.call_count == 4
        assert [a.args[0] for a in list_direct_descendants_and_parts_mock.call_args_list] == ["0", "1", "2", "3"]

        # Test graph structure
        self.graph_structure_test(graph, all_nodes)

    def get_graph(self) -> AGraph:
        graph = AGraph()
        edges = [
            ("0_", "3_"),
            ("0_", "00_ (organoid)"),
            ("1_", "4_"),
            ("2_", "4_"),
            ("2_", "5_"),
            ("2_", "6_"),
            ("3_", "7_"),
            ("4_", "8_"),
            ("6_", "9_"),
        ]
        [graph.add_edge(x, y) for x, y in edges]

        #
        # graph:
        #
        #               0    1   2
        #              / \    \ /|\
        #     (organoid)  3    4 5 6
        #                  \    \   \
        #                   7    8   9
        #

        return graph

    @pytest.mark.parametrize(
        "entity, expected_ancestors_set",
        [
            ("0_", {"0_"}),
            ("3_", {"3_", "0_"}),
            ("7_", {"7_", "3_", "0_"}),
            ("8_", {"8_", "4_", "1_", "2_"}),
            ("9_", {"9_", "6_", "2_"}),
            ("5_", {"5_", "2_"}),
        ],
    )
    def test_build_ancestors_set(self, entity, expected_ancestors_set):
        ancestors_set = set()
        build_ancestors_set(str(entity), self.get_graph(), ancestors_set)
        assert ancestors_set == expected_ancestors_set

    @pytest.mark.parametrize(
        "entity, expected_descendants_set",
        [
            ("0_", {"3_", "7_", "00_ (organoid)"}),
            ("3_", {"7_"}),
            ("7_", set()),
            ("1_", {"8_", "4_"}),
            ("2_", {"9_", "6_", "5_", "4_", "8_"}),
            ("4_", {"8_"}),
        ],
    )
    def test_build_descendants_set(self, entity, expected_descendants_set):
        descendants_set = set()
        graph = self.get_graph()
        build_descendants_set(str(entity), graph, descendants_set)
        assert descendants_set == expected_descendants_set

    @pytest.mark.parametrize(
        "entities, expected_ancestors_dict",
        [
            (["0_"], {"0:": ["0:"]}),
            (["3_", "0_"], {"3:": ["3:", "0:"], "0:": ["0:"]}),
            (["7_", "6_", "2_"], {"7:": ["7:", "3:", "0:"], "6:": ["6:", "2:"], "2:": ["2:"]}),
            (["8_", "5_"], {"8:": ["8:", "4:", "1:", "2:"], "5:": ["5:", "2:"]}),
            (["9_"], {"9:": ["9:", "6:", "2:"]}),
        ],
    )
    @patch("scripts.compute_mappings.compute_tissue_and_cell_type_mappings.write_to_file")
    def test_write_ancestors_by_entity(self, write_to_file_mock, entities, expected_ancestors_dict):
        write_ancestors_by_entity(entities, self.get_graph(), "filename")
        assert write_to_file_mock.call_count == 1
        assert write_to_file_mock.call_args.args[1] == "filename"
        sorted_result_dict = {k: sorted(v) for k, v in write_to_file_mock.call_args.args[0].items()}
        sorted_expected_dict = {k: sorted(v) for k, v in expected_ancestors_dict.items()}
        assert sorted_result_dict == sorted_expected_dict

    @pytest.mark.parametrize(
        "entity_hierarchy, expected_descendants_dict",
        [
            (
                [["1_", "2_"], ["00_", "3_", "4_", "5_"], ["7_", "8_"]],
                {"1:": ["4:", "8:"], "2:": ["4:", "5:", "8:"], "3:": ["7:"], "4:": ["8:"]},
            ),
            ([["0_"], ["00_ (organoid)", "3_", "4_"], []], {"0:": ["3:", "00: (organoid)"]}),
        ],
    )
    @patch("scripts.compute_mappings.compute_tissue_and_cell_type_mappings.write_to_file")
    def test_write_descendants_by_entity(self, write_to_file_mock, entity_hierarchy, expected_descendants_dict):
        write_descendants_by_entity(entity_hierarchy, self.get_graph(), "filename")
        assert write_to_file_mock.call_count == 1
        assert write_to_file_mock.call_args.args[1] == "filename"
        sorted_result_dict = {k: sorted(v) for k, v in write_to_file_mock.call_args.args[0].items()}
        sorted_expected_dict = {k: sorted(v) for k, v in expected_descendants_dict.items()}
        assert sorted_result_dict == sorted_expected_dict
