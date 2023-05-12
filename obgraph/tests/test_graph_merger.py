from obgraph.graph_merger import merge_graphs
from obgraph.graph import Graph


def test_simple():

    graph1 = Graph.from_dicts(
        {1: "ACTG", 2: "A", 3: "C", 4: "ACT"},
        {1: [2, 3], 2: [4], 3: [4]},
        [1, 2, 4]
    )

    graph2 = Graph.from_dicts(
        {1: "AAAA", 2: "A", 3: "C", 4: "ACT"},
        {1: [2, 3], 2: [4], 3: [4]},
        [1, 2, 4]
    )

    merged_graph = merge_graphs([graph1, graph2])

    print(merged_graph)
    print(merged_graph.edges)
    assert list(merged_graph.get_edges(1)) == [2, 3]
    assert list(merged_graph.get_edges(6)) == [7, 8]
    assert merged_graph.get_node_sequence(6) == "AAAA"
    assert merged_graph.get_node_sequence(7) == "A"
    assert merged_graph.get_node_sequence(8) == "C"

    assert merged_graph.get_node_at_ref_offset(0) == 1
    assert merged_graph.get_node_sequence(merged_graph.get_node_at_ref_offset(8)) == "AAAA"
    assert merged_graph.get_node_sequence(merged_graph.get_node_at_ref_offset(4)) == "A"
    assert len(merged_graph.get_edges(merged_graph.get_node_at_ref_offset(8))) == 2
    assert len(merged_graph.get_edges(merged_graph.get_node_at_ref_offset(11))) == 2

    assert merged_graph.get_ref_offset_at_node(6) == 8
    assert 7 in merged_graph.linear_ref_nodes()
    assert merged_graph.get_ref_offset_at_node(7) == 12

    assert list(merged_graph.chromosome_start_nodes) == [1, 6]

    merged_graph.to_file("merged_graph.npz")
    merged_graph2 = Graph.from_file("merged_graph.npz")

test_simple()