import logging
logging.basicConfig(level=logging.INFO)
from obgraph.gfa import create_graph_from_gfa_file


def test_graph_from_gfa():
    gfa_lines = [
        "S 1 ACT",
        "S 2 G",
        "S 3 GGGG",
        "S 4 AAAA",
        "L 1 + 2 + *",
        "L 1 + 3 + *",
        "L 2 + 4 + *",
        "L 3 + 4 + *",
        "P 1 1+,2+,4+ *"
    ]

    with open("gfa.tmp", "w") as f:
        f.writelines((l + "\n" for l in gfa_lines))

    graph = create_graph_from_gfa_file("gfa.tmp")

    assert list(graph.get_edges(1)) == [2, 3]
    assert list(graph.get_edges(2)) == [4]
    assert list(graph.get_edges(3)) == [4]
    assert list(graph.get_edges(4)) == []
    assert graph.get_node_sequence(3) == "GGGG"
    assert graph.get_node_sequence(4) == "AAAA"

    assert graph.is_linear_ref_node(1)
    assert not graph.is_linear_ref_node(3)
    assert graph.is_linear_ref_node(2)
    assert graph.is_linear_ref_node(4)


test_graph_from_gfa()