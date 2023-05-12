import numpy as np
from obgraph.position_id import PositionId
from obgraph.graph import Graph


class DummyGraph:
    def __init__(self, nodes):
        self.nodes = nodes

def test():
    graph = DummyGraph(np.array([1, 3, 1, 0, 5, 4]))
    position_id = PositionId.from_graph(graph)

    assert position_id.get(0, 0) == 0
    assert position_id.get(1, 0) == 1
    assert position_id.get(1, 1) == 2
    assert position_id.get(1, 2) == 3
    assert position_id.get(2, 0) == 4


def test2():
    graph = Graph.from_dicts(
        {1: "CCCCC", 2: "G", 3: "", 4: "ACT", 5: "", 6: "GC", 7: "A", 8: "T", 9: "G", 10: "GGG"},
        {1: [2, 3], 2: [4], 3: [4], 4: [5, 6], 5: [7], 6: [7], 7: [8, 9], 8: [10], 9: [10]},
        [1, 2, 4, 7, 8, 10]
    )
    position_id = PositionId.from_graph(graph)
    ids = set()
    for node in graph.get_all_nodes():
        if graph.get_node_size(node):
            for offset in range(0, graph.get_node_size(node)):
                id = position_id.get(node, offset)
                assert id not in ids
                ids.add(id)

test()
test2()

