import logging
logging.basicConfig(level=logging.INFO)
from obgraph.graph import Graph
import numpy as np


def test_from_dicts():
    g = Graph.from_dicts(
        {1: "ACTG", 2: "A", 3: "G", 4: "AAA"},
        {1: [2, 3],
         2: [4],
         3: [4]},
        [1, 2, 4]
    )

    assert g.get_node_size(1) == 4
    assert g.get_node_size(2) == 1
    assert g.get_node_size(3) == 1
    assert g.get_node_size(4) == 3

    assert g.get_node_sequence(2) == "A"

    assert list(g.get_edges(1)) == [2, 3]


    assert list(g.get_numeric_node_sequence(2)) == [0]

    print(g.get_numeric_node_sequence(np.array([1, 2, 3, 4])))
    assert list(g.get_numeric_node_sequences(np.array([1, 2, 3, 4]))) == [0, 1, 3, 2, 0, 2, 0, 0, 0]

    assert g.get_node_at_ref_offset(0) == 1
    assert g.get_node_at_ref_offset(3) == 1
    assert g.get_node_at_ref_offset(4) == 2
    assert g.get_node_at_ref_offset(5) == 4

    assert g.get_ref_offset_at_node(4) == 5
    assert g.get_ref_offset_at_node(2) == 4

    assert list(g.chromosome_start_nodes) == [1]


def test_sparse_graph():
    g = Graph.from_dicts(
        {1: "AGGG", 4: "CACCT"},
        {1: [4]},
        [1, 4]
    )
    assert list(g.get_edges(1)) == [4]
    assert g.get_node_at_ref_offset(4) == 4
    assert g.get_nodes_sequence([1, 4]) == "AGGGCACCT"


def test2():
    g = Graph.from_dicts(
        {0: "ACT", 1: "A", 2: "", 3: "GGG"},
        {0: [1, 2], 1: [3], 2: [3]},
        [0, 3]
    )

    assert list(g.chromosome_start_nodes) == [0]

    assert set(g.get_edges(0)) == set([1, 2])



def test_get_numeric_node_sequences_by_chromosome():
    g = Graph.from_dicts({
            1: "T", 2: "A", 1000: "A", 20: "ACAC", 23: "G", 50: "AAA", 51: "G"
        },
        {1: [2], 2: [1000], 1000: [20], 20: [23], 23: [50], 50: [51]},
        [1, 2, 1000, 20, 23, 50, 51],
        chromosome_start_nodes=[1, 20, 50]
    )

    sequences = list(g.get_numeric_node_sequences_by_chromosomes([1, 2, 1000, 20, 23, 50, 51]))

    assert np.all(sequences[0] == [3, 0, 0])
    assert np.all(sequences[1] == [0, 1, 0, 1, 2])
    assert np.all(sequences[2] == [0, 0, 0, 2])



test_from_dicts()
test_sparse_graph()
test2()
