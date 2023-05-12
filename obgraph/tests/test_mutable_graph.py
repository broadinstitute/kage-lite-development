from obgraph import MutableGraph

def test_create():
    graph = MutableGraph({1: 4, 2: 3, 3: 1, 4: 1}, {1: "ACTG", 2: "A", 3: "C", 4: "AAAA"}, {1: [2, 3], 3: [4], 2: [4]}, [1, 2, 4])

    assert graph.get_node_size(1) == 4
    assert graph.get_edges(1) == [2, 3]
    assert graph.get_nodes_before(2) == [1]
    assert 2 in graph.get_nodes_before(4)
    assert 3 in graph.get_nodes_before(4)
    assert graph.get_node_sequence(2) == "A"


def test_get_nodes_matching_sequence_single_node():
    graph = MutableGraph({1: 4, 2: 3, 3: 1, 4: 1}, {1: "ACTG", 2: "A", 3: "C", 4: "AAAA"}, {1: [2, 3], 3: [4], 2: [4]}, )

    nodes, all_paths = graph.find_nodes_from_node_that_matches_sequence(1, "A", None, [], [])
    print(all_paths)
    assert nodes == [2]

def test_get_nodes_matching_sequence_double_deletion_and_snp():
    graph = MutableGraph(
        {1: 4, 2: 1, 3: 1, 4: 1, 5: 3, 6: 1},
        {1: "ACTG", 2: "A", 3: "C", 4: "T", 5: "AAA", 6: "G"},
        {
            1: [2, 5, 6],
            2: [3, 4],
            3: [5, 6],
            4: [5, 6],
            5: [6]
        }
    )

    nodes, all_paths = graph.find_nodes_from_node_that_matches_sequence(1, "AT", None, [], [])
    print(all_paths)
    assert nodes == [2, 4]

    nodes, all_paths = graph.find_nodes_from_node_that_matches_sequence(1, "AAA", None, [], [])
    print(all_paths)
    assert nodes == [5]

def test_get_nodes_matching_multiple_paths():
    graph = MutableGraph(
        {1: 1, 2: 1, 3: 1, 4: 1, 5: 2},
        {1: "A", 2: "G", 3: "G", 4: "A", 5: "GG"},
        {
            1: [2, 4, 5],
            2: [3],
            3: [5],
            4: [5],
        }
    )

    nodes, all_paths = graph.find_nodes_from_node_that_matches_sequence(1, "GG", "", [], [])
    print(all_paths)
    assert len(all_paths) == 2
    assert [2, 3] in all_paths
    assert [5] in all_paths


test_create()
test_get_nodes_matching_sequence_single_node()
test_get_nodes_matching_sequence_double_deletion_and_snp()
test_get_nodes_matching_multiple_paths()
