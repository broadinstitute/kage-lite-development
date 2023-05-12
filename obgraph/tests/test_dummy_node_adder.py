from obgraph.graph import Graph, DummyNodeAdder
from obgraph.graph_construction import GraphConstructor
from obgraph.variants import VcfVariants, VcfVariant

def test_simple_insertion():
    graph = Graph.from_dicts(
        {1: "ACTG", 2: "C", 3: "AAAA"},
        {1: [2, 3], 2: [3]},
        [1, 3]
    )

    variants = VcfVariants([VcfVariant(1, 4, "G", "GC", type="INSERTION")])
    dummy_node_adder = DummyNodeAdder(graph, variants)
    new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes()

    assert new_graph.node_has_edges(5, [3])
    assert new_graph.node_has_edges(1, [2, 5])
    assert new_graph.node_has_edges(2, [3])
    #print(new_graph)


def test_double_deletion_with_snp_inside_first_deletion():

    graph = Graph.from_dicts(
        {1: "ACTG", 2: "A", 3: "C", 4: "T", 5: "AAA", 6: "G"},
        {
            1: [2, 5, 6],
            2: [3, 4],
            3: [5, 6],
            4: [5, 6],
            5: [6]
        },
        [1, 2, 4, 6]
    )

    variants = VcfVariants([VcfVariant(1, 4, "GAT", "G", type="DELETION"), VcfVariant(1, 6, "TAAA", "T", type="DELETION")])
    dummy_node_adder = DummyNodeAdder(graph, variants)
    new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes()
    print(new_graph)

def test_double_deletion_with_snp_inside_first_deletiod_and_false_deletion_path():

    repeated_sequence = "AGGTCCCAGGTCCATCT"
    graph = Graph.from_dicts(
        {1: "TTTT", 2: "AGGTCC", 3: "C", 4: "A", 5: repeated_sequence, 6: repeated_sequence},
        {
            1: [2, 5, 6],
            2: [3, 4],
            3: [5, 6],
            4: [5, 6],
            5: [6]
        },
        [1, 2, 3, 5, 6]
    )

    variants = VcfVariants([VcfVariant(1, 4, "TAGGTCCC", "T", type="DELETION"), VcfVariant(1, 11, "CAGGTCCCAGGTCCATCT", "C", type="DELETION")])
    dummy_node_adder = DummyNodeAdder(graph, variants)
    new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes()

    print(new_graph)
    assert list(new_graph.get_edges(1)) == [2, 8]
    assert list(new_graph.get_edges(2)) == [3, 4]
    assert list(new_graph.get_edges(3)) == [5, 9]
    assert list(new_graph.get_edges(4)) == [5, 9]
    assert list(new_graph.get_edges(9)) == [6]
    assert list(new_graph.get_edges(8)) == [5, 9]

def test_insertion_with_multiple_paths():

    graph = Graph.from_dicts(
        {1: "AAAG", 2: "GAGT", 3: "GA", 4: "C", 5: "G", 6: "T"},
        {
            1: [2, 3],
            2: [3],
            3: [4, 5],
            4: [6],
            5: [6]
        },
        [1, 3, 5, 6]
    )

    variants = VcfVariants([VcfVariant(1, 4, "G", "GGAGT", type="INSERTION")])
    dummy_node_adder = DummyNodeAdder(graph, variants)
    new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes()
    assert list(new_graph.get_edges(1)) == [2, 8]
    assert list(new_graph.get_edges(8)) == [3]
    print(new_graph)


def test_tricky_case_nested_deletions():
    graph = Graph.from_dicts(
        {1: "TATAT", 2: "AT", 3: "A", 4: "T", 5: "A", 6: "A", 7: "T", 8: "A", 9: "GG" },
        {
            1: [2, 6],
            2: [3, 6],
            3: [4, 5],
            4: [6],
            5: [6],
            6: [7, 8],
            7: [9],
            8: [9]
        },
        [1, 2, 3, 5, 6, 8, 9]
    )

    variants = VcfVariants(
        [
            VcfVariant(1, 5, "TATAA", "T", type="DELETION"),
            VcfVariant(1, 7, "TAA", "T", type="DELETION"),
            VcfVariant(1, 5, "A", "T", type="SNP"),
        ]
    )

    dummy_node_adder = DummyNodeAdder(graph, variants)
    new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes()
    print(new_graph)

    assert list(new_graph.get_edges(1)) == [2, 11]
    assert list(new_graph.get_edges(2)) == [3, 12]
    assert list(new_graph.get_edges(11)) == [6]
    assert list(new_graph.get_edges(12)) == [6]



def test_overlapping_deletions():
    graph = Graph.from_dicts(
        {1: "AA", 2: "TCTG", 3: "TCT", 4: "G", 5: "A", 6: "GG"},
        {
            1: [2, 3],
            2: [3, 6],
            3: [4, 5],
            4: [6],
            5: [6]
        },
        [1, 2, 3, 5, 6]
    )

    variants = VcfVariants(
        [
            VcfVariant(1, 2, "ATCTG", "A", type="DELETION"),
            VcfVariant(1, 6, "GTCTA", "T", type="DELETION"),
            VcfVariant(1, 10, "A", "G", type="SNP")
        ]
    )
    dummy_node_adder = DummyNodeAdder(graph, variants)
    new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes()

    assert list(new_graph.get_edges(1)) == [2, 8]
    assert list(new_graph.get_edges(8)) == [3, 9]
    assert list(new_graph.get_edges(2)) == [3, 9]
    assert list(new_graph.get_edges(9)) == [6]

    ref_node, var_node = new_graph.get_variant_nodes(variants[1])
    assert ref_node == 3
    assert var_node == 9
    print(new_graph)



def test_insertion_with_identical_false_path():
    graph = Graph.from_dicts(
        {1: "AA", 2: "TCTG", 3: "TCTG", 4: "GG"},
        {
            1: [2, 3],
            2: [3],
            3: [4],
        },
        [1, 3, 4]
    )

    variants = VcfVariants(
        [
            VcfVariant(1, 2, "A", "ATCTG", type="INSERTION"),
        ]
    )
    dummy_node_adder = DummyNodeAdder(graph, variants)
    new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes()
    print(new_graph)

    assert list(new_graph.get_edges(1)) == [2, 6]
    assert list(new_graph.get_edges(6)) == [3]
    assert list(new_graph.get_edges(2)) == [3]

    ref_node, var_node = new_graph.get_variant_nodes(variants[0])
    assert ref_node == 6
    assert var_node == 2


def test_full_create_with_two_insertions_at_same_position():
    variants = VcfVariants(
        [
            VcfVariant(1, 3, "A", "ATATT", type="INSERTION"),
            VcfVariant(1, 3, "A", "ATATTTATT", type="INSERTION")
        ]
    )

    ref_sequence = "ATAGGGGGG"
    constructor = GraphConstructor(ref_sequence, variants)
    without_dummy = constructor._graph.to_mutable_graph()
    print(without_dummy)
    print(without_dummy.nodes)

    assert set(without_dummy.edges[1]) == set([2, 3, 4])
    assert set(without_dummy.edges[2]) == set([4])
    assert set(without_dummy.edges[3]) == set([4])

    graph = constructor.get_graph_with_dummy_nodes()
    assert set(graph.edges[1]) == set([2, 3, 6, 7])
    assert set(graph.edges[2]) == set([4])
    assert set(graph.edges[6]) == set([4])


