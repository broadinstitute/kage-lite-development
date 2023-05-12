from obgraph.haplotype_nodes import HaplotypeNodes, NodeToHaplotypes, HaplotypeToNodes
import numpy as np

def test_node_to_haplotypes():
    nodes = np.zeros((2, 10))
    nodes[0,0:3] = np.array([1, 2, 3])
    nodes[1,0:3] = np.array([1, 2, 4])
    haplotype_nodes = HaplotypeNodes(nodes, None)

    lookup = NodeToHaplotypes.from_haplotype_nodes(haplotype_nodes)
    lookup.to_file("testfile")
    lookup2 = NodeToHaplotypes.from_file("testfile")

    assert list(lookup2.get_haplotypes_on_node(1)) == [0, 1]
    assert list(lookup2.get_haplotypes_on_node(2)) == [0, 1]
    assert list(lookup2.get_haplotypes_on_node(3)) == [0]
    assert list(lookup2.get_haplotypes_on_node(4)) == [1]

test_node_to_haplotypes()


def test_haplotype_to_nodes():
    haplotype_to_nodes = HaplotypeToNodes.from_flat_haplotypes_and_nodes(
        [0, 0, 0, 4, 4, 2, 2, 2, 2],
        [1, 2, 3, 1, 3, 1, 2, 3, 4]
    )


    haplotype_to_nodes.to_file("testfile")
    new = HaplotypeToNodes.from_file("testfile")

    assert list(new.get_nodes(0)) == [1, 2, 3]
    assert list(new.get_nodes(4)) == [1, 3]
    assert list(new.get_nodes(2)) == [1, 2, 3, 4]


test_haplotype_to_nodes()
