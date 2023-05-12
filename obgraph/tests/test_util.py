from obgraph.util import create_coordinate_map
from obgraph.graph import Graph
import numpy as np
from obgraph.util import phased_genotype_matrix_to_haplotype_matrix


def test_create_coordinate_map():

    g = Graph.from_dicts(
        {1: "AAAA", 2: "G", 3: "GGGGGGG", 4: "ACA", 5: "C"},
        {1: [2, 3], 2: [4], 3: [4], 4: [5]},
        [1, 2, 4, 5]
    )

    path = np.array([1, 3, 4, 5])

    m = create_coordinate_map(path, g, chromosome_index=0)

    assert m[0] == 0
    assert m[11] == 5
    assert m[14] == 8


def test_phased_genotype_matrix_to_haplotype_matrix():
    matrix = np.array([
        [0, 1, 2, 3],
        [3, 2, 1, 0]
    ])

    true = np.array([
        [0, 0, 0, 1, 1, 0, 1, 1],
        [1, 1, 1, 0, 0, 1, 0, 0]
    ])


    assert np.all(
        phased_genotype_matrix_to_haplotype_matrix(matrix) == true
    )

test_create_coordinate_map()

