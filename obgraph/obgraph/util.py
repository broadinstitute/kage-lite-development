import logging
import resource

import numpy as np

from .cython_traversing import fill_zeros_increasingly


def encode_chromosome(chromosome):
    if chromosome.startswith("chr"):
        return encode_chromosome(chromosome.replace("chr", ""))

    if chromosome.upper() == "X":
        return 23
    elif chromosome.upper() == "Y":
        return 24
    else:
        return int(chromosome)


def phased_genotype_matrix_to_haplotype_matrix(genotype_matrix):
    haplotype_matrix = np.zeros((genotype_matrix.shape[0], genotype_matrix.shape[1]*2), dtype=genotype_matrix.dtype)
    haplotype_matrix[:,::2] = (genotype_matrix // 2) == 1
    haplotype_matrix[:,1::2] = (genotype_matrix % 2) == 1
    return haplotype_matrix


def get_number_of_variants_and_individuals_from_vcf(file_name):
    raise NotImplementedError("")
    file = open(file_name)
    for line in file:
        continue


def fill_zeros_with_last(arr, initial=0):
    ind = np.nonzero(arr)[0]
    cnt = np.cumsum(np.array(arr, dtype=bool))
    return np.where(cnt, arr[ind[cnt-1]], initial)


def create_coordinate_map(path_nodes, graph, chromosome_index=0):
    path_node_sizes = graph.nodes[path_nodes]
    path_offsets = np.concatenate([[0], np.cumsum(path_node_sizes)[:-1]]).astype(int)
    linear_ref_offsets = graph.node_to_ref_offset[path_nodes] - graph.get_ref_offset_at_node(graph.chromosome_start_nodes[chromosome_index])

    print(path_offsets)
    print(linear_ref_offsets)

    # lookup is from path_offsets to linear_ref_offsets
    lookup = np.zeros(int(path_offsets[-1])+1, dtype=np.int64)
    lookup[path_offsets] = linear_ref_offsets
    print(lookup)
    fill_zeros_increasingly(lookup)

    """
    # should fill zeros with 1 more than last, if not large nodes will not give accurate mapping
    lookup = fill_zeros_with_last(lookup)

    # Add node offsets
    path_node_to_dist = np.zeros(int(path_offsets[-1])+1)
    path_node_to_dist[path_offsets] = path_offsets
    path_node_to_dist = fill_zeros_with_last(path_node_to_dist)
    print("Path node to dist: %s" % path_node_to_dist)

    node_offsets = np.arange(0, len(path_node_to_dist)) - path_node_to_dist
    print("Lookup before adding: %s" % lookup)
    lookup += node_offsets

    """
    print("Final lookup: %s" % lookup)
    return lookup



def log_memory_usage_now(logplace=""):
    memory = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1000000
    logging.info("Memory usage (%s): %.4f GB" % (logplace, memory))


