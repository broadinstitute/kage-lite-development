import logging
import numpy as np
cimport numpy as np
cimport cython


def fill_zeros_increasingly(np.int64_t[:] array):
    cdef int i = 0
    for i in range(1, len(array)):
        if array[i] == 0:
            array[i] = array[i-1] + 1


@cython.wraparound(False)
def traverse_graph_by_following_nodes(graph, np.uint8_t[:] follow_nodes, split_into_chromosomes=False):
    logging.info("Traversing with cython")

    #assert type(follow_nodes) == np.ndarray

    logging.info("N nodes in follow set: %d" % len(follow_nodes))

    cdef np.ndarray[np.uint32_t] nodes_found = np.zeros(len(graph.nodes), dtype=np.uint32)
    cdef int node_index = 0

    # accessing RaggedArray internal stuff to speed things up later
    cdef np.uint32_t[:] edges = graph.edges.ravel()
    cdef np.int64_t[:] node_to_n_edges = graph.edges._shape.lengths
    cdef np.int64_t[:] node_to_edge_index = graph.edges._shape.starts

    cdef np.uint8_t[:] linear_nodes_index = graph.linear_ref_nodes_and_dummy_nodes_index

    cdef np.uint32_t[:] next_nodes
    cdef int edge_i = 0
    cdef int next_node = -1
    cdef int j


    #linear_ref_nodes = graph.linear_ref_nodes()
    cdef int current_node

    #print(type(current_node))
    chromosome_index_positions = []

    for current_node in graph.chromosome_start_nodes:
        if split_into_chromosomes:
            chromosome_index_positions.append(node_index)

        while True:
            #if current_node in follow_nodes:
                #nodes_found.append(current_node)
            nodes_found[node_index] = current_node
            node_index += 1

            #next_nodes = graph.get_edges(current_node)

            edge_i = node_to_edge_index[current_node]
            next_nodes = edges[edge_i:edge_i+node_to_n_edges[current_node]]

            if len(next_nodes) == 0:
                break

            next_node = -1
            if len(next_nodes) == 1:
                next_node = next_nodes[0]
            else:
                for j in range(next_nodes.shape[0]):
                #for potential_next in next_nodes:
                    if follow_nodes[next_nodes[j]] == 1 or (next_node == -1 and linear_nodes_index[next_nodes[j]] == 1):
                        next_node = next_nodes[j]

                assert next_node != -1, "Could not find next node from node %d" % current_node

            current_node = next_node

    #if current_node in follow_nodes:
    #nodes_found.append(current_node)
    #nodes_found[node_index] = current_node
    #node_index += 1

    if split_into_chromosomes:
        return nodes_found[0:node_index], chromosome_index_positions

    return nodes_found[0:node_index]
