import logging


def traverse_graph_by_following_nodes(graph, follow_nodes, only_add_nodes_in_follow_nodes=False):
    assert type(follow_nodes) == set

    logging.info("N nodes in follow set: %d" % len(follow_nodes))

    nodes_found = []
    #linear_ref_nodes = graph.linear_ref_nodes()
    current_node = graph.get_first_node()
    logging.info("Starting traversing")
    chromosome_index = 0
    i = 0
    while True:
        if current_node in follow_nodes and only_add_nodes_in_follow_nodes:
            nodes_found.append(current_node)
        else:
            nodes_found.append(current_node)

        next_nodes = graph.get_edges(current_node)
        if len(next_nodes) == 0:
            # Check if we are at the end of a chromosome, if so, jump to next chromosome
            chromosome_index += 1
            if chromosome_index < len(graph.chromosome_start_nodes):
                next_nodes = [graph.chromosome_start_nodes[chromosome_index]]
                logging.info("Jumping to next chromosome with start node %d" % next_nodes[0])
            else:
                logging.info("Stopping at node %d. Did not find any next chromosome among nodes %s" % (current_node, graph.chromosome_start_nodes))
                break

        if i % 100000 == 0:
            logging.info("%i nodes traversed" % (i))

        next_node = None
        if len(next_nodes) == 1:
            next_node = next_nodes[0]
        else:
            for potential_next in next_nodes:
                if potential_next in follow_nodes or (next_node is None and graph.is_linear_ref_node(potential_next)):
                    next_node = potential_next

            # IF did not find any nodes, choose an empty node if there is one
            if next_node is None:
                next_node = [n for n in next_nodes if graph.get_node_size(n) == 0]
                assert len(next_node) == 1, "There are multiple empty nodes from node %d: %s" % (current_node, next_nodes)
                next_node = next_node[0]


        assert next_node is not None, "Could not find any new nodes from %d. Edges are %s" % (current_node, str(next_nodes))
        current_node = next_node
        i += 1

    if current_node in follow_nodes or not only_add_nodes_in_follow_nodes:
        nodes_found.append(current_node)

    return nodes_found
