import logging
import pickle

import numpy as np

from .graph import Graph
from .mutable_graph import MutableGraph


def convert_gfa_ids_to_numeric(gfa_file_name, out_base_name):
    gfa = open(gfa_file_name)
    out_gfa = open(out_base_name + ".gfa", "w")
    id_mapping = {}

    current_id = 0
    for line in gfa:
        l = line.strip().split()
        if l[0] == "S" or l[0] == "A":
            id = l[1]
            if id in id_mapping:
                new_id = str(id_mapping[id])
            else:
                new_id = str(current_id)
                id_mapping[id] = current_id
                current_id += 1

            l[1] = new_id
        elif l[0] == "L":
            assert l[1] in id_mapping and l[3] in id_mapping, "Gfa may have L lines before S lines"
            l[1] = str(id_mapping[l[1]])
            l[3] = str(id_mapping[l[3]])

        new_line = "\t".join(l)
        out_gfa.write(new_line + "\n")

    with open(out_base_name + ".id_mapping", "wb") as f:
        pickle.dump(id_mapping, f)
        logging.info("Wrote id mapping to %s" % out_base_name + ".id_mapping")


def create_graph_from_gfa_file(file_name):

    graph = MutableGraph()

    path_lines = []

    with open(file_name) as f:
        for line in f:
            l = line.split()
            if line.startswith("S"):
                id = int(l[1])
                sequence = l[2]
                graph.add_node(id, sequence)
            elif line.startswith("L"):
                from_node = int(l[1])
                to_node = int(l[3])

                if l[2] != "+" or l[4] != "+":
                    logging.warning("Only links from positive side to positive side are supported. Found link with negative side")

                if l[5] != "*":
                    logging.warning("Overlaps between segments not supported. Segments %d-%d will be linked by edge"
                                    % (from_node, to_node))

                graph.add_edge(from_node, to_node)
            elif line.startswith("P"):
                if not l[1].startswith("_alt"):
                    path_lines.append(l)

    #assert len(path_lines) == 1, "Supports only one path, which should be reference path. Found %d" % len(path_lines)

    # find reference node from path
    linear_ref_nodes = []
    chromosome_start_nodes = []
    for path_line in path_lines:
        path_name = path_line[1]
        logging.info("Processing path %s" % path_name)
        nodes = path_line[2].split(",")
        nodes = [int(node.replace("+", "")) for node in nodes]
        chromosome_start_nodes.append(nodes[0])
        linear_ref_nodes.extend(nodes)
        logging.info("There are %d linear ref nodes in path %s"  % (len(nodes), path_name))

    graph.set_linear_ref_nodes(linear_ref_nodes)
    graph.chromosome_start_nodes = np.array(chromosome_start_nodes)

    logging.info("Chromosome start nodes: %s" % graph.chromosome_start_nodes)
    return Graph.from_mutable_graph(graph)


