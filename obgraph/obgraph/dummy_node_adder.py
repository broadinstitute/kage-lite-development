# Modifies a graph so that deletions/insertions has parallel dummy nodes
import logging
from itertools import product

from .graph import Graph, VariantNotFoundException


class DummyNodeAdder:
    def __init__(self, graph,  variants):
        self.graph = graph
        self.variants = variants
        self.mutable_graph = None
        self.current_new_node_id = len(self.graph.nodes) + 1
        self._n_variants_failed = 0
        self._n_variants_fixed = 0
        self._old_edges_to_new_node_mapping = {}

    def create_new_graph_with_dummy_nodes(self, use_mutable_graph=None):
        if use_mutable_graph is None:
            logging.info("Creating a mutable graph that can be changed")
            self.mutable_graph = self.graph.to_mutable_graph()
        else:
            self.mutable_graph = use_mutable_graph

        logging.info("Adding dummy nodes")
        for i, variant in enumerate(self.variants):
            if i % 10000 == 0:
                logging.info("%d variants processed" % i)

            if variant.type == "SNP" or variant.type == "SUBSTITUTION":
                continue

            try:
                self._add_dummy_edges_around_indel(variant)
            except VariantNotFoundException as e:
                logging.error("Could not find variant: %s" % str(e))
                self._n_variants_failed += 1
                continue
            self._n_variants_fixed += 1

        logging.info("Creating a new immutable graph from the mutable graph")
        #print(self.mutable_graph.get_all_nodes())
        #print(self.mutable_graph.node_sequences)
        logging.info("%d variants were not found in graph" % self._n_variants_failed)
        logging.info("%d variants were found and dummy nodes were added for these" % self._n_variants_fixed)
        return Graph.from_mutable_graph(self.mutable_graph)

    def get_edge_mapping(self):
        return self._old_edges_to_new_node_mapping

    def get_nodes_for_inserted_sequence_at_ref_pos(self, variant):
        inserted_sequence = variant.get_inserted_sequence()
        node_before_inserted_nodes = self.graph.get_node_at_chromosome_and_chromosome_offset(variant.chromosome, variant.position-1)
        inserted_nodes = self.mutable_graph.find_nodes_from_node_that_matches_sequence(node_before_inserted_nodes, inserted_sequence, variant.type, [], [])
        if len(inserted_nodes) == 0:
            raise VariantNotFoundException("Could not find inserted nodes for sequence %s, variant %s. Node before is %d" % (inserted_sequence, variant, node_before_inserted_nodes))

        return inserted_nodes

    def _add_dummy_edges_around_indel(self, variant):
        inserted_sequence = variant.get_inserted_sequence()
        _, possible_inserted_node_paths = self.get_nodes_for_inserted_sequence_at_ref_pos(variant)

        #logging.info("===========")
        #logging.info("Variant %s. Inserted nodes: %s" % (variant, possible_inserted_node_paths))

        did_add = False

        # If there are more than one path, some of them may be invalid (just by chance same sequence)
        # We only accespt paths where there is an edge from before the path to after the path
        possible_inserted_node_paths_orig = possible_inserted_node_paths.copy()
        if len(possible_inserted_node_paths) > 1:
            correct_paths = []
            for possible_path in possible_inserted_node_paths:
                in_nodes = self.mutable_graph.get_nodes_before(possible_path[0])
                out_nodes = self.mutable_graph.get_edges(possible_path[-1])

                if len([(from_node, to_node) for from_node, to_node in product(in_nodes, out_nodes) if to_node in self.graph.get_edges(from_node)]) > 0:

                    if variant.type == "DELETION":
                        # If deletion, the last position in the last node on the path must be on the linear reference and match the length of the deleted sequence
                        deletion_end_pos = variant.position + len(inserted_sequence)
                        if deletion_end_pos != self.graph.get_chromosome_ref_offset_at_node(variant.chromosome, possible_path[-1]) + self.graph.get_node_size(possible_path[-1]):
                            logging.info("Ignoring deletion path %s because ref pos at end is not correct" % possible_path)
                            continue

                    correct_paths.append(possible_path)

            possible_inserted_node_paths = correct_paths

        if len(possible_inserted_node_paths) == 0:
            logging.error("Did not find any valid inserted node paths for variant %s" % variant)
            logging.error("Possible paths are: %s" % possible_inserted_node_paths)
            logging.error("Original paths are: %s" % possible_inserted_node_paths_orig)
            raise VariantNotFoundException("Could not parse variatn")

        for inserted_nodes in possible_inserted_node_paths:

            # Find edges going from any node that goes into the inserted nodes and that ends at any node going out from the end of the inserted nodes
            # These edges are bypassing the inserted nodes
            in_nodes = self.mutable_graph.get_nodes_before(inserted_nodes[0])
            out_nodes = self.mutable_graph.get_edges(inserted_nodes[-1])
            #logging.info(" Out nodes: %s" % out_nodes)
            # This rule cannot be used since there might be an insertion right after a deletion that makes to_node not to be in edges from from_node
            # bypass_edges = [(from_node, to_node) for from_node, to_node in product(in_nodes, out_nodes) if to_node in self.mutable_graph.get_edges(from_node)]

            bypass_edges = set([(from_node, to_node) for from_node, to_node in product(in_nodes, out_nodes)])
            if len(bypass_edges) == 0:
                #logging.warning("Found no edges bypassing the inserted nodes %s. Nodes before are %s, nodes after are %s" % (inserted_nodes, in_nodes, out_nodes))
                #logging.warning("Variant is %s" % variant)
                continue

            if 178112 in inserted_nodes:
                logging.info("Adding bypass edges: %s" % bypass_edges)
                logging.info("Inserted node paths are %s" % possible_inserted_node_paths)

            self._add_new_dummy_node_for_edges(bypass_edges)
            did_add = True

        if not did_add:
            logging.warning("Not able to process variant %s" % variant)
            logging.warning("Possible inserted nodes in variant are %s" % possible_inserted_node_paths)
            raise VariantNotFoundException("No paths worked")

    def _add_new_dummy_node_for_edges(self, edges):
        dummy_node = self.current_new_node_id
        self.mutable_graph.add_node(dummy_node)

        # Remove all edges
        for from_node, to_node in edges:
            if to_node in self.mutable_graph.get_edges(from_node):
                self.mutable_graph.remove_edge(from_node, to_node)
                self._old_edges_to_new_node_mapping[(from_node, to_node)] = dummy_node

        # For all unique from_nodes add an edge to the dummy node
        for from_node in set((from_node for from_node, to_node in edges)):
            self.mutable_graph.add_edge(from_node, dummy_node)

        # For all unique to_nodes, add an edge from dummy node to to node
        for to_node in set((to_node for from_node, to_node in edges)):
            #print("Adding node from dummy node %d to node %d" % (dummy_node, to_node))
            self.mutable_graph.add_edge(dummy_node, to_node)

        self.current_new_node_id += 1

