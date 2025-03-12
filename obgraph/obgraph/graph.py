import json
import logging
import re
from collections import defaultdict

import numpy as np
from npstructures import RaggedArray, HashTable

from shared_memory_wrapper import to_file, from_file
from .mutable_graph import MutableGraph
from .nplist import NpList
from .util import encode_chromosome

import sys
sys.setrecursionlimit(20000)

class VariantNotFoundException(Exception):
    pass

def remap_array(array, from_values, to_values):
    index = np.digitize(array.ravel(), from_values, right=True)
    try:
        res = to_values[index].reshape(array.shape)
    except IndexError:
        logging.error("Invalid byte values in sequence. Max value: %d" % np.max(array))
        raise
    return res

def convert_sequence_array_to_numeric(sequence):
    assert type(sequence) == np.ndarray
    sequences_as_byte_values = sequence.astype("|S1").view(np.int8)
    # from byte values
    from_values = np.array([65, 67, 71, 84, 97, 99, 103, 116], dtype=np.uint64)  # NB: Must be increasing
    # to internal base values for a, c, t, g
    to_values = np.array([0, 1, 2, 3, 0, 1, 2, 3], dtype=np.uint64)
    return remap_array(sequences_as_byte_values, from_values, to_values).astype(np.uint8)

numeric_to_letter_sequence = np.array(["A", "C", "G", "T"])

class Graph:
    def __init__(self, nodes,
                 sequences: RaggedArray,
                 edges: RaggedArray,
                 node_to_ref_offset=None, ref_offset_to_node=None, chromosome_start_nodes=None, allele_frequencies=None, linear_ref_nodes_index=None,
                 numeric_node_sequences=None, linear_ref_nodes_and_dummy_nodes_index=None):
        self.nodes = nodes
        self.edges = edges
        self.sequences = sequences

        self.node_to_ref_offset = node_to_ref_offset
        self.ref_offset_to_node = ref_offset_to_node
        self._linear_ref_nodes_cache = None
        self.chromosome_start_nodes = chromosome_start_nodes
        self.allele_frequencies = allele_frequencies
        self.numeric_node_sequences = numeric_node_sequences

        if self.nodes is not None and linear_ref_nodes_index is None:
            logging.info("Makng a linear ref lookup, since it is not provided")
            self.linear_ref_nodes_index = np.zeros(len(self.nodes), dtype=np.uint8)
            self.linear_ref_nodes_index[np.array(list(self.linear_ref_nodes()))] = 1
            logging.info("Done making linear ref nodes index")
        else:
            self.linear_ref_nodes_index = linear_ref_nodes_index
        
        if self.nodes is not None and linear_ref_nodes_and_dummy_nodes_index is None:
            self.linear_ref_nodes_and_dummy_nodes_index = None
            logging.info("Making linear ref nodes and dummy nodes index")
            self.make_linear_ref_node_and_ref_dummy_node_index()
        else:
            self.linear_ref_nodes_and_dummy_nodes_index = linear_ref_nodes_and_dummy_nodes_index


    def __eq__(self, other):
        t = np.all(self.nodes == other.nodes) and np.all(self.edges == other.edges) and np.all(self.sequences == other.sequences)
        if self.chromosome_start_nodes is not None:
            if other.chromosome_start_nodes is None:
                return False

            t &= np.all(self.chromosome_start_nodes == other.chromosome_start_nodes)

        return t

    def node_has_edges(self, node, edges):
        e = self.get_edges(node)
        if len(edges) != len(e):
            return False

        for edge in edges:
            if edge not in edges:
                return False

        return True

    def get_all_nodes(self):
        node_to_n_edges = self.edges.shape.lengths   # a bit hacky, accessing length of RaggedShape in npstructures
        return np.union1d(np.where(self.nodes != 0)[0], np.where(node_to_n_edges > 0)[0])
        # This does not return nodes of size 0
        #return np.where(self.nodes != 0)[0]

    def __str__(self):
        return "Obgraph"

    def __repr__(self):
        return self.__str__()

    def to_mutable_graph(self):
        all_nodes = self.get_all_nodes()
        nodes = {node: self.get_node_size(node) for node in all_nodes}
        node_sequences = {node: self.get_node_sequence(node) for node in all_nodes}
        edges = {node: list(self.get_edges(node)) for node in all_nodes}
        return MutableGraph(nodes, node_sequences, edges, list(self.linear_ref_nodes()), self.node_to_ref_offset, self.ref_offset_to_node, self.chromosome_start_nodes, self.allele_frequencies)

    @classmethod
    def from_mutable_graph(cls, mutable_graph):
        return cls.from_dicts(mutable_graph.node_sequences, mutable_graph.edges, mutable_graph.linear_ref_nodes,
                              mutable_graph.chromosome_start_nodes)

    def get_node_size(self, node):
        return self.nodes[node]

    def get_node_allele_frequencies(self, nodes):
        if self.allele_frequencies is None:
            return np.ones(len(nodes))
        return self.allele_frequencies[nodes]

    def get_node_subsequence(self, node, start, end):
        return self.get_node_sequence(node)[int(start):int(end)]

    def get_node_sequence(self, node):
        return ''.join(numeric_to_letter_sequence[self.sequences[node]])

    def max_node_id(self):
        return len(self.nodes)-1

    def get_numeric_node_sequence(self, node):
        return self.sequences[node]

        if self.get_node_size(node) == 0:
            return np.array([-1])

        index_pos = self.node_to_sequence_index[node]
        return self.numeric_node_sequences[index_pos:index_pos+self.get_node_size(node)]

    def get_numeric_base_sequence(self, node, offset):
        if self.get_node_size(node) == 0:
            return -1

        return self.get_numeric_node_sequence(int(node))[int(offset)]

    def get_numeric_node_sequences_by_chromosomes(self, nodes):
        assert nodes[0] == self.chromosome_start_nodes[0], "First node must be first chromosome node"

        # yields one sequence per chromosome, assumes the nodes are sorted so that
        # nodes on same chromosome are together
        split_positions = list(np.where(np.in1d(nodes, self.chromosome_start_nodes))[0])
        #split_positions = list(np.searchsorted(nodes, self.chromosome_start_nodes))
        split_positions.append(nodes[-1])

        for start, end in zip(split_positions[0:-1], split_positions[1:]):
            yield self.get_numeric_node_sequences(nodes[start:end])

    def get_numeric_node_sequences(self, nodes):
        return self.sequences[nodes].ravel()

    def get_node_offset_at_chromosome_and_chromosome_offset(self, chromosome, offset):
        chromosome_position = chromosome - 1
        try:
            chromosome_start_node = self.chromosome_start_nodes[chromosome_position]
        except IndexError:
            chromosome_start_node = self.chromosome_start_nodes[0]

        chromosome_offset = self.get_ref_offset_at_node(chromosome_start_node)
        real_offset = int(chromosome_offset + offset)
        node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, offset)
        node_offset = self.get_ref_offset_at_node(node)
        return real_offset - node_offset

    def get_node_at_chromosome_and_chromosome_offset(self, chromosome, offset):
        chromosome_position = chromosome - 1
        try:
            chromosome_start_node = self.chromosome_start_nodes[chromosome_position]
        except IndexError:
            if len(self.chromosome_start_nodes) == 1:
                #logging.info("Assuming chromosome start node is the only start node %d" % self.chromosome_start_nodes[0])
                chromosome_start_node = self.chromosome_start_nodes[0]
            else:
                logging.error(
                    "Could not find chromosome start position for chromosome %d. Chromosome start nodes are %s" % (
                        chromosome, self.chromosome_start_nodes))
                raise


        # Shift offset with chromosome start offset
        chromosome_offset = self.get_ref_offset_at_node(chromosome_start_node)
        real_offset = int(chromosome_offset + offset)
        return self.get_node_at_ref_offset(real_offset)

    def get_edges(self, node):
        if node >= len(self.edges):
            return []

        return self.edges[node]

    def linear_ref_nodes(self):
        if self._linear_ref_nodes_cache is not None:
            return self._linear_ref_nodes_cache
        else:
            nodes = set(np.unique(self.ref_offset_to_node))
            self._linear_ref_nodes_cache = nodes
            return nodes

    def make_linear_ref_node_and_ref_dummy_node_index(self):
        linear_ref_nodes_and_dummy_nodes_index = np.zeros(len(self.linear_ref_nodes_index), dtype=np.uint8)
        for node in range(len(linear_ref_nodes_and_dummy_nodes_index)):
            if node % 100000 == 0:
                logging.info("%d nodes processed" % node)

            if self.is_linear_ref_node_or_linear_ref_dummy_node(node):
                linear_ref_nodes_and_dummy_nodes_index[node] = 1

        self.linear_ref_nodes_and_dummy_nodes_index = linear_ref_nodes_and_dummy_nodes_index

    def is_linear_ref_node_or_linear_ref_dummy_node(self, node):
        if self.linear_ref_nodes_and_dummy_nodes_index is not None:
            if self.linear_ref_nodes_and_dummy_nodes_index[node] == 1:
                return True
            else:
                return False

        if self.is_linear_ref_node(node):
            return True

        if self.get_node_size(node) == 0:
            # is a linear if it has a linear ref node as next node and the previous node on the linear ref from this linear ref node has an edge to node
            if len([next_node for next_node in self.get_edges(node) if self.is_linear_ref_node(next_node) and node in self.get_edges(self.get_node_at_ref_offset(int(self.get_ref_offset_at_node(next_node)-1)))]) > 0:
                return True

        return False

    def is_linear_ref_node(self, node):
        if self.linear_ref_nodes_index[node] != 0:
            return True
        return False

    def linear_ref_length(self):
        return len(self.ref_offset_to_node)

    def get_node_at_ref_offset(self, ref_offset):
        try:
            return self.ref_offset_to_node[ref_offset]
        except IndexError:
            logging.error("Invalid ref offset %s" % str(ref_offset))
            raise

    def get_ref_offset_at_node(self, node):
        return self.node_to_ref_offset[node]

    def get_chromosome_ref_offset_at_node(self, chromosome, node):
        chromosome_position = chromosome - 1
        try:
            chromosome_start_node = self.chromosome_start_nodes[chromosome_position]
        except IndexError:
            if len(self.chromosome_start_nodes) == 1:
                chromosome_start_node = self.chromosome_start_nodes[0]
            else:
                logging.error(
                "Could not find chromosome start position for chromosome %d. Chromosome start nodes are %s" % (
                    chromosome, self.chromosome_start_nodes))
                raise
        chromosome_offset = self.get_ref_offset_at_node(chromosome_start_node)
        return self.get_ref_offset_at_node(node) - chromosome_offset

    def get_node_offset_at_ref_offset(self, ref_offset):
        node = self.get_node_at_ref_offset(ref_offset)
        offset_at_node = self.get_ref_offset_at_node(node)
        return ref_offset - offset_at_node

    def to_file(self, file_name):
        return to_file(self, file_name)

    @classmethod
    def from_file(cls, file_name):
        return from_file(file_name)

    @classmethod
    def from_dicts(cls, node_sequences, edges, linear_ref_nodes, chromosome_start_nodes=None):
        assert linear_ref_nodes is not None
        logging.info("Making graph from dicts")
        nodes = np.sort(list(node_sequences.keys())).astype(np.uint32)
        node_ids = nodes
        max_node = nodes[-1]
        logging.info("Max node id is %d" % max_node)

        logging.info("Making sequences")
        node_sequences = [node_sequences[node] for node in nodes]
        node_sizes = [len(seq) for seq in node_sequences]
        node_sizes_array = np.zeros(nodes[-1]+1, dtype=np.uint32)
        node_sizes_array[nodes] = node_sizes
        nodes = node_sizes_array

        #from_nodes = []

        to_nodes = []
        n_edges = []

        to_nodes = np.concatenate([edges[n] for n in range(max_node+1) if n in edges]).astype(np.uint32)
        n_edges = np.array([len(edges[n]) if n in edges else 0 for n in range(max_node+1)], dtype=np.uint8)

        to_nodes_set = set(to_nodes)
        logging.info("Done preparing data from dicts")
        logging.info("There are %d node sequences" % len(node_sequences))

        # sequences
        sequence = np.array(list(''.join(node_sequences)))
        numeric_sequence = convert_sequence_array_to_numeric(sequence)
        sequences = RaggedArray(numeric_sequence, nodes, dtype=np.uint8)

        # linear ref index
        node_to_ref_offset = np.zeros(max_node + 1, np.uint64)
        ref_node_sizes = nodes[linear_ref_nodes]

        # print("Ref node sizes: %s"  % ref_node_sizes)
        ref_offsets = np.cumsum(ref_node_sizes)
        node_to_ref_offset[linear_ref_nodes[1:]] = ref_offsets[:-1]


        # Find last node to add to linear ref size
        last_node_in_graph = np.argmax(node_to_ref_offset)
        last_node_size = nodes[last_node_in_graph]
        logging.info("Last node in graph is %d with size %d" % (last_node_in_graph, last_node_size))

        ref_offset_to_node = np.zeros(int(np.max(node_to_ref_offset)) + last_node_size)
        index_positions = np.cumsum(ref_node_sizes)[:-1]
        # logging.info("INdex positions: %s" % index_positions)
        # logging.info("Linear ref nodes: %s" % linear_ref_nodes)
        # logging.info("Diff linear ref nodes: %s" % np.diff(linear_ref_nodes))
        ref_offset_to_node[index_positions] = np.diff(linear_ref_nodes)
        ref_offset_to_node[0] = linear_ref_nodes[0]
        # logging.info("Ref offset to node 1: %s" % ref_offset_to_node)
        ref_offset_to_node = np.cumsum(ref_offset_to_node, dtype=np.uint32)

        if chromosome_start_nodes is None:
            chromosome_start_nodes = [node for node in node_ids if node not in to_nodes_set]
            logging.info("Found chromosome start nodes to be %s" % chromosome_start_nodes)
        else:
            logging.info("Chromosome start nodes already set to %s" % chromosome_start_nodes)

        return Graph(node_sizes_array, sequences, RaggedArray(to_nodes, n_edges, dtype=np.uint32), node_to_ref_offset, ref_offset_to_node, chromosome_start_nodes)


    def get_snp_nodes(self, ref_offset, reference_bases, variant_bases, chromosome=1):
        node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset)
        node_offset = self.get_node_offset_at_chromosome_and_chromosome_offset(chromosome, ref_offset)

        if node_offset != 0:
            raise VariantNotFoundException("Node offset is %d at ref_offset %d. Variant bases: %s. Chromosome %d" % (node_offset, ref_offset, variant_bases, chromosome))

        prev_node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset - 1)

        _, possible_snp_nodes = self.find_nodes_from_node_that_matches_sequence(prev_node, variant_bases, [], [])

        if len(possible_snp_nodes) == 0:
            raise VariantNotFoundException("Did not find any possible snp nodes for variant at ref offset %s with variant sequence %s" % (ref_offset, variant_bases))

        next_ref_pos = ref_offset + len(reference_bases)
        next_ref_node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, next_ref_pos)

        for possible_snp_node in possible_snp_nodes:
            # Not true if deletion before snp
            #assert len(possible_snp_node) == 1, "SNP nodes should only be one node"
            if len([n for n in possible_snp_node if self.get_node_size(n) > 0]) != 1:
                logging.error("There should only be one non-empty SNP node")
                logging.info("Found nodes: %s" % (possible_snp_node))
                logging.info("prev node is %s, variant bases is %s" % (prev_node, variant_bases))
                raise Exception("error")

            potential_next = possible_snp_node[-1]  # Get last node, this will be the snp node if there are more nodes
            assert self.get_node_size(potential_next) > 0
            # Also require that this node shares next node with our ref snp node, if not this could be a false match against an indel node
            # this assumption does not work for subsitutions
            # why does this not work for substitutions??
            next_from_snp_node = self.get_edges(node)
            if len([n for n in next_from_snp_node if n in self.get_edges(potential_next)]) > 0:
                return node, potential_next

            # logging.info("Node after deletion: %d" % next_ref_node)
            if self.get_node_offset_at_chromosome_and_chromosome_offset(chromosome, next_ref_pos) != 0:
                raise VariantNotFoundException("Could not find next ref node")

            #logging.info("Next ref node: %d. Edges from potential: %s" % (next_ref_node, str(self.get_edges(potential_next))))
            if next_ref_node in self.get_edges(potential_next):
                return node, potential_next

        logging.error("Could not parse substitution at offset %d with bases %s. Next ref pos is %d. Next ref node is %d" % (ref_offset, variant_bases, next_ref_pos, next_ref_node))
        logging.info("Reference node for snp/subsitution is %d" % node)
        logging.error("Possible nodes are: %s" % possible_snp_nodes)
        raise VariantNotFoundException("Error")

    def get_deletion_nodes(self, ref_offset, deletion_length, deleted_sequence, chromosome=1):
        ref_offset += 1
        # Node is the reference node that is deleted (the first one if there are multiple)
        node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset)
        node_offset = self.get_node_offset_at_chromosome_and_chromosome_offset(chromosome, ref_offset)
        #logging.info("Processing deletion at ref pos %d with size %d. Node inside deletion: %d" % (ref_offset, deletion_length, node))

        if node_offset != 0:
            raise VariantNotFoundException("Node offset is %d, not 0" %  (node_offset))

        prev_node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset - 1)
        #logging.info("Node before deletion: %d" % prev_node)

        # Find next reference node with offset corresponding to the number of deleted base pairs
        next_ref_pos = ref_offset + deletion_length
        next_ref_node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset + deletion_length)
        #logging.info("Node after deletion: %d" % next_ref_node)
        if self.get_node_offset_at_chromosome_and_chromosome_offset(chromosome, next_ref_pos) != 0:
            #logging.error("Offset %d is not at beginning of node" % next_ref_pos)
            #logging.error("Node at %d: %d" % (next_ref_pos, next_ref_node))
            #logging.error("Ref length in deletion: %s" % deletion_length)
            #logging.info("Ref pos beginning of deletion: %d" % ref_offset)

            raise VariantNotFoundException("Deletion at ref offset %d not in graph" % ref_offset)

        # Find an empty node between prev_node and next_ref_node
        # This solution should handle cases where the reference has multiple nodes inside deleted segment (e.g. a SNP inside a deleted sequence)
        # This is not true anymore, since a deltion can be followed by an insertion, leading to two or more dummy nodes before reaching the linear ref again
        #deletion_nodes = [node for node in self.get_edges(prev_node) if next_ref_node in self.get_edges(node) and self.get_node_size(node) == 0]
        deletion_nodes = [node for node in self.get_edges(prev_node) if self.get_node_size(node) == 0]
        #debug = [(node, self.get_edges(node), self.get_node_size(node)) for node in self.get_edges(prev_node)]
        #logging.info("%s" % debug)




        if len(deletion_nodes) != 1:

            # There can be more than one deletion node if there is an insertion right after a deletion
            if len(deletion_nodes) > 1:
                # If there are more than 1 deletion node, we find all paths that matches the deleted sequence, and try to find the deletion node that starts before and ends after this path
                one_path, possible_deleted_paths = self.find_nodes_from_node_that_matches_sequence(prev_node, deleted_sequence, [], [])
                assert len(one_path) >  0 and len(possible_deleted_paths) > 0, "Found no paths from node %d that matches the deleted sequence %s" % (prev_node, deleted_sequence)
                deletion_nodes_orig = deletion_nodes.copy()

                #print(possible_deleted_paths)
                deletion_nodes = [n for n in deletion_nodes if len([path for path in possible_deleted_paths if len([edge_from_deletion_node for edge_from_deletion_node in self.get_edges(n) if edge_from_deletion_node in self.get_edges(path[-1])]) > 0]) > 0]

                if len(deletion_nodes) == 0:
                    logging.warning("Error on deletion at pos %s. Node: %d. Deletion nodes: %s. Orig deletion nodes: %s. Edges going out: %s" % (ref_offset, node, deletion_nodes, deletion_nodes_orig, self.get_edges(node)))
                    logging.error("Could not find any nodes among possible deletion nodes %s that wrapped around any of the paths %s" % (deletion_nodes_orig, possible_deleted_paths))
                    raise VariantNotFoundException("Deletion not found")
                elif len(deletion_nodes) > 1:
                    logging.warning("Deletion at ref offset %s has multiple empty nodes going out from %d: %s. They all wrap around a deleted path in %s" % (ref_offset, prev_node, deletion_nodes, possible_deleted_paths))
                    if deletion_nodes[0] not in self.get_edges(node):
                        logging.warning("Error on deletion. Deletion nodes: %s. Edges going out: %s" % (deletion_nodes, self.get_edges(node)))
                        raise VariantNotFoundException("Deletion not found")
                    # Assert that there is an empty node between this deletion node and next ref node
                    # This assertion won't always hold, e.g. if there are two insertions after a deltion, so just skip it
                    #assert len([n for n in self.get_edges(deletion_nodes[0]) if next_ref_node in self.get_edges(n) and self.get_node_size(n) == 0]), "There is more than one deletion node, and the next node
            else:
                #logging.warning("Could not find deletion at ref offset %d in graph" % ref_offset)
                raise VariantNotFoundException("There should be only one deletion node between %d and %d for variant at ref pos %d. There are %d. Edges out from %d are: %s. Edges out from %d are %s" % (prev_node, next_ref_node, ref_offset, len(deletion_nodes), prev_node, self.get_edges(prev_node), node, self.get_edges(node)))

        deletion_node = deletion_nodes[0]
        # Not true anymore, e.g. if there is a snp right after a deletion
        #assert len(self.get_edges(deletion_node)) == 1, "Deletion node %d has not exactly 1 edge. Variant at ref offset %d. Edges are %s" % (deletion_node, ref_offset, self.get_edges(deletion_node))
        return (node, deletion_node)

    def get_insertion_nodes(self, variant, chromosome=1):
        ref_offset = variant.position
        #logging.info("Ref offset: %d" % ref_offset)
        # Find node right before insertion node
        node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset-1)
        #logging.info("Node at ref offset: %d" % node)
        node_offset = self.get_node_offset_at_chromosome_and_chromosome_offset(chromosome, ref_offset-1)
        node_size = self.get_node_size(node)
        insertion_length = len(variant.variant_sequence) - 1
        insertion_sequence = variant.variant_sequence[1:]  # First base is not included

        if node_offset != node_size - 1:
            raise VariantNotFoundException("Node offset %d is not at end of node %d which has size %d. Insertion not found in graph." % (node_offset, node, node_size))

        # Find out which next node matches the insertion
        # NB! There might be more than one node with the same sequence going out before the variant, since a SNP can have same sequence as the insertion
        next_ref_node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset)
        variant_node = None
        for potential_next in self.get_edges(node):
            #if potential_next in self.linear_ref_nodes():
            if self.is_linear_ref_node(potential_next):
                continue  # Next should not be in linear ref

            #print("Processing insertion at ref offset %d with base %s" % (ref_offset, base))
            #print("  Node %d with offset %d. Node size: %d" % (node, node_offset, self.blocks[node].length()))

            if self.get_node_size(potential_next) == 0:
                continue  # This is an empty deletion node

            #next_bases = self.get_node_sequence(potential_next)[0:insertion_length].upper()
            next_bases = self.get_node_subsequence(potential_next, 0, insertion_length).upper()
            if next_bases == insertion_sequence.upper():
                if next_ref_node in self.get_edges(potential_next):
                    variant_node = potential_next
                    break

        if variant_node is None:
            raise VariantNotFoundException("Could not find insertion node. No nodes going out from %d has sequence %s" % (node, variant.variant_sequence))

        # Find linear reference node (empty dummy node)
        assert next_ref_node in self.get_edges(variant_node), "Failed parsing insertion %s. Found %d as next ref node after variant node, but there is no edge from variantnode %d to %d" % (variant, next_ref_node, variant_node, next_ref_node)
        # If there are multiple insertions, one could find the correct dummy node by choosing the one that goes to next node with lowest ref pos
        dummy_nodes = [node for node in self.get_edges(node) if next_ref_node in self.get_edges(node) and self.get_node_size(node) == 0]
        if len(dummy_nodes) != 1:
            logging.error("Edges from %d are %s" % (node, self.get_edges(node)))
            logging.error("Error when parsing insertion %s. There are not exactly 1 insertion node between %d and %d. Nodes of length 0 between are %s" % (variant, node, next_ref_node, dummy_nodes))
            raise Exception("")

        insertion_node = dummy_nodes[0]
        assert next_ref_node in self.get_edges(insertion_node), "Failed parsing insertion %s. Found %d as next ref node after dummy node, but there is no edge from dumy node %d to %d" % (variant, next_ref_node, insertion_node, next_ref_node)

        return (insertion_node, variant_node)

    def get_variant_nodes(self, variant):
        if variant.type == "SNP":
            return self.get_snp_nodes(variant.position-1, variant.get_reference_sequence(), variant.variant_sequence, variant.chromosome)
        elif variant.type == "SUBSTITUTION":
            return self.get_snp_nodes(variant.position, variant.get_reference_sequence(), variant.get_variant_sequence(), variant.chromosome)
        elif variant.type == "DELETION":
            return self.get_deletion_nodes(variant.position-1, len(variant.ref_sequence)-1, variant.get_deleted_sequence(), variant.chromosome)
        elif variant.type == "INSERTION":
            return self.get_insertion_nodes(variant, variant.chromosome)

        raise Exception("Invalid variant %s. Has no type set." % variant)

    def set_allele_frequencies_from_variants(self, variants, use_chromosome=None):
        frequencies = np.zeros(len(self.nodes), dtype=np.float16) + 1  # Set all to 1.0 initially

        for i, variant in enumerate(variants):
            if use_chromosome is not None:
                variant.chromosome = encode_chromosome(use_chromosome)

            if i % 1000000 == 0:
                logging.info("%d variants processed" % i)
            try:
                reference_node, variant_node = self.get_variant_nodes(variant)
                reference_af = variant.get_reference_allele_frequency()
            except VariantNotFoundException:
                logging.info("Variant %s not found in graph" % variant)
                continue
            variant_af = variant.get_variant_allele_frequency()
            frequencies[reference_node] = reference_af
            frequencies[variant_node] = variant_af

        self.allele_frequencies = frequencies

    def get_variant_nodes_in_region(self, chromosome, linear_ref_start, linear_ref_end):
        nodes = set(self.get_node_at_chromosome_and_chromosome_offset(chromosome, pos) for pos in range(linear_ref_start, linear_ref_end))
        return [self.get_edges(node) for node in nodes if len(self.get_edges(node)) == 2]

    def get_first_node(self):
        # chromosome start nodes should always be set
        #return self.chromosome_start_nodes[0]

        for candidate in [0, 1]:
            if len(self.get_edges(candidate)) > 0:
                return candidate

        raise Exception("Couldn't find first graph node")

    def convert_chromosome_ref_offset_to_graph_ref_offset(self, chromosome_offset, chromosome):
        # Add the graph ref offset at this chromosome to the chromosome offset
        chromosome_position = 0 # chromosome - 1
        return int(self.node_to_ref_offset[self.chromosome_start_nodes[chromosome_position]] + chromosome_offset)

    def find_nodes_from_node_that_matches_sequence(self, from_node, sequence, nodes_found, all_paths_found):
        # hack to fix issue with "n"
        sequence = sequence.replace("n", "a").replace("N", "A")

        #logging.info("== Searching from node %d with sequence %s. Variant type %s. Nodes found: %s. All paths found: %s" % (from_node, sequence, variant_type, nodes_found, all_paths_found))
        if sequence == "":
            # All of the sequence is used successfully, return
            all_paths_found.append(nodes_found)
            #logging.info("   RETURNING. All paths found: %s" % all_paths_found)
            return nodes_found, all_paths_found

        next_nodes = self.get_edges(from_node)
        result = (nodes_found, all_paths_found)
        for possible_next in next_nodes:

            #print("Checking next node %d" % possible_next)

            node_size = self.get_node_size(possible_next)
            if node_size == 0 or self.get_node_sequence(possible_next).lower() == sequence[0:node_size].lower():
                # This node is a match, we can continue searching
                new_sequence = sequence[node_size:]
                new_nodes_found = nodes_found.copy()
                new_nodes_found.append(possible_next)

                #print("   Matching sequence. New sequence is now %s" % new_sequence)
                result = self.find_nodes_from_node_that_matches_sequence(possible_next, new_sequence, new_nodes_found, all_paths_found)
                if not result:
                    continue
            else:
                #print("   No sequence match")
                continue

        return result

