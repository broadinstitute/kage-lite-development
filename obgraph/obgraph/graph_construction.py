import logging
from collections import defaultdict

from .dummy_node_adder import DummyNodeAdder
from .graph import Graph
from .mutable_graph import MutableGraph
from .variants import VcfVariants


class EmptyNodeException(Exception):
    pass

class GraphConstructor:
    def __init__(self, reference_sequence, variants: VcfVariants):
        self.reference_sequence = reference_sequence
        self.variants = variants
        self.breakpoints = []
        self.make_sorted_breakpoints()

        self._ref_pos_to_node_after = defaultdict(list)
        self._node_to_ref_pos_after = defaultdict(list)
        self._node_to_ref_pos_before = defaultdict(list)
        self._ref_pos_to_node_start = defaultdict(list)

        self._node_ids = []
        self._node_sizes = []
        self._node_sequences = []
        self._reference_nodes = []
        self._current_node_id = 1
        self._edges_from = []
        self._edges_to = []

        self._edges_from_ref_pos_to_node = []
        self._deletions = defaultdict(list)  # A lookup from a ref pos to another ref pos representing a deltion
        self._edges_added = set()

        self._mutable_graph = MutableGraph()


        self.make_nodes()
        self.make_edges()
        self._graph = None
        self._graph_with_dummy_nodes = None
        logging.info("Making a graph")
        self._graph = self.get_graph()
        logging.info("Adding dummy nodes for indels")
        self.add_dummy_nodes()

    def get_graph(self):
        return Graph.from_mutable_graph(self._mutable_graph)

    def get_graph_with_dummy_nodes(self):
        return self._graph_with_dummy_nodes

    def add_dummy_nodes(self):
        dummy_node_adder = DummyNodeAdder(self._graph, self.variants)
        new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes(self._mutable_graph)
        self._graph_with_dummy_nodes = new_graph

    def _make_edge(self, from_node, to_node, ref_pos_before_to_node):
        if (from_node, to_node) in self._edges_added:
            return

        assert to_node != from_node, "Trying to make edge from same node to same node %d to %d" % (from_node, to_node)

        self._edges_from.append(from_node)
        self._edges_to.append(to_node)
        self._edges_added.add((from_node, to_node))
        self._mutable_graph.add_edge(from_node, to_node)

        # Check if there is a deletion going from start of to_node, if so we should also have an edge past to node
        # Check if ref_pos_after comes in at a deletion, if so also add edges to after deletion
        for next_ref_pos in self._deletions[ref_pos_before_to_node+1]:
            for node_after_deletion in self._ref_pos_to_node_after[next_ref_pos - 1]:
                #logging.info("    Adding deletion edge from %d to %d" % (from_node, node_after_deletion))
                self._make_edge(from_node, node_after_deletion, next_ref_pos-1)


    def _make_node(self, ref_position_before_node, ref_position_after_node, sequence, is_ref_node=False, is_deletion=False):
        if sequence == "":
            logging.error("Ref pos before node: %d" % ref_position_before_node)
            raise EmptyNodeException("Empty sequence for node")

        #size = len(sequence)
        #self._node_ids.append(self._current_node_id)
        #self._node_sequences.append(np.array(list(sequence)))
        #self._node_sizes.append(size)

        #if is_ref_node and len(sequence) > 0:  # We never want empty dummy nodes as reference node (by definition)
        #    self._reference_nodes.append(self._current_node_id)

        self._ref_pos_to_node_after[ref_position_before_node].append(self._current_node_id)
        self._ref_pos_to_node_start[ref_position_before_node+1].append(self._current_node_id)
        #logging.info("         Addding node %d to ref_pos_to_node_after for ref pos %d" % (self._current_node_id, ref_position_before_node))
        self._node_to_ref_pos_after[self._current_node_id].append(ref_position_after_node)
        self._node_to_ref_pos_before[self._current_node_id].append(ref_position_before_node)
        node_id = self._current_node_id

        self._mutable_graph.add_node(node_id, sequence, is_ref_node=is_ref_node)

        self._current_node_id += 1
        return node_id

    def make_nodes(self):
        logging.info("Making nodes")
        prev_ref_node_end = -1

        i = 0
        for breakpoint_position, variant in self.breakpoints:
            i += 1

            if i % 100000 == 0:
                logging.info("%d/%d breakpoints processed" % (i, len(self.breakpoints)))

            #logging.info("At breakpoint %s, %s" % (breakpoint_position, variant))
            # Always make a ref variant from prev_ref_node_end to this breakpoints
            if breakpoint_position > prev_ref_node_end:
                #logging.info("   Making a reference node %d. Pos before/after: %d/%d." % (
                #self._current_node_id, prev_ref_node_end, breakpoint_position + 1))
                try:
                    sequence_between_breakpoints = self.reference_sequence[prev_ref_node_end+1:breakpoint_position+1]
                    if sequence_between_breakpoints == "":
                        logging.error("Empty sequence on ref node. Prev ref node ended at %d, breakpoint position: %s" % (prev_ref_node_end, breakpoint_position))
                        logging.error("Is your reference genome matching the vcf?")
                        assert False

                    prev_ref_node_id = self._make_node(prev_ref_node_end, breakpoint_position+1,
                                                       sequence_between_breakpoints, is_ref_node=True)
                except EmptyNodeException:
                    logging.error("Tried making a node between pos %d and %d" % (prev_ref_node_end, breakpoint_position+1))
                    logging.info("Variant is %s" % variant)
                    raise

                prev_ref_node_end = breakpoint_position
                #logging.info("    Setting prev ref node end to %d" % prev_ref_node_end)

            # If breakpoint is a new variant, make a variant node
            if variant is not None:
                if variant.type != "DELETION":
                    #logging.info("   Making a variant node %d. Pos before/after: %d/%d. Variant node sequence: %s" % (
                    #self._current_node_id, variant.get_reference_position_before_variant(),
                    #variant.get_reference_position_after_variant(), variant.get_variant_sequence()))
                    self._make_node(variant.get_reference_position_before_variant(), variant.get_reference_position_after_variant(), variant.get_variant_sequence())
                else:
                    #logging.info("Adding deletion info from ref pos %d to ref pos %d" % (variant.get_reference_position_before_variant()+1, variant.get_reference_position_after_variant()))
                    self._deletions[variant.get_reference_position_before_variant()+1].append(variant.get_reference_position_after_variant())
                    #self._node_to_ref_pos_after[prev_ref_node_id].append(variant.get_reference_position_after_variant())

        # Add one reference node for the rest of the reference sequence that is left
        self._make_node(prev_ref_node_end, len(self.reference_sequence), self.reference_sequence[prev_ref_node_end+1:len(self.reference_sequence)], is_ref_node=True)

    def make_edges(self):
        #print("Ref pos to node after: %s" % self._ref_pos_to_node_after)
        #print("Ref pos to node start: %s" % self._ref_pos_to_node_start)
        #print("Node to ref pos after: %s" % self._node_to_ref_pos_after)
        #print("Node to ref pos before: %s" % self._node_to_ref_pos_before)
        #logging.info("Adding edges")
        i = 0
        # we want an edge from the node to nodes that come after this node
        # we find nodes coming after by looking at which ref position this node goes to
        # then we find all nodes where the ref position before that ref pos goes to the node
        for from_node, ref_positions_after in self._node_to_ref_pos_after.items():
            if i % 500000 == 0:
                logging.info("%d edge nodes processed" % i)
            i += 1
            for ref_pos_after in ref_positions_after:
                # Find all nodes where the ref pos before ref_pos_after goes to the node
                # ignore nodes with same ref pos before annd after as this node
                # such nodes are parallel nodes to this node and we don't want to connect them
                for to_node in self._ref_pos_to_node_after[ref_pos_after-1]:
                    if to_node == from_node:
                        continue
                    if self._node_to_ref_pos_after[from_node] == self._node_to_ref_pos_after[to_node] \
                            and self._node_to_ref_pos_before[from_node] == self._node_to_ref_pos_before[to_node]:
                        #print("Skipping edge from %d to %d" % (from_node, to_node))
                        continue
                    #print("Adding edge from %d to %d" % (from_node, to_node))
                    self._make_edge(from_node, to_node, ref_pos_after-1)


    def make_sorted_breakpoints(self):
        # Assume variants are already sorted
        logging.info("Traversing variants")
        for i, variant in enumerate(self.variants):
            if i % 10000 == 0:
                logging.info("%d variants processed" % i)

            # Breakpoints are last base pair in a reference node
            self.breakpoints.append((variant.get_reference_position_before_variant(), variant))
            self.breakpoints.append((variant.get_reference_position_after_variant()-1, None))

        logging.info("Sorting breakpoints")
        self.breakpoints = sorted(self.breakpoints, key=lambda b: b[0])
        logging.info("Done sorting breakpoints")

        #logging.info("Breakpoints: %s" % self.breakpoints)
