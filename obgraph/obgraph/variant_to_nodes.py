import logging

import dill
import numpy as np

from .graph import VariantNotFoundException


class VariantToNodes:
    properties = {"ref_nodes", "var_nodes"}
    def __init__(self, ref_nodes=None, var_nodes=None):
        self.ref_nodes = ref_nodes
        self.var_nodes = var_nodes

    @classmethod
    def from_file(cls, file_name):
        try:
            with open(file_name, 'rb') as f:
                data = dill.load(f)
        except FileNotFoundError:
            with open(file_name + ".pkl", 'rb') as f:
                data = dill.load(f)

        logging.info("Loaded from %s" % file_name)
        return data

    def to_file(self, file_name):
        if not file_name.endswith(".pkl"):
            file_name = file_name + ".pkl"

        with open(file_name, 'wb') as f:
            dill.dump(self, f)
        logging.info("Saved to %s" % file_name)

    def slice(self, from_variant, to_variant):
        return VariantToNodes(self.ref_nodes[from_variant:to_variant], self.var_nodes[from_variant:to_variant])

    @classmethod
    def from_graph_and_variants(cls, graph, variants):
        n_variants = len(variants)
        var_nodes = np.zeros(n_variants, dtype=np.uint32)
        ref_nodes = np.zeros(n_variants, dtype=np.uint32)

        max_graph_node = graph.max_node_id()
        for i, variant in enumerate(variants):
            if i % 100000 == 0:
                logging.info("%d variants processed" % i)
            try:
                ref_node, var_node = graph.get_variant_nodes(variant)
            except VariantNotFoundException as e:
                logging.error(str(e))
                logging.error("Could not find variant, aborting")
                raise

            var_nodes[i] = var_node
            ref_nodes[i] = ref_node

            assert var_node <= max_graph_node
            assert ref_node <= max_graph_node, "Found ref node %d which is not <= max graph node %d. Variant %s" % (ref_node, max_graph_node, variant)

        return cls(ref_nodes, var_nodes)


    def __iter__(self):
        return zip(self.ref_nodes, self.var_nodes)

    def len(self):
        return len(self.ref_nodes)


class NodeToVariants:
    properties = {"index"}
    def __init__(self, index):
        self.index = index

    @classmethod
    def from_file(cls, file_name):
        try:
            with open(file_name, 'rb') as f:
                data = dill.load(f)
        except FileNotFoundError:
            with open(file_name + ".npz", 'rb') as f:
                data = dill.load(f)

        logging.info("Loaded from %s" % file_name)
        return data

    def to_file(self, file_name):
        if not file_name.endswith(".pkl"):
            file_name = file_name + ".pkl"

        with open(file_name, 'wb') as f:
            dill.dump(self, f)
        logging.info("Saved to %s" % file_name)

    def get_variant_at_node(self, node):
        if self.index[node] == -1:
            return None
        return self.index[node]

    @classmethod
    def from_variant_to_nodes(cls, variant_to_nodes):
        n_nodes = max(np.max(variant_to_nodes.ref_nodes), np.max(variant_to_nodes.var_nodes))
        index = np.zeros(n_nodes+1, dtype=np.int32) - 1

        for variant in range(len(variant_to_nodes.ref_nodes+1)):
            if variant % 100000 == 0:
                logging.info("%d variants processed" % variant)

            ref_node = variant_to_nodes.ref_nodes[variant]
            var_node = variant_to_nodes.var_nodes[variant]

            index[ref_node] = variant
            index[var_node] = variant

        return cls(index)
