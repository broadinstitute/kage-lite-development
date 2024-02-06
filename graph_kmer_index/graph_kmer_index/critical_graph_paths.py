import logging

import numpy as np


class CriticalGraphPaths:
    def __init__(self, nodes, offsets, index=None):
        self.nodes = nodes
        self.offsets = offsets
        self._index = index

    def _make_index(self):
        if len(self.nodes) == 0:
            self._index = np.zeros(0)
            return

        self._index = np.zeros(np.max(self.nodes)+1, dtype=np.uint16)
        self._index[self.nodes] = self.offsets

    @classmethod
    def empty(cls):
        return cls(np.array([]), np.array([]), np.array([]))

    def is_critical(self, node, offset):
        if self._index is None:
            logging.info("Making critical paths index")
            logging.info("NOdes: %s, offsets: %s" % (self.nodes, self.offsets))
            self._make_index()

        if node >= len(self._index):
            return False

        return self._index[node] == offset

    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        return ((node, offset) for node, offset in zip(self.nodes, self.offsets))
