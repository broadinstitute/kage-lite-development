import numpy as np

class PositionId:
    def __init__(self, index):
        self._index = index

    def get(self, nodes, offsets):
        return self._index[nodes]+offsets

    @classmethod
    def from_graph(cls, graph):
        node_sizes = graph.nodes
        node_ids = np.where(node_sizes!=0)[0]
        index = np.zeros(len(node_sizes)+1, dtype=np.int64)
        index[node_ids[1:]] = np.cumsum(node_sizes)[node_ids[1:]-1]
        index[node_ids[0]] = 0
        return cls(index)
