import numpy as np


class CoordinateConverter:
    def __init__(self, chromosome_offsets):
        self.chromosome_offsets = chromosome_offsets

    @classmethod
    def from_graph(cls, graph):
        chromosome_start_offsets = graph.node_to_ref_offset[graph.chromosome_start_nodes]
        return cls(chromosome_start_offsets)

    def convert_chromosome_offset(self, chromosome, offset):
        if chromosome == "X":
            chromosome = 23
        elif chromosome == "Y":
            chromosome = 24
        id = int(chromosome) - 1
        if id >= len(self.chromosome_offsets):
            return None

        return int(self.chromosome_offsets[id] + offset)

    def convert_offset_to_chromosome_and_offset(self, offset):
        index = np.searchsorted(self.chromosome_offsets, offset) - 1
        chromosome = str(index + 1)
        offset = offset - self.chromosome_offsets[index]

        return chromosome, int(offset)

    @classmethod
    def from_file(cls, file_name):
        data = np.load(file_name)
        return cls(data)

    def to_file(self, file_name):
        np.save(file_name, self.chromosome_offsets)