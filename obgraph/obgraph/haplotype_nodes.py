import logging

import numpy as np

from shared_memory_wrapper import from_file, to_file
from .util import phased_genotype_matrix_to_haplotype_matrix


class DiscBackedRaggedArray:
    def __init__(self, file_name, offsets, lengths):
        self._file_name = file_name
        self._offsets = offsets
        self._lengths = lengths

    def __getitem__(self, row: int):
        assert isinstance(row, int)
        out = np.fromfile(self._file_name,
                          offset=self._offsets[row] * 8,  # *8 because bytes
                          count=self._lengths[row],
                          dtype=np.int64)
        return out

    @classmethod
    def from_iter(cls, file_name, data):
        with open(file_name, "wb") as out_file:

            offsets = []
            lengths = []
            offset = 0

            for row_data in data:
                assert row_data.dtype == np.int64
                offsets.append(offset)
                length = len(row_data)
                lengths.append(length)
                offset += length
                row_data.tofile(out_file)

            return cls(file_name, np.array(offsets, dtype=np.uint32), np.array(lengths, dtype=np.uint32))

class DiscBackedHaplotypeToNodes:
    def __init__(self, data: DiscBackedRaggedArray):
        self.data = data

    def n_haplotypes(self):
        return len(self.data._offsets)

    def get_nodes(self, haplotype):
        return self.data[haplotype]

    def __getitem__(self, item):
        return self.data[item]


    @classmethod
    def from_file(cls, file_name):
        o = from_file(file_name)
        # hack with file name to get relative path correct
        o.data._file_name = file_name + ".haplotype_nodes"
        return o

    def to_file(self, file_name):
        to_file(self, file_name)

    @classmethod
    def from_phased_genotype_matrix(cls, genotype_matrix, variant_to_nodes, out_file_name):

        n_variants = genotype_matrix.shape[0]
        logging.info("%d variants" % n_variants)

        def get_nodes():
            for individual in range(genotype_matrix.shape[1]):
                has_variant = phased_genotype_matrix_to_haplotype_matrix(
                    genotype_matrix[:, individual].reshape(n_variants, 1))
                assert has_variant.shape[1] == 2
                for haplotype in [0, 1]:
                    has_node = np.nonzero(has_variant[:, haplotype])[0]
                    nodes = variant_to_nodes.var_nodes[has_node].astype(np.int64)
                    assert np.max(nodes) <= np.max(variant_to_nodes.var_nodes)
                    assert nodes.dtype == np.int64
                    yield nodes

        ra = DiscBackedRaggedArray.from_iter(out_file_name + ".haplotype_nodes", (nodes for nodes in get_nodes()))
        return cls(ra)
