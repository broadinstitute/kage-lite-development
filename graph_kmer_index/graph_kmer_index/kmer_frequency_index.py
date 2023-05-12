import logging

import dill
import numpy as np


class KmerFrequencyIndex:
    def __init__(self, kmers, frequencies):
        self._kmers = kmers
        self._frequencies = frequencies

    def get(self, kmer):
        index = np.searchsorted(self._kmers, kmer, side="right")
        if self._kmers[index] == kmer:
            return self._frequencies[index]

        logging.warning("No hit for kmer %d" % kmer)
        return 0

    @classmethod
    def from_kmers(cls, kmers):
        logging.info("Sorting")
        sorting = np.argsort(kmers)
        kmers = kmers[sorting]
        logging.info("Counting")
        unique, frequencies = np.unique(kmers, return_counts=True)
        return cls(unique, frequencies)

    @classmethod
    def from_file(cls, file_name):
        try:
            with open(file_name + ".npz", 'rb') as f:
                data = dill.load(f)
        except FileNotFoundError:
            with open(file_name, 'rb') as f:
                data = dill.load(f)

        return data

    def to_file(self, file_name):
        with open(file_name, 'wb') as f:
            dill.dump(self, f)
