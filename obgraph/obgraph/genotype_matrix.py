import logging

from shared_memory_wrapper import to_file, from_file


class GenotypeMatrix:
    properties = {"matrix"}

    def __init__(self, matrix=None):
        if matrix is not None:
            logging.info("Type of matrix: %s" % matrix.dtype)
            logging.info("Size of matrix: %3.f MB" % (int(matrix.nbytes)/1000000))
        self.matrix = matrix

    def to_file(self, file_name):
        to_file(self, file_name)

    @classmethod
    def from_file(cls, file_name):
        return from_file(file_name)
