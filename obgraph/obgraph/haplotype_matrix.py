import logging
from multiprocessing import Pool

import numpy as np

from shared_memory_wrapper.shared_memory import to_shared_memory, from_shared_memory


class HaplotypeMatrix:
    properties = {"matrix"}
    def __init__(self, matrix=None):
        self.matrix = matrix

    def get_allele_frequency_for_nodes(self, nodes, node_to_variants, variant_to_nodes):
        # determine variant ids and which haplotype (ref or var) these nodes represent
        variant_ids = []
        haplotypes = []

        for node in nodes:
            variant_id = node_to_variants.get_variant_at_node(node)
            if variant_id is not None:
                variant_ids.append(variant_id)
                if variant_to_nodes.ref_nodes[variant_id] == node:
                    # ref node
                    haplotypes.append(1)
                else:
                    assert variant_to_nodes.var_nodes[variant_id] == node
                    haplotypes.append(2)

        if len(variant_ids) == 0:
            # no variants, meaning allele frequency is 1.0
            return 1.0

        return self.get_allele_frequency_for_variants_and_haplotypes(variant_ids, haplotypes)

    def get_allele_frequency_for_variants_and_haplotypes(self, variants, haplotypes):
        assert len(variants) == len(haplotypes)
        n_total_haplotypes = self.matrix.shape[0]
        # find the number of individual that have the given haplotypes at the given variants
        submatrix = self.matrix[:, variants]
        has_haplotypes = len(np.where(np.sum(submatrix == haplotypes, axis=1) == len(variants))[0])
        return has_haplotypes / n_total_haplotypes

    @classmethod
    def from_variants(cls, variants, n_individuals, n_variants, n_threads=10, chunk_size=10000):
        matrix = np.zeros((n_individuals*2, n_variants), dtype=np.uint8)
        matrix = cls(matrix)
        logging.info("Putting genotype matrix in shared memory")
        to_shared_memory(matrix, "haplotype_matrix")

        logging.info("Getting variant chunks")
        variant_chunks = variants.get_chunks(chunk_size=chunk_size)

        pool = Pool(n_threads)

        i = 0
        for result in pool.imap(HaplotypeMatrix.fill_shared_memory_matrix_with_variants, variant_chunks):
            i += 1
            logging.info("Done with %d variant chunks" % i)

        logging.info("Done with all variant chunks")
        matrix = from_shared_memory(HaplotypeMatrix, "haplotype_matrix")
        return cls(matrix.matrix)

    def get_alleles_with_haplotype_at_variant(self, variant_id, haplotype):
        assert haplotype [1, 2]
        return np.where(self.matrix[:, variant_id] == haplotype)[0]

    @staticmethod
    def fill_shared_memory_matrix_with_variants(variants):
        matrix = from_shared_memory(HaplotypeMatrix, "haplotype_matrix")

        for variant in variants:
            variant_number = variant.vcf_line_number
            if variant_number % 10000 == 0:
                logging.info("%d variants processeed" % variant_number)

            for individual_id, haplotypes in variant.get_individuals_and_numeric_haplotypes():
                matrix.matrix[individual_id*2, variant_number] = haplotypes[0]
                matrix.matrix[individual_id*2+1, variant_number] = haplotypes[1]

    def to_file(self, file_name):
        np.save(file_name, self.matrix)

    @classmethod
    def from_file(cls, file_name):
        try:
            data = np.load(file_name)
        except FileNotFoundError:
            data = np.load(file_name + ".npy")

        return cls(data)
