import itertools
import logging
from multiprocessing.pool import Pool

import numpy as np

from kage.models.helper_model import make_helper_model_from_genotype_matrix
from obgraph.genotype_matrix import GenotypeMatrix
from obgraph.variant_to_nodes import VariantToNodes
from shared_memory_wrapper import from_shared_memory, to_shared_memory, to_file
from shared_memory_wrapper.util import interval_chunks
from .helper_model import HelperVariants, CombinationMatrix


def create_helper_model_single_thread(data):
    interval, args = data
    from_variant, to_variant = interval

    genotype_matrix = from_shared_memory(
        GenotypeMatrix, "genotype_matrix" + args.shared_memory_unique_id
    )

    # read genotype matrix etc. from shared memory
    submatrix = GenotypeMatrix(genotype_matrix.matrix[from_variant:to_variant, :])
    logging.info(
        "Creating helper model for %d individuals and %d variants"
        % (submatrix.matrix.shape[1], submatrix.matrix.shape[0])
    )

    subhelpers, subcombo = make_helper_model_from_genotype_matrix(
        submatrix.matrix, None, window_size=args.window_size
    )

    # variant ids in results are now from 0 to (to_variant-from_variant)
    subhelpers += from_variant
    return from_variant, to_variant, subhelpers, subcombo


def create_helper_model(args):
    args.shared_memory_unique_id = str(np.random.randint(0, 1e15))
    pool = Pool(args.n_threads)
    logging.info("Made pool")

    variant_to_nodes = VariantToNodes.from_file(args.variant_to_nodes)
    genotype_matrix = GenotypeMatrix.from_file(args.genotype_matrix)
    # NB: Transpose
    genotype_matrix.matrix = genotype_matrix.matrix.transpose().astype(np.int64)

    n_variants = len(variant_to_nodes.ref_nodes)
    n_threads = args.n_threads
    while n_variants < n_threads * 50 and n_threads > 2:
        n_threads -= 1
        logging.info("Lowered n threads to %d so that not too few variants are analysed together" % n_threads)

    variant_intervals = interval_chunks(0, n_variants, n_threads)
    logging.info("Will process variant intervals: %s" % variant_intervals)

    helpers = np.zeros(n_variants, dtype=np.int64)
    genotype_matrix_combo = np.zeros((n_variants, 3, 3), dtype=np.float64)

    logging.info("Putting data in shared memory")
    # put data in shared memory
    to_shared_memory(
        genotype_matrix, "genotype_matrix" + args.shared_memory_unique_id
    )
    to_shared_memory(
        variant_to_nodes, "variant_to_nodes" + args.shared_memory_unique_id
    )

    logging.info("Put data in shared memory")

    for from_variant, to_variant, subhelpers, subcombo in pool.imap(
            create_helper_model_single_thread, zip(variant_intervals, itertools.repeat(args))
    ):
        logging.info(f"Done with one chunk: {from_variant}, {to_variant}")
        helpers[from_variant:to_variant] = subhelpers
        genotype_matrix_combo[from_variant:to_variant] = subcombo

    genotype_matrix_combo = genotype_matrix_combo.astype(np.float64)

    to_file(HelperVariants(helpers), args.out_file_name + ".pkl")
    logging.info("Saved helper model to file: %s" % args.out_file_name + ".pkl")
    to_file(CombinationMatrix(genotype_matrix_combo), args.out_file_name + "_combo_matrix.pkl")
    logging.info(
        "Saved combo matrix to file %s" % args.out_file_name + "_combo_matrix.pkl"
    )
