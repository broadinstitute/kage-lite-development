import logging

import numpy as np

from kage.models.joint_distribution import create_combined_matrices
from shared_memory_wrapper import from_file

MAIN = -1
HELPER = -2
M = MAIN
H = HELPER


class HelperVariants:
    properties = {"helper_variants"}

    def __init__(self, helper_variants):
        self.helper_variants = helper_variants

    @classmethod
    def from_file(cls, file_name):
        return from_file(file_name)


class CombinationMatrix:
    properties = {"matrix"}

    def __init__(self, matrix):
        self.matrix = matrix

    def __getitem__(self, variant_id):
        return self.matrix[variant_id]

    @classmethod
    def from_file(cls, file_name):
        return from_file(file_name)


def calc_likelihood(count_matrix):
    dummy_count = 1.
    count_matrix = count_matrix + dummy_count
    count_matrix = count_matrix.astype(np.float32)

    return (
        np.log(count_matrix[:, 0, 0])
        + np.log(count_matrix[:, 1, 1])
        + np.log(count_matrix[:, 2, 2])
    )


def calc_argmax(count_matrix):
    return np.sum(np.max(count_matrix, axis=M), axis=-1) / count_matrix.sum(axis=(M, H))


def get_helper_posterior(genotype_combo_matrix, global_helper_weight=5):
    helper_sum = np.sum(genotype_combo_matrix, axis=M, keepdims=True).astype(np.uint32)
    assert helper_sum[0].shape == (3, 1)
    global_helper_prior = (
        np.mean(helper_sum, axis=0, keepdims=True) + 1 / genotype_combo_matrix.shape[0]
    )
    logging.debug("Global helper prior: \n%s" % global_helper_prior)
    logging.debug("Helper sum: \n%s" % helper_sum)
    assert global_helper_prior.shape == (1, 3, 1)
    helper_posterior = (
        global_helper_prior / global_helper_prior.sum() * global_helper_weight
        + helper_sum
    )
    helper_posterior = helper_posterior / helper_posterior.sum(axis=H, keepdims=True)
    assert np.allclose(helper_posterior.sum(axis=H), 1), helper_posterior
    return helper_posterior


def get_population_priors(
    genotype_combo_matrix,
    weight=150,
    weight_diagonal=0,
    weight_left_column=0,
    weight_global=1,
):
    """n_variants x helper x main"""
    prior = np.eye(3, dtype=np.float32) * weight_diagonal
    prior[:, 0] = weight_left_column  # going to 0/0 is high
    prior += weight_global
    logging.debug("Weights added to population priors: \n%s" % prior)
    mean = np.sum(genotype_combo_matrix, axis=0) + prior
    logging.debug("Population prior before weighted: \n%s" % mean)
    weighted = mean.astype(np.float32) / np.sum(mean, axis=M, keepdims=True) * weight  # helper_sum*weight
    logging.debug("Population prior after weighted: \n%s" % weighted)
    return weighted


def make_helper_model_from_genotype_matrix(
    genotype_matrix,
    most_similar_variant_lookup=False,
    score_func=calc_likelihood,
    window_size=1000,
):
    # genotype_matrix = convert_genotype_matrix(genotype_matrix)
    logging.info("Finding best helper")

    if most_similar_variant_lookup is not None:
        logging.info("Making from most similar variant lookup")
        helpers = most_similar_variant_lookup.lookup_array
    else:
        logging.info(
            "Making raw from genotype matrix with window size %d" % window_size
        )
        combined = create_combined_matrices(genotype_matrix, window_size)
        helpers = find_best_helper(
            combined,
            score_func,
            len(genotype_matrix),
            with_model=score_func != calc_likelihood,
        )

    helper_counts = genotype_matrix[helpers] * 3
    flat_idx = (genotype_matrix + helper_counts).astype(np.uint32)

    genotype_combo_matrix = np.array(
        [np.sum(flat_idx == k, axis=1) for k in np.arange(9, dtype=np.uint8)]
    ).T.reshape(-1, 3, 3).astype(np.uint8)

    logging.debug("########")
    logging.debug("Genotype combo matrix raw:")
    logging.debug(genotype_combo_matrix)
    logging.debug("Genotype combo matrix dtype:")
    logging.debug(genotype_combo_matrix.dtype)

    population_prior = get_population_priors(genotype_combo_matrix)
    helper_posterior = get_helper_posterior(genotype_combo_matrix)

    logging.debug("Population prior: \n%s" % population_prior)
    logging.debug("Helper posterior:\n%s", helper_posterior)
    logging.debug("Combo matrix before population priors:\n%s" % genotype_combo_matrix)
    population_posterior = genotype_combo_matrix + population_prior
    logging.debug("Combo matrix after population priors:\n%s" % population_posterior)
    population_posterior = (
        population_posterior
        / population_posterior.sum(axis=M, keepdims=True)
        * helper_posterior
    )
    logging.debug("Genotype combo matrix posterior: ")
    logging.debug(population_posterior)

    assert len(helpers) == genotype_matrix.shape[0]

    return helpers, population_posterior


def find_best_helper(combined, score_func, N, with_model=False):
    best_idx, best_score = np.zeros(N, dtype=np.uint32), -np.inf * np.ones(N, dtype=np.float32)
    for j, counts in enumerate(combined, 1):
        if j < 4:
            continue
        if j % 50 == 0:
            logging.info("Window %d" % j)
        scores = score_func(counts, j) if with_model else score_func(counts)

        # TODO result of this floating-point comparison depends on whether CPU has AVX512 instructions and
        # leads to reproducibility issues on GitHub Actions runners; we disable AVX512 in the Docker environment
        # by setting the environment variable NPY_DISABLE_CPU_FEATURES appropriately
        do_update = scores > best_score[j:]
        best_score[j:][do_update] = scores[do_update]
        best_idx[j:][do_update] = np.flatnonzero(do_update)
        # reverse
        rev_scores = (
            score_func(counts.swapaxes(-2, -1), -j)
            if with_model
            else score_func(counts.swapaxes(-2, -1))
        )
        do_update = rev_scores >= best_score[:-j]
        best_score[:-j][do_update] = rev_scores[do_update]
        best_idx[:-j][do_update] = np.flatnonzero(do_update) + j
    return best_idx
