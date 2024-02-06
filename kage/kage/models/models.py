import logging
import time

import numpy as np
from scipy.special import logsumexp

from shared_memory_wrapper.shared_memory import run_numpy_based_function_in_parallel
from ..util import log_memory_usage_now

MAIN = -1
HELPER = -2
M = MAIN
H = HELPER


class Model:
    count_probs = None  # for outside access after computing

    def predict(self, k1, k2, return_probs=False):
        scores = self.score(k1, k2)
        result = np.argmax(scores, axis=-1)

        if return_probs:
            return result, scores

    def score(self, k1, k2):
        count_probs = np.array([self.logpmf(k1, k2, g) for g in [0, 1, 2]]).T
        self.count_probs = count_probs
        return count_probs

    def logpmf(self, k1, k2, genotype):
        return NotImplemented


class ComboModelBothAlleles(Model):
    def __init__(self, model_ref, model_alt):
        self._model_ref = model_ref
        self._model_alt = model_alt
        self._logpmf_cache = {}

    def compute_logpmfs(self, k1, k2):
        for genotype in [0, 1, 2]:
            self.logpmf(k1, k2, genotype)

    def clear(self):
        self._model_ref = None
        self._model_alt = None

    def logpmf(self, k1, k2, genotype, base_lambda=1.0, gpu=False, n_threads=16):
        if genotype in self._logpmf_cache:
            return self._logpmf_cache[genotype]

        logging.debug("Using base lambda %.3f in combo model both alleles" % base_lambda)
        logging.info("Model is %s" % type(self._model_ref))
        ref_probs = self._model_ref.logpmf(k1, 2 - genotype, base_lambda=base_lambda, gpu=gpu, n_threads=n_threads)
        alt_probs = self._model_alt.logpmf(k2, genotype, base_lambda=base_lambda, gpu=gpu, n_threads=n_threads)
        prob = ref_probs + alt_probs

        self._logpmf_cache[genotype] = prob
        return prob


class HelperModel(Model):
    def __init__(
        self, model, helper_variants, genotype_combo_matrix, tricky_variants=None, base_lambda=1.0,
            gpu=False, n_threads=16
    ):
        self._model = model
        self._helper_variants = helper_variants
        t = time.perf_counter()
        log_memory_usage_now("Before making _genotype_probs")
        self._genotype_probs = np.log(
            genotype_combo_matrix
            / genotype_combo_matrix.sum(axis=(-1, -2), keepdims=True)
        ).astype(np.float16)
        logging.info(
            "Computing genotype probs in HelperModel init took %.4f sec"
            % (time.perf_counter() - t)
        )
        self._tricky_variants = tricky_variants
        self.count_probs = None
        self._base_lambda = base_lambda
        self._gpu = gpu
        self._n_threads = n_threads

    def score(self, k1, k2):
        t0 = time.perf_counter()
        count_probs = np.array([self._model.logpmf(k1, k2, g, self._base_lambda, gpu=self._gpu, n_threads=self._n_threads) for g in [0, 1, 2]]).T
        logging.info("Getting count probs for all combinations of genotypes took %.5f sec" % (time.perf_counter()-t0))
        self.count_probs = count_probs

        if self._tricky_variants is not None:
            logging.info("Using tricky variants in HelperModel.score. There are %d tricky variants" % np.sum(self._tricky_variants))
            count_probs = np.where(self._tricky_variants.reshape(-1, 1), np.log(1 / 3), count_probs)

        time_start = time.perf_counter()
        t0 = time.perf_counter()
        log_probs = (
            self._genotype_probs
            + count_probs[self._helper_variants].reshape(-1, 3, 1)
            + count_probs.reshape(-1, 1, 3)
        )
        logging.info("Computing log_probs took %.5f sec" % (time.perf_counter()-t0))


        logging.info(
            "Time spent on log_probs in HelperModel.score: %.4f"
            % (time.perf_counter() - time_start)
        )
        time_start = time.perf_counter()
        # result = logsumexp(log_probs, axis=H)
        # result = result - logsumexp(result, axis=-1, keepdims=True)
        t0 = time.perf_counter()
        result = run_numpy_based_function_in_parallel(
            lambda p: logsumexp(p, axis=H), self._n_threads, [log_probs]
        )
        logging.info("Computinng logsumex across helper axis took %.5f sec" % (time.perf_counter() - t0))

        t0 = time.perf_counter()
        result = run_numpy_based_function_in_parallel(
            lambda result: result - logsumexp(result, axis=-1, keepdims=True),
            self._n_threads,
            [result],
        )
        logging.info("Computinng result-logsumexp(result) took %.5f sec" % (time.perf_counter() - t0))
        logging.info(
            "Time spent to compute probs using helper probs in HelperModel.score: %.4f"
            % (time.perf_counter() - time_start)
        )
        return result

    def logpmf(self, ref_counts, alt_counts, genotype):
        #return NotImplemented
        count_probs = np.array(
            [self._model.logpmf(ref_counts, alt_counts, g) for g in [0, 1, 2]]
        ).T
        log_probs = (
            self._genotype_probs
            + count_probs[self._helper_variants].reshape(-1, 3, 1)
            + count_probs.reshape(-1, 1, 3)
        )
        return logsumexp(log_probs, axis=H)[..., genotype]
