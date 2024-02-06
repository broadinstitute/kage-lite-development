import logging
import time
from dataclasses import dataclass
from typing import List

import numpy as np
import scipy
from scipy.special import logsumexp

from kage.models.models import Model
from shared_memory_wrapper.shared_memory import run_numpy_based_function_in_parallel


def fast_poisson_logpmf(k, r):
    # should return the same as scipy.stats.logpmf(k, r), but is  ~4x faster
    #return k  * np.log(r) - r - np.log(scipy.special.factorial(k))
    return k * np.log(r) - r - scipy.special.gammaln(k+1)


@dataclass
class LimitedFrequencySamplingComboModel(Model):
    # stores only number of individuals having counts up to a given limit
    diplotype_counts: List[np.ndarray]  # list of matrices, each matrix is n_variants x max count supported

    def limit_to_n_individuals(self, n):
        for i, count in enumerate(self.diplotype_counts):
            self.diplotype_counts[i] = count[:, 0:n].copy()

    def __post_init__(self):
        #self.diplotype_counts = [c.astype(float) for c in self.diplotype_counts]
        return

    def astype(self, dtype):
        for i in range(3):
            #logging.info(self.diplotype_counts[i])
            self.diplotype_counts[i] = self.diplotype_counts[i].astype(dtype)

    def add_error_rate(self, error_rate=0.1):
        pass

    @classmethod
    def create_empty(cls, n_variants, max_count=3):
        logging.info("Creating empty limited freq model with dimensions %d x %d "% (n_variants, max_count))
        ret = cls([np.zeros((n_variants, max_count), dtype=np.uint16) for i in range(3)])
        return ret

    def __add__(self, other):
        for i in range(3):
            self.diplotype_counts[i] += other.diplotype_counts[i]
        return self

    def __getitem__(self, item):
        return self.diplotype_counts[item]

    @classmethod
    def create_naive(cls, n_variants, max_count=3, prior=0):
        empty = np.zeros((n_variants, max_count))
        counts0 = empty.copy()
        counts0[:,0] = 1
        counts0[:,0] += prior
        counts1 = empty.copy()
        counts1[:,1] = 1
        counts1[:,1] += prior
        counts2 = empty.copy()
        counts2[:,2] = 1
        counts2[:,2] += prior
        return cls([counts0, counts1, counts2])

    def subset_on_nodes(self, nodes):
        return self.__class__([matrix[nodes,:] for matrix in self.diplotype_counts])

    def __eq__(self, other):
        return all(np.all(m1 == m2) for m1, m2 in zip(self.diplotype_counts, other.diplotype_counts))

    @staticmethod
    def _logpmf(observed_counts, counts, base_lambda, error_rate):
        #counts = counts.astype(np.float16)
        sums = np.sum(counts, axis=-1)[:, None]
        frequencies = np.log(counts / sums)
        poisson_lambda = (np.arange(counts.shape[1])[None,:] + error_rate) * base_lambda
        poisson_lambda = poisson_lambda.astype(np.float16)
        prob = fast_poisson_logpmf(observed_counts[:,None].astype(np.float16), poisson_lambda)
        prob = logsumexp(frequencies + prob, axis=-1)
        return prob

    @staticmethod
    def _logpmf2(observed_counts, counts, base_lambda, error_rate):
        # Identical to _logpmf2, some lower memory usage by combining stuff
        poisson_lambda = (np.arange(counts.shape[1])[None, :] + error_rate) * base_lambda
        poisson_lambda = poisson_lambda.astype(np.float16)

        prob = logsumexp(np.log(counts / np.sum(counts, axis=-1)[:, None]) +
                         fast_poisson_logpmf(observed_counts[:, None].astype(np.float16), poisson_lambda)
                         , axis=-1)
        return prob

    @staticmethod
    def _gpu_logpmf(observed_counts, counts, base_lambda, error_rate):
        try:
            import custats
            import cupy as cp
        except ImportError:
            logging.error("Cucounter and cupy must be installed to run GPU genotyping")
            raise

        logging.info("Putting in cuda memory")

        t0 = time.perf_counter()
        observed_counts = observed_counts.astype(np.int32)
        logging.info("As type took %.3f sec" % (time.perf_counter()-t0))
        t0 = time.perf_counter()
        observed_counts = cp.asanyarray(observed_counts)
        logging.info("Moving observed counts to cuda memory took %.3f sec" % (time.perf_counter()-t0))

        t0 = time.perf_counter()
        counts = counts.astype(np.float32)
        logging.info("changing counts dtype took %.3f sec" % (time.perf_counter()-t0))

        t0 = time.perf_counter()
        counts = cp.asanyarray(counts)
        logging.info("SHAPE/DTYPE counts: %s/%s" % (counts.shape, counts.dtype))
        logging.info("Moving counts to cuda memory took %.3f sec" % (time.perf_counter()-t0))

        res = custats.functions.experimental_logpmf(observed_counts, counts, base_lambda, error_rate)
        return cp.asnumpy(res)

    def logpmf(self, observed_counts, d, base_lambda=1.0, error_rate=0.01, gpu=False, n_threads=16):
        logging.info("Will use %d threads" % n_threads)
        logging.debug("base lambda in LimitedFreq model is %.3f" % base_lambda)
        logging.debug("Error rate is %.3f" % error_rate)
        logging.info("Will use GPU? %s" % gpu)
        t0 = time.perf_counter()
        counts = self.diplotype_counts[d]  # [:, 0:5]
        counts = counts.astype(np.float16)
        if gpu:
            logging.info("USING GPU")
            prob = LimitedFrequencySamplingComboModel._gpu_logpmf(observed_counts, counts, base_lambda, error_rate)
        else:
            prob = run_numpy_based_function_in_parallel(
                LimitedFrequencySamplingComboModel._logpmf2, n_threads, (observed_counts, counts, base_lambda, error_rate)
            )
        logging.info("Logpmf took %.4f sec" % (time.perf_counter()-t0))
        return prob

    def describe_node(self, node):
        description = "\n"
        for count in range(3):
            description += "Having %d copies: " % count
            description += ', '.join("%d: %.3f" % (i, self.diplotype_counts[count][node][i]) for i in np.nonzero(self.diplotype_counts[count][node])[0])
            #description += ', '.join("%d: %d" % (c, f) for c, f in np.unique(self.diplotype_counts[count][node], return_counts=True))
            description += "\n"

        return description

    def fill_empty_data(self, prior=0.1):
        t = time.perf_counter()
        logging.info("Prior is %.4f" % prior)
        expected_counts = []
        for diplotype in [0, 1, 2]:
            logging.info("Computing expected counts for genotype %d" % diplotype)
            m = self.diplotype_counts[diplotype]
            not_missing = np.where(np.sum(m, axis=-1) > 0)[0]
            expected = np.zeros(m.shape[0], dtype=float)
            expected[not_missing] = np.sum(np.arange(m.shape[1]) * m[not_missing], axis=-1) / np.sum(m[not_missing], axis=-1)
            expected_counts.append(expected)

        # add priors
        # assume naively counts for 1 is 1 more than expected as 0. Counts at 2 is 2 more than expected at 1
        e0, e1, e2 = expected_counts
        m0, m1, m2 = self.diplotype_counts
        max_count = m0.shape[1]-1
        logging.debug("Adding priors")
        logging.debug("Size of e: %s" % (e1.shape))
        positions = np.round(e1).astype(int) - 1
        logging.debug("Found positions")
        positions = np.maximum(0, positions)
        logging.debug("Found max")
        rows = np.arange(0, m0.shape[0])
        m0[rows, positions] += prior
        m1[rows,np.minimum(max_count, np.round(e0).astype(int) + 1)] += prior
        m2[rows,np.minimum(max_count, np.round(e0).astype(int) + 2)] += prior

        assert all([np.all(np.sum(c, axis=-1) > 0) for c in self.diplotype_counts])
        logging.debug("Filling empty data took %.2f sec " % (time.perf_counter()-t))

    def has_no_data(self, idx):
        missing = [np.sum(c[idx]) == 0 for c in self.diplotype_counts]
        # if 2 out of 3 is missing data, return True
        return sum(missing) >= 2

    def has_duplicates(self, idx):
        counts = [c[idx] for c in self.diplotype_counts]
        for i in range(3):
            # if there are counts outside position i, there are duplicates
            # also if any elements outside i, there are duplicates
            #if np.sum(counts[i]) > counts[i][i] * 5:
            if np.sum(counts[i][i+1:]) > 0:
            #if np.argmax(counts[i]) != i:
                return True

        return False