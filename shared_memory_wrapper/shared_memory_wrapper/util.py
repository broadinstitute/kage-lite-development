import logging
import resource
from collections.abc import Iterable

import numpy as np

from .shared_memory import get_shared_pool, close_shared_pool
from .shared_memory import remove_shared_memory
from .shared_memory_v2 import object_to_shared_memory, object_from_shared_memory


def interval_chunks(start, end, n_chunks):
    assert end > start
    assert n_chunks > 0
    step = max(1, (end-start)//n_chunks)
    boundaries = list(range(start, end, step))
    boundaries.append(end)
    return [(start, end) for start, end in zip(boundaries[0:-1], boundaries[1:])]


class AddReducer():
    def __init__(self, initial):
        self.result = initial

    def __call__(self):
        return self.get_final_result()

    def add_result(self, result):
        self.result += result

    def get_final_result(self):
        return self.result


def subchunker(iterable, n=10):
    i = 0
    for i, element in enumerate(iterable):
        yield element
        if i >= n-1:
            break

    if i == 0:
        return


def chunker(iterable, n=10):
    while True:
        subchunk = subchunker(iterable, n=n)
        yield subchunk


class FunctionWrapper:
    def __init__(self, function, data_id, backend=None):
        self.function = function
        self.data_id = data_id
        self._backend = backend

    def __call__(self, run_specific_variable):
        data = object_from_shared_memory(self.data_id, backend=self._backend)
        result = self.function(*data, run_specific_variable)
        return result


def parallel_map_reduce_with_adding(function, data, mapper, initial_data, n_threads=8, chunk_size=50):
    reducer = AddReducer(initial_data)
    return parallel_map_reduce(function, data, mapper, reducer, n_threads, chunk_size=chunk_size)


def chunked_imap(pool, function, iterable, chunk_size=10):
    logging.info("running chunked imap with chunk_size %d" % chunk_size)
    # runs pool.imap but on chunks to lower memory
    chunks = chunker(iterable, chunk_size)
    for chunk in chunks:
        logging.info("Processing chunk in chunked_imap")
        # run only imap on this chunk
        for i, result in enumerate(pool.imap(function, chunk)):
            yield result

        if i < chunk_size-1:
            logging.info("No more results, returning")
            return


def parallel_map_reduce(function, data, mapper, reducer=None, n_threads=7, backend="shared_array", chunk_size=50):
    assert reducer is None or isinstance(reducer, AddReducer)
    assert isinstance(mapper, Iterable), "Mapper must be iterable"
    assert isinstance(data, tuple)

    logging.info("Putting data in shared memory")
    data = object_to_shared_memory(data, backend=backend)
    logging.info("Done putting data in shared memory")
    pool = get_shared_pool(n_threads)

    function = FunctionWrapper(function, data, backend=backend)

    #for result in pool.imap(function, mapper):
    for i, result in enumerate(chunked_imap(pool, function, mapper, chunk_size=chunk_size)):
        logging.info("Result %d returned" % i)
        if reducer is not None:
            assert result is not None
            reducer.add_result(result)

    close_shared_pool()

    if reducer is not None:
        out = reducer.get_final_result()
    else:
        out = object_from_shared_memory(data, backend=backend)

    remove_shared_memory(data)
    return out
