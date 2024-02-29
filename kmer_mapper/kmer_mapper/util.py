import gzip
import logging
import resource
import sys
from pathlib import PurePath

import bionumpy as bnp
import numpy as np

from graph_kmer_index.collision_free_kmer_index import CounterKmerIndex, CollisionFreeKmerIndex
from graph_kmer_index.index_bundle import IndexBundle
from shared_memory_wrapper import from_file


def read_fasta(file_name):
    i = 0
    with open(file_name, "rb") as f:
        for line in f:
            if line[0] != 62:
                i += 1
                yield line



def remap_array(array, from_values, to_values):
    index = np.digitize(array.ravel(), from_values, right=True)
    return to_values[index].reshape(array.shape)



def log_memory_usage_now(logplace=""):
    memory = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1000000
    logging.info("Memory usage (%s): %.4f GB" % (logplace, memory))


def _get_kmer_index_from_args(args):
    # allowing multiple ways to specify kmer index, try to get the index
    if isinstance(args.kmer_index, CollisionFreeKmerIndex):
        kmer_index = args.kmer_index
        kmer_index.convert_to_int32()
        kmer_index.remove_ref_offsets()
        return kmer_index

    if args.kmer_index is None:
        if args.index_bundle is None:
            logging.error("Either a kmer index (-i) or an index bundle (-b) needs to be specified")
            sys.exit(1)
        else:
            kmer_index = IndexBundle.from_file(args.index_bundle).indexes["KmerIndex"]
            kmer_index.convert_to_int32()
            kmer_index.remove_ref_offsets()  # not needed, will save us some memory
    else:

        cls = CollisionFreeKmerIndex
        try:
            kmer_index = cls.from_file(args.kmer_index)
            kmer_index.convert_to_int32()
            kmer_index.remove_ref_offsets()  # not needed, will save us some memory
        except:
            kmer_index = from_file(args.kmer_index)
            assert isinstance(kmer_index, CounterKmerIndex)
            logging.info("Kmer index is counter index")

    return kmer_index


def get_kmer_hashes_from_chunk_sequence(chunk_sequence, kmer_size):
    encoded_array = bnp.as_encoded_array(chunk_sequence)
    mask = np.sum(encoded_array == "N", axis=-1) > 0
    encoded_array = encoded_array[~mask]
    encoded_array = bnp.change_encoding(encoded_array, bnp.DNAEncoding)
    # encoded_array = bnp.as_encoded_array(chunk_sequence, bnp.DNAEncoding)
    hashes = bnp.sequence.get_kmers(encoded_array, kmer_size).ravel().raw().astype(np.uint64)
    logging.debug("N hashes: %d" % len(hashes))
    return hashes


def open_file(filename):
    path = PurePath(filename)

    is_gz = path.suffixes[-1] == ".gz"
    suffix = path.suffixes[-2] if is_gz else path.suffixes[-1]

    try:
        buffer_type = bnp.io.files._get_buffer_type(suffix)
    except RuntimeError:
        logging.error("Unsupported file suffix %s" % suffix)
        raise

    if suffix == ".fa" or suffix == ".fasta":
        # override buffer type for some performance gain
        buffer_type = bnp.TwoLineFastaBuffer
        logging.info("Using buffer type TwoLineFastaBuffer")

    open_func = gzip.open if is_gz or suffix == ".bam" else open
    reader = bnp.io.parser.NumpyFileReader(open_func(filename, "rb"), buffer_type)
    reader.set_prepend_mode()
    return reader
