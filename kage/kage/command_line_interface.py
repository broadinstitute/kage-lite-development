import logging
import sys

from .util import vcf_pl_and_gl_header_lines, convert_string_genotypes_to_numeric_array, _write_genotype_debug_data

logging.basicConfig(
    stream=sys.stderr,
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
)
import argparse, time, random
from kage.models.helper_index import create_helper_model
from .configuration import GenotypingConfig
from kage.models.mapping_model import sample_node_counts_from_population_cli, refine_sampling_model
from shared_memory_wrapper import (
    remove_shared_memory_in_session,
    get_shared_pool,
    close_shared_pool, )
from shared_memory_wrapper.shared_memory import remove_all_shared_memory
import numpy as np
from .indexing.tricky_variants import find_tricky_variants
from .indexing.index_bundle import IndexBundle
from .genotyping.combination_model_genotyper import CombinationModelGenotyper

np.random.seed(1)
np.seterr(all="ignore")
np.set_printoptions(suppress=True)


def main():
    run_argument_parser(sys.argv[1:])


def genotype(args):
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.debug("Will show debug")

    start_time = time.perf_counter()
    logging.info("Read coverage is set to %.3f" % args.average_coverage)
    get_shared_pool(args.n_threads)

    logging.info("Reading all indexes from an index bundle")
    t = time.perf_counter()
    index = IndexBundle.from_file(args.index_bundle).indexes
    logging.debug("Reading indexes took %.3f sec" % (time.perf_counter()-t))
    config = GenotypingConfig.from_command_line_args(args)

    node_counts = np.load(args.counts)

    max_variant_id = len(index.variant_to_nodes.ref_nodes) - 1
    logging.info("Max variant id is assumed to be %d" % max_variant_id)
    if args.limit_model_counts > 0:
        logging.info("Making model smaller (ignoring counts > %d)" % args.limit_model_counts)
        for i, count_model in enumerate(index.count_model):
            index.count_model[i].limit_to_n_individuals(args.limit_model_counts)


    genotyper = CombinationModelGenotyper(0, max_variant_id, node_counts, index, config=config)
    genotypes, probs, count_probs = genotyper.genotype()

    t = time.perf_counter()
    numpy_genotypes = convert_string_genotypes_to_numeric_array(genotypes)
    logging.debug("Converting string genotypes to numeric took %.4f sec" % (time.perf_counter()-t))

    if args.debug:
        _write_genotype_debug_data(genotypes, numpy_genotypes, args.out_file_name, index.variant_to_nodes, probs, count_probs)

    numpy_variants = index.numpy_variants

    t = time.perf_counter()
    numpy_variants.to_vcf_with_genotypes(
        args.out_file_name,
        config.sample_name_output,
        numpy_genotypes,
        add_header_lines=vcf_pl_and_gl_header_lines(),
        ignore_homo_ref=config.ignore_homo_ref,
        add_genotype_likelihoods=probs,
    )
    logging.info("Writing to vcf took %.3f sec" % (time.perf_counter() - t))

    close_shared_pool()
    logging.info("Genotyping took %d sec" % (time.perf_counter() - start_time))


def run_argument_parser(args):
    parser = argparse.ArgumentParser(
        description="kage",
        prog="kage",
        formatter_class=lambda prog: argparse.HelpFormatter(
            prog, max_help_position=50, width=100
        ),
    )

    subparsers = parser.add_subparsers()


    subparser = subparsers.add_parser("genotype")
    subparser.add_argument("-c", "--counts", required=False)
    subparser.add_argument("-r", "--reads", required=False)
    subparser.add_argument("-k", "--kmer_size", required=False, type=int, default=31)
    subparser.add_argument("-g", "--gpu", required=False, type=bool, default=False)
    subparser.add_argument("-i", "--index-bundle", required=True)
    subparser.add_argument("-m", "--kmer-index", required=False, help="Can be specified to override kmer index in index bundle for mapping.")
    subparser.add_argument("-o", "--out-file-name", required=True, help="Will write genotyped variants to this file")
    subparser.add_argument("-t", "--n-threads", type=int, required=False, default=8)
    subparser.add_argument("-a", "--average-coverage", type=float, default=15, help="Expected average read coverage", )
    subparser.add_argument("-q", "--min-genotype-quality", type=float, default=0.0,
        help="Min prob of genotype being correct. Genotypes with prob less than this are set to homo ref.")
    subparser.add_argument( "-s", "--sample-name-output", required=False, default="DONOR", help="Sample name that will be used in the output vcf")
    subparser.add_argument( "-u", "--use-naive-priors", required=False, type=bool, default=False,
        help="Set to True to use only population allele frequencies as priors.")
    subparser.add_argument("-l", "--limit-model-counts", default=0, type=int, help="If larger than 0, model will ignore counts larger than this. Can be used to use lower memory, but will make model less accurate.")
    subparser.add_argument("-I", "--ignore-helper-model", required=False, type=bool, default=False)
    subparser.add_argument("-b", "--ignore-homo-ref", required=False, type=bool, default=False, help="Set to True to not write homo ref variants to output vcf")
    subparser.add_argument("-B", "--do-not-write-genotype-likelihoods", required=False, type=bool, default=False, help="Set to True to not write genotype likelihoods to output vcf")
    subparser.add_argument("-d", "--debug", type=bool, default=False)
    subparser.set_defaults(func=genotype)


    subparser = subparsers.add_parser("find_tricky_variants")
    subparser.add_argument("-v", "--variant-to-nodes", required=True)
    subparser.add_argument("-m", "--node-count-model", required=True)
    subparser.add_argument("-r", "--reverse-kmer-index", required=True)
    subparser.add_argument("-M", "--max-counts-model", required=False, type=int, default=3, help="If model count exceeds this number, variant is tricky")
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.add_argument("-u", "--only-allow-unique", required=False, type=bool, help="Only allow variants where all kmers are unique")
    subparser.set_defaults(func=find_tricky_variants)


    def remove_shared_memory_command_line(args):
        remove_all_shared_memory()

    subparser = subparsers.add_parser("free_memory")
    subparser.set_defaults(func=remove_shared_memory_command_line)


    subparser = subparsers.add_parser("create_helper_model")
    subparser.add_argument("-g", "--genotype-matrix", required=False)
    subparser.add_argument("-v", "--variant-to-nodes", required=True)
    subparser.add_argument("-m", "--most-similar-variants", required=False)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.add_argument("-t", "--n-threads", required=False, default=1, type=int)
    subparser.add_argument("-u", "--use-duplicate-counts", required=False, type=bool, default=False)
    subparser.add_argument("-w", "--window-size", required=False, default=50, type=int, help="Number of variants before/after considered as potential helper variant")
    subparser.set_defaults(func=create_helper_model)

    def make_index_bundle(args):
        bundle = IndexBundle.from_args(args)
        bundle.to_file(args.out_file_name, compress=True)
        logging.info("Wrote index bundle to file %s" % args.out_file_name)

    subparser = subparsers.add_parser("make_index_bundle")
    subparser.add_argument("-g", "--variant-to-nodes", required=True)
    subparser.add_argument("-v", "--numpy-variants", required=True)
    subparser.add_argument("-A", "--count-model", required=True, help="Node count model")
    subparser.add_argument("-o", "--out-file-name", required=True, help="Will write genotyped variants to this file")
    subparser.add_argument("-x", "--tricky-variants", required=True)
    subparser.add_argument("-f", "--helper-model", required=True)
    subparser.add_argument("-F", "--helper-model-combo-matrix", required=True)
    subparser.add_argument("-i", "--kmer-index", required=True)
    subparser.set_defaults(func=make_index_bundle)

    subparser = subparsers.add_parser("sample_node_counts_from_population")
    subparser.add_argument("-g", "--graph", required=True)
    subparser.add_argument("-i", "--kmer-index", required=True)
    subparser.add_argument("-k", "--kmer-size", required=False, type=int, default=31)
    subparser.add_argument("-H", "--haplotype-to-nodes", required=True)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.add_argument("-t", "--n-threads", required=False, type=int, default=1)
    subparser.add_argument("-M", "--max-count", required=False, type=int, default=30)
    subparser.add_argument("-l", "--limit-to-n-individuals", required=False, type=int, default=0)
    subparser.set_defaults(func=sample_node_counts_from_population_cli)

    subparser = subparsers.add_parser("refine_sampling_model")
    subparser.add_argument("-s", "--sampling_model", required=True)
    subparser.add_argument("-v", "--variant-to-nodes", required=True)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.set_defaults(func=refine_sampling_model)

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)
    args.func(args)
    remove_shared_memory_in_session()
