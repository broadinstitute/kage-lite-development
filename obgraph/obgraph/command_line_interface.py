import logging

logging.basicConfig(level=logging.INFO, format='%(module)s %(asctime)s %(levelname)s: %(message)s')
import pyximport; pyximport.install()
import sys
import argparse
from .graph import Graph
from .variants import VcfVariants
from pyfaidx import Fasta
from .graph_construction import GraphConstructor
from .graph_merger import merge_graphs
from shared_memory_wrapper import remove_shared_memory_in_session, to_file, from_file
from .variant_to_nodes import VariantToNodes


def main():
    run_argument_parser(sys.argv[1:])


def run_argument_parser(args):
    parser = argparse.ArgumentParser(
        description='Obgraph.',
        prog='obgraph',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=50, width=100))

    def make(args):
        if args.vcf is not None:
            logging.info("Will create from vcf file")
            reference = Fasta(args.reference_fasta_file)

            chromosome = args.chromosome

            ref_sequence = str(reference[args.chromosome])
            logging.info("Extracted sequence for chromosome %s. Length is: %d" % (chromosome, len(ref_sequence)))
            variants = VcfVariants.from_vcf(args.vcf, limit_to_chromosome=chromosome)
            logging.info("There are %d variants in chromosome" % len(variants))

            constructor = GraphConstructor(ref_sequence, variants)
            graph = constructor.get_graph_with_dummy_nodes()
            graph.to_file(args.out_file_name)
        else:
            logging.info("Will create from files %s" % args.vg_json_files)
            graph = Graph.from_vg_json_files(args.vg_json_files)
            graph.to_file(args.out_file_name)

    subparsers = parser.add_subparsers()
    subparser = subparsers.add_parser("make")
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-j", "--vg-json-files", nargs='+', required=False)
    subparser.add_argument("-v", "--vcf", required=False)
    subparser.add_argument("-r", "--reference_fasta_file", required=False)
    subparser.add_argument("-c", "--chromosome", required=False)
    subparser.set_defaults(func=make)


    def add_allele_frequencies(args):
        logging.info("Reading graph")
        graph = Graph.from_file(args.graph_file_name)
        variants = VcfVariants.from_vcf(args.vcf_file_name, limit_to_chromosome=args.chromosome, skip_index=True)
        graph.set_allele_frequencies_from_variants(variants,
                                                   use_chromosome=args.chromosome)  # Use chromosome 1 because we always assume this is a single-chromosome graph
        graph.to_file(args.graph_file_name)
        logging.info("Wrote modified graph to the same file %s" % args.graph_file_name)

    subparser = subparsers.add_parser("add_allele_frequencies")
    subparser.add_argument("-g", "--graph-file-name", required=True)
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.add_argument("-c", "--chromosome", required=False, help="If vcf contains multiple chromosomes, use this to limit to the chromosome that the graph is made from")
    subparser.set_defaults(func=add_allele_frequencies)


    def make_haplotype_to_nodes_bnp(args):
        from .haplotype_nodes import make_ragged_haplotype_to_nodes
        phased_genotype_matrix = from_file(args.phased_genotype_matrix).matrix
        variant_to_nodes = VariantToNodes.from_file(args.variant_to_nodes)

        if args.make_disc_backed:
            from .haplotype_nodes import DiscBackedHaplotypeToNodes
            logging.info("Making disc backed")
            result = DiscBackedHaplotypeToNodes.from_phased_genotype_matrix(phased_genotype_matrix, variant_to_nodes, args.out_file_name)
            result.to_file(args.out_file_name)
        else:
            result = make_ragged_haplotype_to_nodes(variant_to_nodes, phased_genotype_matrix, args.n_threads)
            to_file(result, args.out_file_name)

    subparser = subparsers.add_parser("make_haplotype_to_nodes_bnp")
    subparser.add_argument("-g", "--variant-to-nodes", required=True)
    subparser.add_argument("-v", "--phased-genotype-matrix", required=True)
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-t", "--n-threads", type=int, default=8, required=False)
    subparser.add_argument("-n", "--n-haplotypes", type=int, required=False)
    subparser.add_argument("-d", "--make-disc-backed", type=bool, required=False, default=False, help="Uses less memory")
    subparser.set_defaults(func=make_haplotype_to_nodes_bnp)


    def merge_graphs_command(args):
        graphs = [Graph.from_file(graph) for graph in args.graphs]
        logging.info("Done reading graphs")

        merged_graph = merge_graphs(graphs)
        merged_graph.to_file(args.out_file_name)

    subparser = subparsers.add_parser("merge_graphs")
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-g", "--graphs", nargs="+", required=True)
    subparser.set_defaults(func=merge_graphs_command)


    def make_variant_to_nodes(args):
        from .variant_to_nodes import VariantToNodes
        graph = Graph.from_file(args.graph)
        variants = VcfVariants.from_vcf(args.vcf, skip_index=True)
        variant_to_nodes = VariantToNodes.from_graph_and_variants(graph, variants)
        variant_to_nodes.to_file(args.out_file_name)
        logging.info("Wrote to file %s" % args.out_file_name)

    subparser = subparsers.add_parser("make_variant_to_nodes")
    subparser.add_argument("-g", "--graph", required=True)
    subparser.add_argument("-v", "--vcf", required=True)
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.set_defaults(func=make_variant_to_nodes)


    def make_numpy_variants(args):
        from .numpy_variants import NumpyVariants
        n = NumpyVariants.from_vcf(args.vcf)
        n.to_file(args.out_file_name)

    subparser = subparsers.add_parser("make_numpy_variants")
    subparser.add_argument("-v", "--vcf", required=True)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.set_defaults(func=make_numpy_variants)


    def make_position_id(args):
        from .position_id import PositionId
        graph = Graph.from_file(args.graph)
        position_id = PositionId.from_graph(graph)
        to_file(position_id, args.out_file_name)

    subparser = subparsers.add_parser("make_position_id")
    subparser.add_argument("-g", "--graph", required=True)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.set_defaults(func=make_position_id)


    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)
    args.func(args)
    remove_shared_memory_in_session()