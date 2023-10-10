import logging
import sys
import argparse
import numpy as np
import allel
from shared_memory_wrapper import to_file, from_file


logging.basicConfig(
    stream=sys.stderr,
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
)


np.random.seed(1)
np.seterr(all="ignore")
np.set_printoptions(suppress=True)


def main():
    run_argument_parser(sys.argv[1:])


def run_argument_parser(args):
    parser = argparse.ArgumentParser(
        description="lite_utils",
        prog="lite_utils",
        formatter_class=lambda prog: argparse.HelpFormatter(
            prog, max_help_position=50, width=100
        ),
    )

    subparsers = parser.add_subparsers()

    def make_genotype_matrix(args):
        from obgraph.genotype_matrix import GenotypeMatrix

        logging.info(f'Reading genotype matrix from {args.vcf_file_name}...')
        gt_vsp = allel.read_vcf(args.vcf_file_name, region=args.chromosome, fields='calldata/GT', log=sys.stderr)['calldata/GT']

        # convert no-calls to ref
        gt_vsp[gt_vsp == -1] = 0

        # convert to [0, 1, 2, 3] TODO assumes biallelic variants
        encoded_gt_vs = 2 * gt_vsp[:, :, 0] + gt_vsp[:, :, 1]

        assert np.all(encoded_gt_vs >= 0)
        assert np.all(encoded_gt_vs <= 3)

        matrix = GenotypeMatrix(encoded_gt_vs)
        matrix.to_file(args.out_file_name)

    subparser = subparsers.add_parser("make_genotype_matrix")
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.add_argument("-c", "--chromosome", required=True)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.set_defaults(func=make_genotype_matrix)

    def merge_genotype_matrices_and_convert_to_unphased(args):
        from obgraph.genotype_matrix import GenotypeMatrix

        genotype_matrices = [from_file(file).matrix for file in args.genotype_matrices]
        matrix = GenotypeMatrix(np.concatenate(genotype_matrices).transpose())

        # phased (hom ref = 0, het = 1 or 2, hom alt = 3) to unphased (hom ref = 0, het = 1, hom alt = 2)
        matrix.matrix[matrix.matrix == 2] = 1
        matrix.matrix[matrix.matrix == 3] = 2

        matrix.to_file(args.out_file_name)

    subparser = subparsers.add_parser("merge_genotype_matrices")
    subparser.add_argument("-g", "--genotype-matrices", nargs='+', required=True)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.set_defaults(func=merge_genotype_matrices_and_convert_to_unphased)

    def merge_chromosome_variant_to_nodes(args):
        from obgraph.variant_to_nodes import VariantToNodes

        ref_nodes = []
        var_nodes = []
        num_nodes = []
        for i, chromosome_variant_to_nodes_file in enumerate(args.variant_to_nodes):
            chromosome_variant_to_nodes = VariantToNodes.from_file(chromosome_variant_to_nodes_file)

            chromosome_ref_nodes = chromosome_variant_to_nodes.ref_nodes
            chromosome_var_nodes = chromosome_variant_to_nodes.var_nodes
            chromosome_num_nodes = np.max([chromosome_ref_nodes, chromosome_var_nodes]) + 1

            ref_nodes.append(chromosome_ref_nodes + int(np.sum(num_nodes[:i])))
            var_nodes.append(chromosome_var_nodes + int(np.sum(num_nodes[:i])))
            num_nodes.append(chromosome_num_nodes)

        ref_nodes = np.concatenate(ref_nodes)
        var_nodes = np.concatenate(var_nodes)

        variant_to_nodes = VariantToNodes(ref_nodes, var_nodes)
        variant_to_nodes.to_file(args.output_prefix + '.variant_to_nodes.pkl')
        to_file(num_nodes, args.output_prefix + '.num_nodes.pkl')

    subparser = subparsers.add_parser("merge_chromosome_variant_to_nodes")
    subparser.add_argument("--variant-to-nodes", nargs='+', required=True)
    subparser.add_argument("--output-prefix", required=True)
    subparser.set_defaults(func=merge_chromosome_variant_to_nodes)

    def merge_chromosome_haplotype_to_nodes(args):
        from obgraph.haplotype_nodes import DiscBackedHaplotypeToNodes, DiscBackedRaggedArray

        num_nodes_per_chromosome = from_file(args.num_nodes)
        num_chromosomes = len(num_nodes_per_chromosome)
        num_haplotypes = from_file(args.haplotype_to_nodes[0]).n_haplotypes()

        assert len(args.haplotype_to_nodes) == num_chromosomes

        out_file_name = args.output_prefix + ".haplotype_to_nodes.pkl"
        haplotype_nodes_out_file_name = out_file_name  + ".haplotype_nodes"
        haplotype_nodes_out_file = open(haplotype_nodes_out_file_name, "wb")

        offsets = []
        lengths = []
        offset = 0

        for h in range(num_haplotypes):
            logging.info("Individual %d/%d" % (h + 1, num_haplotypes))
            offsets.append(offset)
            length = 0
            for c in range(num_chromosomes):
                nodes = DiscBackedHaplotypeToNodes.from_file(args.haplotype_to_nodes[c]).get_nodes(h) + int(np.sum(num_nodes_per_chromosome[:c]))
                print(nodes)
                haplotype_nodes_out_file.write(nodes.astype(np.int64))
                length += len(nodes)
            logging.info("%d nodes" % length)
            lengths.append(length)
            offset += length

        haplotype_nodes_out_file.close()

        haplotype_to_nodes = DiscBackedHaplotypeToNodes(DiscBackedRaggedArray(haplotype_nodes_out_file_name, np.array(offsets, dtype=np.uint32), np.array(lengths, dtype=np.uint32)))
        haplotype_to_nodes.to_file(out_file_name)

    subparser = subparsers.add_parser("merge_chromosome_haplotype_to_nodes")
    subparser.add_argument("--haplotype-to-nodes", nargs='+', required=True)
    subparser.add_argument("--num-nodes", required=True)
    subparser.add_argument("--output-prefix", required=True)
    subparser.set_defaults(func=merge_chromosome_haplotype_to_nodes)

    def merge_flat_kmers(args):
        from graph_kmer_index.flat_kmers import FlatKmers

        num_nodes = from_file(args.num_nodes) if args.num_nodes is not None else None
        chromosome_lengths = np.atleast_2d(np.loadtxt(args.reference_fasta_fai, dtype=str))[:, 1].astype(np.int64) if args.reference_fasta_fai is not None else None

        hashes = []
        nodes = []
        ref_offsets = []
        allele_frequencies = []
        for i, flat_kmers_file in enumerate(args.flat_kmers):
            flat_kmers = from_file(flat_kmers_file)
            hashes.extend(flat_kmers._hashes)
            chromosome_start_node = int(np.sum(num_nodes[:i])) if num_nodes is not None else 0
            nodes.extend(flat_kmers._nodes + chromosome_start_node)
            if flat_kmers._ref_offsets is not None:
                chromosome_ref_offset = np.sum(chromosome_lengths[:i]) if chromosome_lengths is not None else 0
                ref_offsets.extend(flat_kmers._ref_offsets + chromosome_ref_offset)

            allele_frequencies.extend(flat_kmers._allele_frequencies)

        if len(ref_offsets) == 0:
            ref_offsets = None
        else:
            ref_offsets = np.array(ref_offsets, np.uint64)

        merged_kmers = FlatKmers(np.array(hashes, dtype=np.uint64), np.array(nodes, np.uint32), ref_offsets,
                                 np.array(allele_frequencies, dtype=np.single))

        merged_kmers.to_file(args.out_file_name)
        logging.info("Wrote merged index to %s" % merged_kmers)

    subparser = subparsers.add_parser("merge_flat_kmers")
    subparser.add_argument("--flat-kmers", nargs='+', required=True)
    subparser.add_argument("--num-nodes", required=False)
    subparser.add_argument("--reference-fasta-fai", required=False)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.set_defaults(func=merge_flat_kmers)

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)


    args = parser.parse_args(args)
    args.func(args)


