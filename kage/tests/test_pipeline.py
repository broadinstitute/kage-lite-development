import pytest
import obgraph.command_line_interface as obgraph_cli
import graph_kmer_index.command_line_interface as graph_kmer_index_cli
import kage.command_line_interface as kage_cli
import lite_utils.command_line_interface as lite_utils_cli
import kmer_mapper.kmer_mapper.command_line_interface as kmer_mapper_cli

test_resources_dir = 'resources'
output_dir = '/tmp'

# tests replicate individual tasks or portions of the WDL workflow using the same test data
# expected results are not currently checked

def test_MakeSitesOnlyVcfAndNumpyVariants():
    vcf = f'{test_resources_dir}/MakeSitesOnlyVcfAndNumpyVariants/inputs/test.sites.vcf'
    output = f'{output_dir}/test.numpy_variants.pkl'

    obgraph_cli.run_argument_parser([
        'make_numpy_variants',
        '-v', vcf,
        '-o', output])

def test_MakeChromosomeGenotypeMatrix():
    vcf_gz = f'{test_resources_dir}/MakeChromosomeGenotypeMatrix/inputs/test.preprocessed.bi.vcf.gz'
    chromosome = 'chr1'
    output = f'{output_dir}/test.{chromosome}.genotype_matrix.pkl'

    lite_utils_cli.run_argument_parser([
        'make_genotype_matrix',
        '-v', vcf_gz,
        '-c', chromosome,
        '-o', output])

def test_MakeChromosomeGraph():
    vcf_gz = f'{test_resources_dir}/MakeChromosomeGraph/inputs/test.preprocessed.bi.vcf.gz'
    reference_fasta = f'{test_resources_dir}/MakeChromosomeGraph/inputs/GRCh38_full_analysis_set_plus_decoy_hla.chr1-1Mbp-chr2-1Mbp.fa'
    chromosome = 'chr1'
    output_prefix = f'{output_dir}/test'
    output_graph = f'{output_prefix}.{chromosome}.obgraph.pkl'
    output_position_id_index = f'{output_prefix}.{chromosome}.position_id_index.pkl'

    obgraph_cli.run_argument_parser([
        'make',
        '-v', vcf_gz,
        '-c', chromosome,
        '-r', reference_fasta,
        '-o', output_graph])
    obgraph_cli.run_argument_parser([
        'add_allele_frequencies',
        '-v', vcf_gz,
        '-c', chromosome,
        '-g', output_graph])
    obgraph_cli.run_argument_parser([
        'make_position_id',
        '-g', output_graph,
        '-o', output_position_id_index])

def test_MakeChromosomeVariantToNodes():
    graph = f'{test_resources_dir}/MakeChromosomeVariantToNodes/inputs/test.chr1.obgraph.pkl'
    vcf_gz = f'{test_resources_dir}/MakeChromosomeVariantToNodes/inputs/chr1.sites.vcf.gz'
    chromosome = 'chr1'
    output = f'{output_dir}/test.{chromosome}.variant_to_nodes.pkl'

    obgraph_cli.run_argument_parser([
        'make_variant_to_nodes',
        '-g', graph,
        '-v', vcf_gz,
        '-o', output])

def test_MakeChromosomeHaplotypeToNodes():
        variant_to_nodes = f'{test_resources_dir}/MakeChromosomeHaplotypeToNodes/inputs/test.chr1.variant_to_nodes.pkl'
        genotype_matrix = f'{test_resources_dir}/MakeChromosomeHaplotypeToNodes/inputs/test.chr1.genotype_matrix.pkl'
        chromosome = 'chr1'
        output = f'{output_dir}/test.{chromosome}.haplotype_to_nodes.pkl'

        obgraph_cli.run_argument_parser([
            'make_haplotype_to_nodes_bnp',
            '-d', 'True',
            '-g', variant_to_nodes,
            '-v', genotype_matrix,
            '-o', output])

def test_MergeChromosomeGraphs():
    chr1_graph = f'{test_resources_dir}/MergeChromosomeGraphs/inputs/test.chr1.obgraph.pkl'
    chr2_graph = f'{test_resources_dir}/MergeChromosomeGraphs/inputs/test.chr2.obgraph.pkl'
    output_graph = f'{output_dir}/test.obgraph.pkl'
    output_position_id_index = f'{output_dir}/test.position_id_index.pkl'

    obgraph_cli.run_argument_parser([
        'merge_graphs',
        '-g', chr1_graph, chr2_graph,
        '-o', output_graph])
    obgraph_cli.run_argument_parser([
        'make_position_id',
        '-g', output_graph,
        '-o', output_position_id_index])

def test_MergeChromosomeVariantToNodes():
    chr1_variant_to_nodes = f'{test_resources_dir}/MergeChromosomeVariantToNodes/inputs/test.chr1.variant_to_nodes.pkl'
    chr2_variant_to_nodes = f'{test_resources_dir}/MergeChromosomeVariantToNodes/inputs/test.chr2.variant_to_nodes.pkl'
    output_prefix = f'{output_dir}/test'

    lite_utils_cli.run_argument_parser([
        'merge_chromosome_variant_to_nodes',
        '--variant-to-nodes', chr1_variant_to_nodes, chr2_variant_to_nodes,
        '--output-prefix', output_prefix])

def test_MergeChromosomeHaplotypeToNodes():
    chr1_haplotype_to_nodes = f'{test_resources_dir}/MergeChromosomeHaplotypeToNodes/inputs/test.chr1.haplotype_to_nodes.pkl'
    chr2_haplotype_to_nodes = f'{test_resources_dir}/MergeChromosomeHaplotypeToNodes/inputs/test.chr2.haplotype_to_nodes.pkl'
    num_nodes = f'{test_resources_dir}/MergeChromosomeHaplotypeToNodes/inputs/test.num_nodes.pkl'
    output_prefix = f'{output_dir}/test'

    lite_utils_cli.run_argument_parser([
        'merge_chromosome_haplotype_to_nodes',
        '--haplotype-to-nodes', chr1_haplotype_to_nodes, chr2_haplotype_to_nodes,
        '--num-nodes', num_nodes,
        '--output-prefix', output_prefix])

@pytest.mark.parametrize("num_threads", [2]) # single thread has a bug
def test_MakeHelperModel(num_threads):
    chr1_genotype_matrix = f'{test_resources_dir}/MakeHelperModel/inputs/test.chr1.genotype_matrix.pkl'
    chr2_genotype_matrix = f'{test_resources_dir}/MakeHelperModel/inputs/test.chr2.genotype_matrix.pkl'
    variant_to_nodes = f'{test_resources_dir}/MakeHelperModel/inputs/test.variant_to_nodes.pkl'
    window_size = '100'
    num_threads = str(num_threads)
    output_genotype_matrix = f'{output_dir}/test.genotype_matrix.pkl'
    output_helper_model = f'{output_dir}/test.helper_model.pkl'

    lite_utils_cli.run_argument_parser([
        'merge_genotype_matrices_and_convert_to_unphased',
        '-g', chr1_genotype_matrix, chr2_genotype_matrix,
        '-o', output_genotype_matrix])
    kage_cli.run_argument_parser([
        'create_helper_model',
        '-g', output_genotype_matrix,
        '-v', variant_to_nodes,
        '-w', window_size,
        '-t', num_threads,
        '-o', output_helper_model])

@pytest.mark.parametrize("num_threads", [1, 2])
def test_SampleChromosomeKmersFromLinearReference(num_threads):
    chromosome_reference_fasta = f'{test_resources_dir}/SampleChromosomeKmersFromLinearReference/inputs/chromosome.fa'
    num_threads = str(num_threads)
    spacing = '1'
    kmer_length = '31'
    include_reverse_complement = 'True'
    chromosome = 'chr1'
    chromosome_length = '1000000'
    output = f'{output_dir}/test.{chromosome}.linear_kmers.pkl'

    graph_kmer_index_cli.run_argument_parser([
        'make',
        '-t', num_threads,
        '-s', spacing,
        '-k', kmer_length,
        '--include-reverse-complement', include_reverse_complement,
        '-R', chromosome_reference_fasta,
        '-n', chromosome,
        '-G', chromosome_length,
        '-o', output])

def test_MergeFlatKmers(): # use test data from aliased task MergeChromosomeKmersFromLinearReference
    chr1_flat_kmers = f'{test_resources_dir}/MergeChromosomeKmersFromLinearReference/inputs/test.chr1.linear_kmers.pkl'
    chr2_flat_kmers = f'{test_resources_dir}/MergeChromosomeKmersFromLinearReference/inputs/test.chr2.linear_kmers.pkl'
    num_nodes = f'{test_resources_dir}/MergeChromosomeKmersFromLinearReference/inputs/test.num_nodes.pkl'
    reference_fasta_fai = f'{test_resources_dir}/MergeChromosomeKmersFromLinearReference/inputs/GRCh38_full_analysis_set_plus_decoy_hla.chr1-1Mbp-chr2-1Mbp.fa.fai'
    output = f'{output_dir}/test.linear_kmers.pkl'

    lite_utils_cli.run_argument_parser([
        'merge_flat_kmers',
        '--flat-kmers', chr1_flat_kmers, chr2_flat_kmers,
        '--num-nodes', num_nodes,
        '--reference-fasta-fai', reference_fasta_fai,
        '-o', output])

def test_MakeLinearReferenceKmerCounter():
    linear_kmers = f'{test_resources_dir}/MakeLinearReferenceKmerCounter/inputs/test.linear_kmers.pkl'
    subsample_ratio = '1'
    output = f'{output_dir}/test.linear_kmers_counter.pkl'

    graph_kmer_index_cli.run_argument_parser([
        'count_kmers',
        '--subsample-ratio', subsample_ratio,
        '-f', linear_kmers,
        '-o', output])

@pytest.mark.parametrize("num_threads", [1, 2])
def test_GetChromosomeShortVariantKmers(num_threads):
        graph = f'{test_resources_dir}/GetChromosomeShortVariantKmers/inputs/test.chr1.obgraph.pkl'
        position_id_index = f'{test_resources_dir}/GetChromosomeShortVariantKmers/inputs/test.chr1.position_id_index.pkl'
        variant_to_nodes = f'{test_resources_dir}/GetChromosomeShortVariantKmers/inputs/test.chr1.variant_to_nodes.pkl'
        linear_kmers_counter = f'{test_resources_dir}/GetChromosomeShortVariantKmers/inputs/test.linear_kmers_counter.pkl'
        vcf_gz = f'{test_resources_dir}/GetChromosomeShortVariantKmers/inputs/chr1.sites.vcf.gz'
        use_dense_kmer_finder = 'True'
        num_threads = str(num_threads)
        kmer_length = '31'
        chunk_size = '20000'
        max_variant_nodes = '3'
        chromosome = 'chr1'
        output = f'{output_dir}/test.{chromosome}.short_variant_kmers.pkl'

        graph_kmer_index_cli.run_argument_parser([
            'make_unique_variant_kmers',
            '-D', use_dense_kmer_finder,
            '-t', num_threads,
            '-k', kmer_length,
            '-c', chunk_size,
            '--max-variant-nodes', max_variant_nodes,
            '-g', graph,
            '-p', position_id_index,
            '-V', variant_to_nodes,
            '-I', linear_kmers_counter,
            '-v', vcf_gz,
            '-o', output])

def test_SampleChromosomeStructuralVariantKmers():
    graph = f'{test_resources_dir}/SampleChromosomeStructuralVariantKmers/inputs/test.chr1.obgraph.pkl'
    variant_to_nodes = f'{test_resources_dir}/SampleChromosomeStructuralVariantKmers/inputs/test.chr1.variant_to_nodes.pkl'
    linear_kmers_counter = f'{test_resources_dir}/SampleChromosomeStructuralVariantKmers/inputs/test.linear_kmers_counter.pkl'
    kmer_length = '31'
    chromosome = 'chr1'
    output = f'{output_dir}/test.{chromosome}.structural_variant_kmers.pkl'

    graph_kmer_index_cli.run_argument_parser([
        'sample_kmers_from_structural_variants',
        '-k', kmer_length,
        '-g', graph,
        '-V', variant_to_nodes,
        '-i', linear_kmers_counter,
        '-o', output])

def test_MakeReverseVariantKmerIndex():
    variant_kmers = f'{test_resources_dir}/MakeReverseVariantKmerIndex/inputs/test.variant_kmers.pkl'
    output = f'{output_dir}/test.reverse_variant_kmers.pkl'

    graph_kmer_index_cli.run_argument_parser([
        'make_reverse',
        '-f', variant_kmers,
        '-o', output])

def test_MakeVariantKmerIndexWithReverseComplements():
    variant_kmers = f'{test_resources_dir}/MakeVariantKmerIndexWithReverseComplements/inputs/test.variant_kmers.pkl'
    kmer_length = '31'
    hash_modulo = '200000033'
    add_reverse_complements = 'True'
    skip_frequences = 'True'
    output = f'{output_dir}/test.kmer_index_only_variants_with_revcomp.pkl'

    graph_kmer_index_cli.run_argument_parser([
        'make_from_flat',
        '-f', variant_kmers,
        '--kmer-size', kmer_length,
        '--hash-modulo', hash_modulo,
        '--add-reverse-complements', add_reverse_complements,
        '--skip-frequencies', skip_frequences,
        '-o', output])

@pytest.mark.parametrize("num_threads", [1, 2])
def test_MakeCountModel(num_threads):
    graph = f'{test_resources_dir}/MakeCountModel/inputs/test.obgraph.pkl'
    haplotype_to_nodes = f'{test_resources_dir}/MakeCountModel/inputs/test.haplotype_to_nodes.pkl'
    kmer_index_only_variants_with_revcomp = f'{test_resources_dir}/MakeCountModel/inputs/test.kmer_index_only_variants_with_revcomp.pkl'
    num_threads = str(num_threads)
    output = f'{output_dir}/test.sampling_count_model.pkl'

    kage_cli.run_argument_parser([
        'sample_node_counts_from_population',
        '-t', num_threads,
        '-g', graph,
        '-H', haplotype_to_nodes,
        '-i', kmer_index_only_variants_with_revcomp,
        '-o', output])

def test_RefineCountModel():
    sampling_count_model = f'{test_resources_dir}/RefineCountModel/inputs/test.sampling_count_model.pkl'
    variant_to_nodes = f'{test_resources_dir}/RefineCountModel/inputs/test.variant_to_nodes.pkl'
    output = f'{output_dir}/test.refined_sampling_count_model.pkl'

    kage_cli.run_argument_parser([
        'refine_sampling_model',
        '-s', sampling_count_model,
        '-v', variant_to_nodes,
        '-o', output])

def test_FindTrickyVariants():
    sampling_count_model = f'{test_resources_dir}/FindTrickyVariants/inputs/test.sampling_count_model.pkl'
    variant_to_nodes = f'{test_resources_dir}/FindTrickyVariants/inputs/test.variant_to_nodes.pkl'
    reverse_variant_kmers = f'{test_resources_dir}/FindTrickyVariants/inputs/test.reverse_variant_kmers.pkl'
    max_counts = '1000'
    output = f'{output_dir}/test.tricky_variants.pkl'

    kage_cli.run_argument_parser([
        'find_tricky_variants',
        '-m', sampling_count_model,
        '-v', variant_to_nodes,
        '-r', reverse_variant_kmers,
        '-M', max_counts,
        '-o', output])

def test_MakeIndexBundle():
    variant_to_nodes = f'{test_resources_dir}/MakeIndexBundle/inputs/test.variant_to_nodes.pkl'
    numpy_variants = f'{test_resources_dir}/MakeIndexBundle/inputs/test.numpy_variants.pkl'
    refined_sampling_count_model = f'{test_resources_dir}/MakeIndexBundle/inputs/test.refined_sampling_count_model.pkl'
    tricky_variants = f'{test_resources_dir}/MakeIndexBundle/inputs/test.tricky_variants.pkl'
    helper_model = f'{test_resources_dir}/MakeIndexBundle/inputs/test.helper_model.pkl'
    helper_model_combo_matrix = f'{test_resources_dir}/MakeIndexBundle/inputs/test.helper_model_combo_matrix.pkl'
    kmer_index_only_variants_with_revcomp = f'{test_resources_dir}/MakeIndexBundle/inputs/test.kmer_index_only_variants_with_revcomp.pkl'
    output = f'{output_dir}/test.index.pkl'

    kage_cli.run_argument_parser([
        'make_index_bundle',
        '-g', variant_to_nodes,
        '-v', numpy_variants,
        '-A', refined_sampling_count_model,
        '-x', tricky_variants,
        '-f', helper_model,
        '-F', helper_model_combo_matrix,
        '-i', kmer_index_only_variants_with_revcomp,
        '-o', output])

@pytest.mark.parametrize("num_threads", [1, 2])
def test_Case(num_threads):
    kmer_index_only_variants_with_revcomp = f'{test_resources_dir}/Case/inputs/test.kmer_index_only_variants_with_revcomp.pkl'
    index = f'{test_resources_dir}/Case/inputs/test.index.pkl'
    fasta = f'{test_resources_dir}/Case/inputs/HG00731.final.chr1-1Mbp-chr2-1Mbp.noN.fasta'
    num_threads = str(num_threads)
    chunk_size = '100000000'
    sample_name = 'HG00731'
    average_coverage = '30'
    ignore_helper_model = 'false'
    output_kmer_counts = f'{output_dir}/HG00731.kmer_counts.npy'
    output_vcf = f'{output_dir}/HG00731.vcf'

    kmer_mapper_cli.run_argument_parser([
        'map',
        '-t', num_threads,
        '-c', chunk_size,
        '-i', kmer_index_only_variants_with_revcomp,
        '-f', fasta,
        '-o', output_kmer_counts])
    kage_cli.run_argument_parser([
        'genotype',
        '-c', output_kmer_counts,
        '-i', index,
        '-s', sample_name,
        '--average-coverage', average_coverage,
        '-I', ignore_helper_model,
        '-o', output_vcf])
