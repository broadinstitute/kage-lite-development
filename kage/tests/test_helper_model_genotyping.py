import logging

from kage.indexing.index_bundle import IndexBundle
from kage.models.helper_model import HelperVariants, CombinationMatrix

logging.basicConfig(level=logging.DEBUG)
import numpy as np
from kage.genotyping.combination_model_genotyper import CombinationModelGenotyper
from obgraph.variant_to_nodes import VariantToNodes
from obgraph.variants import VcfVariants
from kage.models.sampling_combo_model import LimitedFrequencySamplingComboModel

def test():
    variant_to_nodes = VariantToNodes(np.array([0, 2]), np.array([1, 3]))
    node_counts = np.array([4, 3, 10, 0])

    combo_matrix = np.array([
        np.zeros((3, 3)),
        [[0.20867266, 0.00922959, 0.00097657],
         [0.00016128, 0.43460627, 0.01792008],
        [0.00022378, 0.00053976, 0.32767001]]
        ])


    helpers = np.array([1, 0])
    model_ref = LimitedFrequencySamplingComboModel.create_naive(2)
    model_var = LimitedFrequencySamplingComboModel.create_naive(2)

    index = IndexBundle({
        "count_model": [model_ref, model_var],
        "helper_variants": HelperVariants(helpers),
        "combination_matrix": CombinationMatrix(combo_matrix),
        "variant_to_nodes": variant_to_nodes
    })

    genotyper = CombinationModelGenotyper(
         0, 3, node_counts,
        index
    )

    genotypes, probs, count_probs = genotyper.genotype()
    genotypes = VcfVariants.translate_numeric_genotypes_to_literal(genotypes)


test()
