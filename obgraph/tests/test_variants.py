from obgraph.graph import Graph
from obgraph.variants import VcfVariants, VcfVariant

def test_find_insertion_nodes():
    g = Graph.from_dicts(
        {
            1: "CTACCA",
            2: "AA",
            3: "TAAATAA",
            4: ""
        },
        {
            1: [2, 4],
            2: [3],
            4: [3]

        },
        [1, 3]
    )
    variant = VcfVariant(1, 6, "A", "AAA", "", "INSERTION")


    ref_node, variant_node = g.get_variant_nodes(variant)
    assert ref_node == 4
    assert variant_node == 2

test_find_insertion_nodes()

def test_position_index():
    variants = VcfVariants(
        [
            VcfVariant(1, 4),
            VcfVariant(1, 7),
            VcfVariant(1, 8),
            VcfVariant(1, 8),
            VcfVariant(3, 8),
            VcfVariant(3, 9),
            VcfVariant(3, 11),
        ]
    )

    variants.make_position_index()

    v = variants.get_variants_in_region(1, 4, 8)
    print(v)


test_position_index()
