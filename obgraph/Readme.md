# Obgraph
Obgraph is a simple Python package for creating sequence-graph from a reference and a set of variants.

An Obgraph is a graph with "dummy" nodes for all deleted sequences, meaning that every insertion and deletion is represented with two nodes: One empty node and one node for the sequence (inserted sequence or reference sequence in case of a deletion). This means that any variant always can be represented with two nodes, which is a format suitable for e.g. genotyping. 

## Installation
1) Clone the repository
2) `cd obgraph`
3) `python3 -m pip install .`

## Typical usage

### Make a graph for a single chromosome
```bash
obgraph make -r ref.fa -c [chromosome] -v variants.vcf -o mygraph
```

### Add allele frequencies to the graph
This changes the graph inplace, and does not create a new graph.
```bash
obgraph add_allele_frequencies -c [chromosome] -v variants.vcf -g mygraph
```

### Merge graphs from multiple chromosomes
Note: Make sure graphs are sent as input in correct chromosome order.
```bash
obgraph merge_graphs -o whole_genome -g graph_chr1 graph_chr2 ...
```

### Other functionality
This package can also make a genotype matrix from a vcf file containing genotype information. The matrix is a basically a numpy matrix with variants as columns and individuals as rows, and can be used to easily query individuals with given genotypes at variant sites.

```bash
obgraph make_genotype_matrix -t N_THREADS -g graph -v variants.vcf -o genotype_matrix -n N_INDIVIDUALS -m N_VARIANTS
```

### Using graphs
```python
from obgraph.graph import Graph
graph = Graph.from_file("graph")

# Get edges
edges = graph.ged_edges(1)  # node id 1, note that all node IDs are numeric
```

