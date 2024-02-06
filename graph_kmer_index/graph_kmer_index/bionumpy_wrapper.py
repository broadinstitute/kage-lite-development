import bionumpy as bnp

def bionumpy_hash(numeric_sequence, k):
    numeric_sequence = numeric_sequence
    encoded_sequence = bnp.EncodedArray(numeric_sequence, bnp.encodings.alphabet_encoding.ACTGEncoding)
    kmers = bnp.sequence.get_kmers(encoded_sequence, k).raw()
    return kmers





