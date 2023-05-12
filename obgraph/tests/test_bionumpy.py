import pytest
import time
import numpy as np
import bionumpy as bnp


@pytest.mark.skip
def test():
    t = time.perf_counter()
    file = bnp.open("variants.vcf", buffer_type=bnp.PhasedVCFMatrixBuffer)
    n = 0
    chunk = file.read_chunk()
    genotypes = chunk.genotypes.raw()
    print(np.unique(genotypes))
    n += chunk.genotypes.shape[0]
    print(n)
    print(time.perf_counter()-t)

test()