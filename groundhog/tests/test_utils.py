
import pytest
import numpy as np

from groundhog import utils


def test_factors():
    facs = utils.factors(32768)
    expc = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 
                     1024, 2048, 4096, 8192, 16384, 32768])
    #assert facs == expc
    np.testing.assert_allclose(np.sort(facs), expc)
    

def test_ruze():
    g0 = 0.71
    eta = utils.ruze(21e-2, g0, 0.)
    assert eta == g0
