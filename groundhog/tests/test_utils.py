
import pytest

from groundhog import utils


def test_ruze():
    g0 = 0.71
    eta = utils.ruze(21e-2, g0, 0.)
    assert eta == g0
