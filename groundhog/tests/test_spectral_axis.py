
import pytest
import numpy as np

from groundhog import sd_fits_io
from groundhog import spectral_axis


def test_compute_freq_axis(sd_fits_table, gbtidl_spec):
    table = sd_fits_table.hdu[1].read()
    freq = spectral_axis.compute_spectral_axis(table)
    np.testing.assert_allclose(freq.to('MHz').value[0], gbtidl_spec[:,0])
    
