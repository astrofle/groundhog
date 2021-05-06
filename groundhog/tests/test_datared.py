import pytest
import numpy as np

from groundhog import sd_fits
from groundhog import datared


def test_get_ps(sd_fits_table, sd_fits_table_hi, gbtidl_spec, gbtidl_spec_hi):
    table, head = sd_fits_table
    sdfits = sd_fits.SDFITS(table, head)
    tsou = datared.get_ps(sdfits, 5, ifnum=4, intnum=1, plnum=0, method='classic')
    np.testing.assert_allclose(tsou, gbtidl_spec[:,1], rtol=1e-4)
    
    table, head = sd_fits_table_hi
    sdfits = sd_fits.SDFITS(table, head)
    tsou = datared.get_ps(sdfits, 6, plnum=0, method='classic')
    np.testing.assert_allclose(tsou, gbtidl_spec_hi[:,1], rtol=0.02)
