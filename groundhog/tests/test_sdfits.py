import pytest
import numpy as np

from groundhog import sd_fits


def test_remove_edge_channels(sd_fits_table_hi):
    table, head = sd_fits_table_hi
    sdfits = sd_fits.SDFITS(table, head)
    sdfits.remove_edge_chans()
    ps_scan = sdfits.get_scans(6, ifnum=0, intnum=0, plnum=0)
    assert ps_scan.data.shape == (1, 26215)
    
