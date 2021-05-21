import pytest
import numpy as np

from groundhog import sd_fits


def test_remove_edge_channels(sd_fits_table_hi):
    table, head = sd_fits_table_hi
    sdfits = sd_fits.SDFITS(table, head)
    # First get the frequency without removing edges.
    ps_scan = sdfits.get_scans(6, ifnum=0, intnum=0, plnum=0)
    freq_tot = ps_scan.freq
    spec_tot = ps_scan.data
    # Now remove the edges.
    sdfits.remove_edge_chans()
    ps_scan = sdfits.get_scans(6, ifnum=0, intnum=0, plnum=0)
    # Does it have the correct shape?
    assert ps_scan.data.shape == (1, 26215)
    freq_crop = ps_scan.freq
    spec_crop = ps_scan.data
    # Compare.
    idx0 = np.argmin(abs(freq_tot[0] - freq_crop[0,0]))
    idxf = np.argmin(abs(freq_tot[0] - freq_crop[0,-1]))
    idx0,idxf = np.sort([idx0,idxf+1])
    np.testing.assert_allclose(freq_tot[0][idx0:idxf], freq_crop[0])
    np.testing.assert_allclose(spec_tot[0][idx0:idxf], spec_crop[0])
