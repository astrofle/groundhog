
import os
import pytest
import numpy as np


from groundhog import sd_fits_io
from groundhog import sd_fits_utils


#def test_make_summary(sd_fits_table):
#    summary = sd_fits_utils.make_summary(sd_fits_table.table[0])
#    assert summary.Scan == [5, 6]
#    assert summary.Proc == ['OffOn', 'OffOn']
#    assert summary.Az == [247.89581746345755, 253.10093738847564]
    

#def test_update_table_column(sd_fits_table):
#    sdfits = sd_fits_table
#    table = sdfits.hdu[1][:]
#    head = sdfits.hdu[1].read_header()
#    new_array = np.zeros((len(table),15000), dtype=float)
#    new_table = sd_fits_utils.update_table_column(table, 'DATA', new_array)
#    np.testing.assert_array_equal(new_table['DATA'], new_array)
#    new_table = sd_fits_utils.update_table_column(new_table, 'IFNUM', np.ones(len(table)))
#    new_table = sd_fits_utils.update_table_column(new_table, 'TCAL', np.ones((len(table),15000)))
#    sd_fits_io._write_sdfits('test.fits', new_table, head, overwrite=False)
#    
#    # Clean up
#    os.remove('test.fits')
