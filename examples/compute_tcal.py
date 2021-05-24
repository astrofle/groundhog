"""
Example code for computing the temperature of the noise diode (Tcal)
from On-Off observations of a source of known flux density.
"""

import numpy as np

# Import groundhog.
from groundhog import datared
from groundhog import sd_fits_io, sd_fits, sd_fits_utils


# Load the data.
sdfits = sd_fits_io.read_sdfits('../groundhog/tests/data/TSCAL_210420_PF.raw.vegas.fits')

# Compute Tcal for the first spectral window and polarization,
# and tie the calibration to the Perley and Butler 2017 flux scale.
tcal = datared.get_tcal(sdfits, 5, ifnum=0, plnum=0, scale="Perley-Butler 2017", units="K")

# If you want to get the frequency as well.
cal_scan = sdfits.get_scans(6, plnum=0, ifnum=[0])
# Average the integrations in the scan to get a single frequency vector.
cal_scan.average()
# The frequency vector lives here.
freq = cal_scan.freq

# Udate the SDFITS with the new Tcal values.
sdfits.update_table_col('TCAL', np.tile(tcal.value, (len(sdfits.table),1)))
# Create a new table with the computed vector of Tcal values.
#new_table = sd_fits_utils.update_table_column(table, 'TCAL', np.tile(tcal.value, (len(table),1)))

# Write to a new SDFITS.
sd_fits_io.write_sdfits("TSCAL_210420_PF.updated.vegas.fits", sdfits, overwrite=False)
