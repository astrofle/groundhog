"""
Example code for computing the temperature of the noise diode (Tcal)
from On-Off observations of a source of known flux density.
"""

# Import groundhog.
from groundhog import datared
from groundhog import sd_fits_io, sd_fits


# Load the data.
table, head = sd_fits_io.read_sdfits('../groundhog/tests/data/TSCAL_210420_PF.raw.vegas.fits')
sdfits = sd_fits.SDFITS(table, head)

# Compute Tcal for the first spectral window and polarization,
# and tie the calibration to the Perley and Butler 2017 flux scale.
tcal = datared.get_tcal(sdfits, 5, ifnum=0, plnum=0, scale="Perley-Butler 2017", units="K")

# If you want to get the frequency as well.
cal_scan = sdfits.get_scans(6, plnum=0, ifnum=[0])
# Average the integrations in the scan to get a single frequency vector.
cal_scan.average()
# The frequency vector lives here.
freq = cal_scan.freq
