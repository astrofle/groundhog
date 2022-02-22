import pytest
import numpy as np

from groundhog import sd_fits
from groundhog import sd_fits_io


# Arrange
@pytest.fixture(scope="module")
def sd_fits_table():
    sdfits = sd_fits.SDFITS()
    sdfits.load('data/AGBT19B_334_04_3C353.raw.vegas.A.fits')
    return sdfits
    #return sd_fits_io.read_sdfits('data/AGBT19B_334_04_3C353.raw.vegas.A.fits')

@pytest.fixture(scope="module")
def sd_fits_table_hi():
    sdfits = sd_fits.SDFITS()
    sdfits.load('data/TGBT20A_506_01.raw.vegas.A.fits')
    return sdfits
    #return sd_fits_io.read_sdfits('data/TGBT20A_506_01.raw.vegas.A.fits')

@pytest.fixture(scope="module")
def gbtidl_spec():
    return np.loadtxt("data/AGBT19B_332_04_scan5_ifnum4_intnum1_plnum0.ascii", skiprows=3)

@pytest.fixture(scope="module")
def gbtidl_spec_hi():
    return np.loadtxt("data/TGBT20A_506_01_scan6_ifnum0_plnum0.ascii", skiprows=3)
