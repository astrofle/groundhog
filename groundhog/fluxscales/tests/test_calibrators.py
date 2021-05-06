
import pytest
import numpy as np

from astropy import units as u

from groundhog.fluxscales import calibrators


def test_compute_sed():
    freq = 1*u.GHz
    scale = 'Perley-Butler 2017'
    source = '3C286'
    snu = calibrators.compute_sed(freq, scale, source, units='Jy')
    ta = calibrators.compute_sed(freq, scale, source, units='K')
    pytest.raises(KeyError, calibrators.compute_sed, *{'freq':freq, 'scale':'Maddalena-Frayer 2020', 'source':source})
    pytest.raises(KeyError, calibrators.compute_sed, *{'freq':freq, 'scale':scale, 'source':'ZwI20'})


def test_compute_sed2():
    scale = 'Perley-Butler 2017'
    source = '3C295'
    expected = np.loadtxt('data/{}.csv'.format(source), delimiter=',')
    freq = expected[:,0]*u.GHz
    flux = expected[:,1]*u.Jy
    snu = calibrators.compute_sed(freq, scale, source, units='Jy')
    np.testing.assert_allclose(snu, flux, atol=2)
