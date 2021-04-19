
import pytest

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
