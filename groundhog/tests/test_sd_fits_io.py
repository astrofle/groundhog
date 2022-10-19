
import pytest

from groundhog import sd_fits_io


def test_read_sdfits(sd_fits_table):
    sdfits = sd_fits_table
    assert (sdfits.hdu[1].read()).shape == (104,)
    

def test_write_sdfits(sd_fits_table, tmp_path):
    d = tmp_path / "sd_fits_io"
    d.mkdir()
    out = d / "test.fits"
    sd_fits_io.write_sdfits(out, sd_fits_table.hdu[1][:], overwrite=True)
    assert len(list(tmp_path.iterdir())) == 1
