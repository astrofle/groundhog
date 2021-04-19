
import pytest

from groundhog import sd_fits_io


def test_read_sdfits(sd_fits_table):
    table, head = sd_fits_table
    assert table.shape == (104,)
    

def test_write_sdfits(sd_fits_table, tmp_path):
    d = tmp_path / "sd_fits_io"
    d.mkdir()
    out = d / "test.fits"
    table, head = sd_fits_table
    sd_fits_io.write_sdfits(out, table, head, overwrite=False)
    assert len(list(tmp_path.iterdir())) == 1
