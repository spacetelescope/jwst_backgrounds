from os.path import abspath, dirname, join

from astropy.io import ascii
from astropy.table import Table
import numpy as np
import pytest

from jwst_backgrounds.jbt import background

TEST_DIR = abspath(join(dirname(__file__)))

@pytest.mark.parametrize("ra, dec, expected", [(82.82, -5.39, "1073/sl_pix_107381.bin"), (10.68471, 41.26875, "0330/sl_pix_033039.bin")])
def test_myfile_from_healpix(ra, dec, expected):
    bkg = background(ra, dec, 2.15)
    healpix_file = bkg.myfile_from_healpix(ra, dec)

    assert healpix_file == expected


def test_backgrounds(thisday=100):
    """ Test backgrounds output from github example."""

    # Column Names
    column_names = ['wavelenght', 'total_background', 'zodical_background', 'nonzodical_background',
                    'stray_light_background', 'thermal_background']
    # Read truth files
    data_file = join(TEST_DIR, 'background.txt')
    truth = ascii.read(data_file, names=column_names)

    # Run background calculation
    bkg = background(261.6833333, -73.33222222, 2.15, thresh=1.1)
    calendar = bkg.bkg_data['calendar']
    thisday_index = np.where(thisday == calendar)[0][0]

    data = []

    for i, wavelength in enumerate(bkg.bkg_data['wave_array']):
        data.append([wavelength, bkg.bkg_data['total_bg'][thisday_index][i],
                    bkg.bkg_data['zodi_bg'][thisday_index][i], bkg.bkg_data['nonzodi_bg'][i],
                    bkg.bkg_data['stray_light_bg'][thisday_index][i], bkg.bkg_data['thermal_bg'][i]])

    # Build background output table
    test_data = Table(rows=data, names=column_names)

    # Compare values of tables.
    assert truth.values_equal(test_data)
