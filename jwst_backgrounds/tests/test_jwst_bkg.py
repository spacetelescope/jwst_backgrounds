from astropy.io import ascii
from astropy.table import Table
import numpy as np

from jwst_backgrounds.jbt import background 


def test_backgrounds(thisday=100):
    """ Test backgrounds output from github example."""

    # Column Names
    column_names = ['wavelenght', 'total_background', 'zodical_background', 'nonzodical_background',
                    'stray_light_background', 'thermal_background']
    # Read truth files
    truth = ascii.read('background.txt', names=column_names)

    # Run background calculation
    bkg = background(261.6833333, -73.33222222, 2.15, thresh=1.1)
    calendar = bkg.bkg_data['calendar']
    thisday_index = np.where(thisday == calendar)[0][0]

    data = []

    for i,wavelength in enumerate(bkg.bkg_data['wave_array']):
        data.append([wavelength, bkg.bkg_data['total_bg'][thisday_index][i], 
                    bkg.bkg_data['zodi_bg'][thisday_index][i],bkg.bkg_data['nonzodi_bg'][i],
                    bkg.bkg_data['stray_light_bg'][thisday_index][i],bkg.bkg_data['thermal_bg'][i]])

    # Build background output table
    test_data = Table(rows=data, names=column_names)

    # Compare values of tables.
    assert truth.values_equal(test_data)