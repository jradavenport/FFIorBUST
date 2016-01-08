import numpy as np
import matplotlib.pylab as plt

from astropy.io import fits
from astropy.stats import mad_std
import glob
from astropy.wcs.utils import pixel_to_skycoord
from astropy.wcs import WCS
from astropy.utils.misc import isiterable

from matplotlib.colors import LogNorm

from photutils import daofind, aperture_photometry, CircularAperture


def do_photometry(hdu, extensions=None, threshold=5, fwhm=2.5):

    if extensions is None:
        extensions = np.arange(1, len(hdu) + 1)
    if not isiterable(extensions):
        extensions = (extensions, )

    output = {}
    for ext in extensions:
        header = hdu[ext].header
        data = hdu[ext].data
        image_wcs = WCS(header)

        background = mad_std(data)

        sources = daofind(data, threshold=threshold * background, fwhm=fwhm)
        positions = (sources['xcentroid'], sources['ycentroid'])
        sky_positions = pixel_to_skycoord(*positions, wcs=image_wcs)

        apertures = CircularAperture(positions, r=2.)
        photometry_table = aperture_photometry(data, apertures)
        photometry_table['sky_center'] = sky_positions

        output[str(ext)] = photometry_table

    return output


input_files = glob.glob('data/kplr*ffi-cal.fits')

for ffi in input_files:
    fits.open(ffi)
    do_photometry(ffi, extentions=1)

