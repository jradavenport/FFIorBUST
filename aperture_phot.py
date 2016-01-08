import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from astropy.wcs import WCS
from astropy.visualization import scale_image
import photutils

'''
Example:

import aperture_phot as ap
ex = ap.exposure(filepath)        # Load the fits file into an exposure class
ex.extension(4)                   # Fun DAOfind on the extension of the given index and add to the source_lists variable
ex.source_lists                   # Print the list of source_lists 

'''

class exposure:
  def __init__(self, filepath, verbose=True):
    '''
    Parameters
    ----------
    filepath: str
      Path to the fits file of a single exposure
    verbose: bool
      Print some info to visually inspect
    
    '''
    self.source_lists = []
    
    # Open the file and print the info
    self.hdulist = pf.open(filepath)
    if verbose:
      self.hdulist.info()
        
  def extension(self, extension_idx, plot=True, verbose=True):
    '''
    Parameters
    ----------
    extension_idx: int
      Index of the extension

    Returns
    -------
    source_list: table
      A source list for the image

    '''

    # Define the data and error arrays
    data = self.hdulist[extension_idx].data.astype(np.float)

    # Extract the data header and create a WCS object
    hdr = self.hdulist[extension_idx].header
    wcs = WCS(hdr)

    # Display the data
    if plot: plt.imshow(scale_image(data, scale='sqrt', percent=99.5))

    # Find the absolute image value above which to select sources and define FWHM
    threshold = 1
    FWHM = 2

    # Generate source_list of all detections
    source_list = photutils.daofind(data, threshold, FWHM)
    
    self.source_lists.append(source_list)         
    
      