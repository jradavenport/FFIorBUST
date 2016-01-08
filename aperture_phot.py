import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from astropy.wcs import WCS
from astropy.visualization import scale_image
import photutils
from astropy.stats import sigma_clipped_stats

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
    Creates an exposure class which has many extensions.
    
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
        
  def extension(self, extension_idx, threshold=3.0, FWHM='', sigma=3.0, plot=True, verbose=True):
    '''
    A method to run aperatue photometry routines on an individual extension and save the results to the exposure class
    
    Parameters
    ----------
    extension_idx: int
      Index of the extension
    threshold: float (optional)
      The absolute image value above which to select sources
    FWHM: float (optional)
      The full width at half maximum  
      

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

    # Estimate the background and background noise
    mean, median, std = sigma_clipped_stats(data, sigma=sigma, iters=5)

    # Subtract background and generate sources list of all detections
    sources = photutils.daofind(data-median, threshold, FWHM or 5.*std)
    self.source_lists.append(sources)     
    
    # Plot the sources
    if plot:
      positions = (sources['xcentroid'], sources['ycentroid'])
      apertures = CircularAperture(positions, r=4.)
      norm = ImageNormalize(stretch=SqrtStretch())
      plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
      apertures.plot(color='blue', lw=1.5, alpha=0.5)
      