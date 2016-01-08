import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from astropy.wcs import WCS
from astropy.visualization import scale_image

class get_exposure:
  def __init__(self, filepath, verbose=True):
    '''
    Parameters
    ----------
    filepath: str
      Path to the fits file of a single exposure
    verbose: bool
      Print some info to visually inspect
    
    '''
    
    # Open the file and print the info
    self.hdulist = pf.open(filepath)
    if verbose:
      self.hdulist.info()
    
  def get_extension(self, extension_idx, plot=True, verbose=True):
    '''
    Parameters
    ----------
    extension_idx: int
      Index of the extension
  
    Returns
    -------
    Nothing yet
  
    '''

    # Define the data and error arrays
    data = self.hdulist[extension_idx].data.astype(np.float)

    # Extract the data header and create a WCS object
    hdr = self.hdulist[extension_idx].header
    wcs = WCS(hdr)

    # Display the data
    if plot: plt.imshow(scale_image(data, scale='sqrt', percent=99.5))
  
    return hdr