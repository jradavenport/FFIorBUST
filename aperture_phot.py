import numpy as np
import matplotlib.pyplot as plt
import photutils, os
import astropy.io.fits as pf
from glob import glob
from astropy.wcs import WCS
from astropy.visualization import scale_image
from astropy.stats import sigma_clipped_stats
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

'''
This is a module to extract light curves from all objects in 53 Kepler Full Frame Image exposures.

Example:

import aperture_phot as ap
ex = ap.exposure(filepath)        # Load the fits file into an exposure class
ex.extension(4)                   # Run DAOfind on the extension of the given index and add to the source_tables list
ex.source_tables                  # Print the list of source_tables 

'''

def light_curves():
  '''
  Plot and save all possible light curves from the Kepler Full Frame Images
  '''
  data = {}
  
  for filepath in glob('./data/*.fits'):
    
    # Create exposure class
    ex = exposure(filepath)
    
    # Get photometry for all extensions
    for idx in range(85): ex.extension(idx)
    
    # Add data to a master dictionary
    data[ex.date_str] = ex.source_tables
  
  # Do some stuff to match objects across exposures
  # CODE CODE CODE
  
  # Pull out all magnitudes for all exposures
  # CODE CODE CODE
  
  # Plot the 53 exposure light curves and save them
  plt.plot(datetime, magnitudes)
  plt.title(source_id)
  plt.savefig('./plots/{}.png'.format(source_id))

class exposure:
  def __init__(self, filepath, verbose=False):
    '''
    Creates an exposure class which has many extensions.
    
    Parameters
    ----------
    filepath: str
      Path to the fits file of a single exposure
    verbose: bool
      Print some info to visually inspect
    
    '''
    self.source_tables = []
    
    # Get the datetime of the exposure from the filename
    self.date_str = os.path.basename(filepath).replace('_ffi-cal.fits','').replace('kplr','')
    self.datetime = time.strptime(date_str, '%Y%j%H%M%S')
    
    # Open the file and print the info
    self.hdulist = pf.open(filepath)
    if verbose:
      self.hdulist.info()
        
  def extension(self, extension_idx, threshold='', FWHM=3.0, sigma=3.0, snr=50., plot=False):
    '''
    A method to run aperatue photometry routines on an individual extension and save the results to the exposure class
    
    Parameters
    ----------
    extension_idx: int
      Index of the extension
    threshold: float (optional)
      The absolute image value above which to select sources
    FWHM: float
      The full width at half maximum
    sigma: float
      Number of standard deviations to use for background estimation
    snr: float
      The signal-to-noise ratio to use in the threshold detection
    plot: bool
      Plot the field with identified sources circled      

    Returns
    -------
    source_list: table
      A source list for the image

    '''

    # Define the data array
    data = self.hdulist[extension_idx].data.astype(np.float)
    
    # Extract the data header and create a WCS object
    hdr = self.hdulist[extension_idx].header
    wcs = WCS(hdr)

    # Estimate the background and background noise
    mean, median, std = sigma_clipped_stats(data, sigma=sigma, iters=5)

    # Calculate the detection threshold and FWHM if not provided
    if not threshold: threshold = np.mean(photutils.detect_threshold(data, snr=snr))
    
    # Print the parameters being used
    for p,v in zip(['mean','median','std','threshold','FWHM'],[mean,median,std,threshold,FWHM]): print '{!s:10}: {:.3f}'.format(p,v)

    # Subtract background and generate sources list of all detections
    sources = photutils.daofind(data-median, threshold, FWHM)
    self.source_tables.append(sources)  
    
    # Plot the sources
    if plot:
      positions = (sources['xcentroid'], sources['ycentroid'])
      apertures = photutils.CircularAperture(positions, r=4.)
      norm = ImageNormalize(stretch=SqrtStretch())
      plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
      apertures.plot(color='blue', lw=1.5, alpha=0.5)
    
    print '{!s:10}: {}'.format('sources',len(sources))
  