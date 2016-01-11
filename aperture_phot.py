import numpy as np
import matplotlib.pyplot as plt
import os, time
import astropy.io.fits as pf
from glob import glob
from photutils import daofind, aperture_photometry, detect_threshold, CircularAperture
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from astropy.visualization import scale_image
from astropy.stats import sigma_clipped_stats
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.table import vstack, Table, Column

'''
This is a module to extract light curves from all objects in 53 Kepler Full Frame Image exposures.

Authors: Joe Filippazzo, Brigitta Sipocz, Jim Davenport, Jennifer Cash (2016)

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
    
    # Get photometry for all 85 extensions
    for idx in range(5): ex.extension(idx)
    
    # Add data to a master dictionary
    data[ex.date_str] = ex.source_table
  
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
    self.source_table = Table(names=['aperture_sum','xcenter','ycenter','ra','dec'])
    
    # Open the file and print the info
    self.hdulist = pf.open(filepath)
    if verbose:
      self.hdulist.info()
    
    # Get the datetime of the exposure from the filename
    self.date_str = os.path.basename(filepath).replace('_ffi-cal.fits','').replace('kplr','')
    self.datetime = time.strptime(self.date_str, '%Y%j%H%M%S')
        
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
    
    # Extract the header and create a WCS object
    hdr = self.hdulist[extension_idx].header
    wcs = WCS(hdr)

    # Estimate the background and background noise
    mean, median, std = sigma_clipped_stats(data, sigma=sigma, iters=5)

    # Calculate the detection threshold and FWHM if not provided
    if not threshold: threshold = np.mean(detect_threshold(data, snr=snr))
    
    # Print the parameters being used
    for p,v in zip(['mean','median','std','threshold','FWHM'],[mean,median,std,threshold,FWHM]): print '{!s:10}: {:.3f}'.format(p,v)

    # Subtract background and generate sources list of all detections
    sources = daofind(data-median, threshold, FWHM)
    
    # Map RA and Dec to pixels
    positions = (sources['xcentroid'], sources['ycentroid'])
    skycoords = pixel_to_skycoord(*positions, wcs=wcs)
    
    # Calculate magnitudes at given source positions
    apertures = CircularAperture(positions, r=2.)
    photometry_table = aperture_photometry(data, apertures)
    
    # 'skycoords' IRCS object is problematic for stacking tables so for now we'll just add the ra and dec
    # photometry_table['sky_center'] = skycoords
    photometry_table['ra'], photometry_table['dec'] = skycoords.ra, skycoords.dec
    
    # Update data in the exposure object
    self.source_table = vstack([self.source_table,photometry_table], join_type='inner')  
    
    # Plot the sources
    if plot:
      norm = ImageNormalize(stretch=SqrtStretch())
      plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
      apertures.plot(color='blue', lw=1.5, alpha=0.5)
    
    print '{!s:10}: {}'.format('sources',len(sources))
  