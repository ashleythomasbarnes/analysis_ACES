from glob import glob
from astropy.io import fits
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs
import numpy as np
from tqdm.auto import tqdm
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import astropy.units as u
import gc
import os 

def get_croppeddata(hdu, region, square=False):
    
    hdu_crop = hdu.copy()  # Create a copy of the input HDU object
    del hdu  # Delete the original HDU object to free up memory
    _ = gc.collect()  # Perform garbage collection

    wcs = WCS(hdu_crop)  # Create a WCS object from the HDU header

    try: 
        if square:
            radius = max([region['width'], region['height']]) # Calculate the radius for the square cutout 
            cutout = Cutout2D(hdu_crop.data, region['position'], radius, wcs=wcs)  # Create a square cutout
        else:
            cutout = Cutout2D(hdu_crop.data, region['position'], [region['width'], region['height']], wcs=wcs)  # Create a rectangular cutout
    except: 
        print('[INFO] NO cross-over')
        return(None)

    hdu_crop.data = cutout.data  # Update the data in the cropped HDU object
    hdu_crop.header.update(cutout.wcs.to_header())  # Update the header of the cropped HDU object with the cutout's WCS information

    del cutout  # Delete the cutout to free up memory
    _ = gc.collect()  # Perform garbage collection

    return hdu_crop  # Return the cropped HDU object

def get_region(l, b, w, h, frame='galactic'): 
    region = {'position': SkyCoord(l=l *u.deg, b=b *u.deg, frame=frame), 
          'width': w *u.deg,
          'height': h *u.deg}
    return(region)

# Get a list of all .fits files in the 'data' directory
files = glob('/Users/abarnes/Dropbox/work/Smallprojects/aces/data/misc/raw/regrid/*.fits')
files.sort()

# os.system('rm -rf regrid_crop/*')

l = 0.8065474
b = -0.1999737
width = 0.2
height = 0.2
region = get_region(l, b, width, height)

# Loop over each file. We're using tqdm to create a progress bar. 
# tqdm automatically determines the total number of iterations from the length of the 'files' list.
for file in tqdm(files, desc="Processing files", unit="file"):
    
    print('[INFO] infile - %s' % file)

    # Open the FITS file
    hdu = fits.open(file)[0]

    hdu_crop = get_croppeddata(hdu, region)

    if hdu_crop is not None: 

        output_file = file.replace('regrid', 'regrid_crop')
        hdu_crop.writeto(output_file, overwrite=True)