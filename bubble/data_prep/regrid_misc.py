from glob import glob
from astropy.io import fits
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs
import numpy as np
from tqdm.auto import tqdm
import os 

# Define the conversion function
def convert_to_float32(hdu):
    hdu.data = hdu.data.astype('float32')
    return(hdu)

def reproject_to_projection(fits_file, output_suffix='regrid', frame='galactic', projection='TAN'):
    """
    Function to reproject a FITS file into a specified coordinate frame and projection.

    Parameters:
    fits_file (str): Input FITS file path.
    output_suffix (str, optional): The suffix for the output file name. Default is 'regrid'.
    frame (str, optional): Coordinate frame for reprojection. Default is 'galactic'.
    projection (str, optional): Projection system for reprojection. Default is 'TAN' (tangent plane).

    """

    print('[INFO] infile - %s' % fits_file)

    # Open the FITS file
    hdu = fits.open(fits_file)[0]
    if hdu.data is None:
        hdu = fits.open(fits_file)[1]
    hdu = fits.PrimaryHDU(hdu.data, hdu.header)
    hdu = convert_to_float32(hdu) 

    # The data array may have degenerate (i.e., length-1) axes, particularly if the original data was 3D or 4D. 
    # Squeeze function is used to remove these before reprojection.
    hdu.data = np.squeeze(hdu.data)

    # We are going to reproject the 2D image data, so the header entries for the third and fourth axes are not needed. 
    # These entries are removed from the header.
    # for key in list(hdu.header.keys()):
    #     if '3' in key or '4' in key:
    #         del hdu.header[key]

    del hdu.header['*3*']
    del hdu.header['*4*']
    del hdu.header['*HISTORY*']
    del hdu.header['*WAVELENG*']

    # Find the optimal World Coordinate System (WCS) for reprojection.
    # The frame parameter specifies the desired coordinate frame (e.g., 'icrs' for equatorial coordinates, 'galactic' for Galactic coordinates).
    # The projection parameter specifies the desired map projection (e.g., 'TAN' for tangent plane, 'CAR' for Cartesian).
    wcs_out, shape_out = find_optimal_celestial_wcs(hdu, frame=frame, projection=projection)

    # Reproject the data using the reproject_interp function. 
    # This function interpolates the data from the input map to the pixel grid of the output map.
    # The output map is defined by the wcs_out and shape_out parameters.
    # The function returns the reprojected data array.
    try: 
        data_r = reproject_interp(hdu, wcs_out.to_header(), shape_out=shape_out, return_footprint=False, parallel=True)
    except: 
        data_r = reproject_interp(hdu, wcs_out.to_header(), shape_out=shape_out, return_footprint=False)

    # Create the output filename by replacing 'raw' in the input filename with the specified output_suffix.
    output_file = fits_file.replace('raw', output_suffix)
    print('[INFO] outfile - %s' % output_file)

    # Create a new FITS HDU (Header Data Unit) with the reprojected data and the new WCS.
    hdu_r = fits.PrimaryHDU(data_r, wcs_out.to_header())

    # Write the new FITS file.
    # If a file with the same name already exists, it will be overwritten.
    hdu_r.writeto(output_file, overwrite=True)


# Get a list of all .fits files in the 'data' directory
files = glob('/Users/abarnes/Dropbox/work/Smallprojects/aces/data/misc/raw/*.fits')
files.sort()

# os.system('rm -rf regrid/*')

# Loop over each file. We're using tqdm to create a progress bar. 
# tqdm automatically determines the total number of iterations from the length of the 'files' list.
for file in tqdm(files, desc="Processing files", unit="file"):
    # Call the reproject function on each file. This function is assumed to handle
    # all the reprojecting tasks like opening the file, manipulating the data,
    # and saving the new reprojected file.
    reproject_to_projection(file)
