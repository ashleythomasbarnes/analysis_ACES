import os
from glob import glob

def regrid_to_galactic(intput_file):

    temp_file = intput_file.replace('.fits', '.img')
    output_file = intput_file.replace('.fits', '_regrid_casa.fits')
    
    importfits(fitsimage=intput_file, imagename=temp_file, overwrite=True)
    
    imregrid(imagename=temp_file, template='Galactic', output=temp_file+'.regrid', overwrite=True)
    
    exportfits(imagename=temp_file+'.regrid', fitsimage=output_file, overwrite=True, velocity=True)
    
    # Step 5: Cleanup Temporary Images
    os.system('rm -rf ' + temp_file)
    os.system('rm -rf ' + temp_file+'.regrid')

    os.system('rm -rf *.log')
    os.system('rm -rf *.last')

# Get a list of all .fits files in the 'data' directory
file = './data_processed/hdu.fits'
regrid_to_galactic(file)