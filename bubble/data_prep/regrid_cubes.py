import os
from glob import glob

def regrid_to_galactic(intput_file):

    temp_file = intput_file.replace('.fits', '.img')
    output_file = intput_file.replace('.fits', '.galactic.fits')
    
    importfits(fitsimage=intput_file, imagename=temp_file, overwrite=True)
    
    imregrid(imagename=temp_file, template='Galactic', output=temp_file+'.regrid', overwrite=True)
    
    exportfits(imagename=temp_file+'.regrid', fitsimage=output_file, overwrite=True, velocity=True)
    
    # Step 5: Cleanup Temporary Images
    os.system('rm -rf ' + temp_file)
    os.system('rm -rf ' + temp_file+'.regrid')

    os.system('rm -rf *.log')
    os.system('rm -rf *.last')

# Get a list of all .fits files in the 'data' directory
# files = glob('/Users/abarnes/Dropbox/work/Smallprojects/aces/data/alma/12m7mtp_bubble/*.fits')
files = glob('/Users/abarnes/Dropbox/work/Smallprojects/aces/data/alma/12m_bubble/*.fits')
files.sort()

# Loop over each file. 
for file in files:

    if 'galactic' in file:
        continue

    regrid_to_galactic(file)