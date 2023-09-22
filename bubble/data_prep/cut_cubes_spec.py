from spectral_cube import SpectralCube
from astropy import units as u
from astropy.io import fits

# Define the conversion function
def convert_to_float32(hdu):
    hdu.data = hdu.data.astype('float32')
    return(hdu)

def crop_cube_velocity_range(fits_file, rest_frequency, v_start, v_end):
    '''
    Crop a FITS cube to a given velocity range for a given rest frequency.

    Parameters:
    fits_file: str
        Path to the FITS file.
    rest_frequency: float
        Rest frequency to use for the conversion from frequency to velocity. Given in GHz.
    v_start, v_end: float
        Velocity range for the cropping. Given in km/s.

    Returns:
    cropped_cube: SpectralCube
        Cropped spectral cube.
    '''

    # Load the cube
    cube = SpectralCube.read(fits_file)
    cube.allow_huge_operations=True

    # Convert the cube to velocity space
    cube = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value=rest_frequency*u.GHz)

    # Define the velocity range for cropping
    vrange = [v_start*u.km/u.s, v_end*u.km/u.s]

    # Crop the cube
    cropped_cube = cube.spectral_slab(*vrange)
    
    # Further crop the cube to the minimal enclosing subcube
    cropped_cube = cropped_cube.minimal_subcube()

    cropped_cube = cropped_cube.to(u.K)

    hdu = cropped_cube.hdu

    hdu = fits.PrimaryHDU(hdu.data, hdu.header)
    hdu = convert_to_float32(hdu) 

    return hdu

fits_file = 'uid___A001_X15a0_X166.s38_0.Sgr_A_star_sci.spw33.cube.I.iter1.image.pbcor.fits'
rest_frequency = 99.02295 
v_start = -0
v_end = 200 

hdu = crop_cube_velocity_range(fits_file, rest_frequency, v_start, v_end)
hdu.writeto('12m_h40a_cube.fits', overwrite=True)