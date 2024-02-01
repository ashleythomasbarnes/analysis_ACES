import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, Gaussian1DKernel, convolve, convolve_fft
from radio_beam import Beam
import astropy.units as au

def get_kernelspectral(res_t, res_c, pix_c):

    """Get smoothing kernel"""

    fwhm_factor = np.sqrt(8*np.log(2))
    gaussian_width = ((res_t**2 - res_c**2)**0.5 / pix_c / fwhm_factor)

    print("[INFO] Current pixel size of %0.1f kms-1" %(pix_c))
    print("[INFO] Current resolution of %0.1f kms-1" %(res_c))
    print("[INFO] Target resolution of %0.1f kms-1" %(res_t))
    print("[INFO] Smoothing with gaussian of %0.1f kms-1" % (gaussian_width * pix_c * fwhm_factor))

    kernel = Gaussian1DKernel(gaussian_width)

    return kernel


def get_spectralsmooth(cube, res_t, res_c=''):

    """Spectrally smooth data

    Parameters
    ----------
    cube : spectral cube object
    res_t : float
        target *resolusion* (in kms-1)
        *not* channel width - do get_spectralregrid
    res_c : float (optional)
        current resolution (in kms-1) if not in header

    Returns
    -------
    smcube = Output cube containing the spectrally smoothed data and updated header
    """

    if res_c == '':
        # Assuming beams are in kms-1
        res_c = np.absolute(cube.header['CDELT3'])
    pix_c = np.absolute(cube.header['CDELT3'])

    kernel = get_kernelspectral(res_t, res_c, pix_c)
    smcube = cube.spectral_smooth(kernel)

    return smcube


def get_spatialsmooth(cube, res_t):

    """Spatailly smooth data

    Parameters
    ----------
    cube : spectral cube object
    res_t : float
        target *resolusion* (in arcsec)
        *not* channel width - do get_spectralregrid

    Returns
    -------
    smcube = Output cube containing the spatailly smoothed data and updated header
    """


    res_c = cube.header['BMAJ'] *3600

    print("[INFO] Current resolution of %0.1f arcsec" %(res_c))
    print("[INFO] Target resolution of %0.1f arcsec" %(res_t))

    res_t_ = Beam(major = res_t*au.arcsec, minor = res_t*au.arcsec, pa = 0*au.deg)
    smcube = cube.convolve_to(res_t_)

    return smcube
