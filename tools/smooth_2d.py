import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft

def get_kernel(res_t, res_c, pix_c):

    """Get smoothing kernel"""

    fwhm_factor = np.sqrt(8*np.log(2))
    gaussian_width = ((res_t**2 - res_c**2)**0.5 / pix_c / fwhm_factor)

    print("[INFO] Current pixel size of %0.1f arcsec" %(pix_c))
    print("[INFO] Current resolution of %0.1f arcsec" %(res_c))
    print("[INFO] Target resolution of %0.1f arcsec" %(res_t))
    print("[INFO] Smoothing with gaussian of %0.1f arcsec" % (gaussian_width * pix_c * fwhm_factor))

    kernel = Gaussian2DKernel(gaussian_width)

    return kernel

def get_spatialsmooth(hdu, res_t, res_c=''):

    """Spatially smooth data
        Input:
            hdu = Input hdu from fits file
            res_t = target resolution (in arcsec)
            res_c = current resolution (in arcsec) if not in header
        Output:
            hdu_out = Output hdu containing the smoothed data and updated header
    """

    data = hdu.data
    header = hdu.header

    if res_c=='':
        beam_maj = header['BMAJ']
        beam_min = header['BMIN']
        res_c = np.sqrt(beam_maj*beam_min)

    pix_ = header['CDELT1']
    pix_c = np.absolute(pix_) *3600. #convert to arcsec

    kernel = get_kernel(res_t, res_c, pix_c)

    data_sm = convolve_fft(data, kernel)

    header['BMAJ'] = res_t /3600. #convert to deg
    header['BMIN'] = res_t /3600. #convert to deg
    header['BPA'] = 0.0

    hdu_out = fits.PrimaryHDU(data_sm, header)

    return hdu_out
