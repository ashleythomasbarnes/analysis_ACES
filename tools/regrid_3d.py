import numpy as np
import astropy.units as au
import math

def get_spectralregrid(cube, res_t):

    """Spectrally smooth data

    Parameters
    ----------
    cube : spectral cube object
    res_t : float
        target *channel width* (in kms-1)
        *not* resolution - do get_spectralsmooth

    Returns
    -------
    hdu_out = Output hdu containing the smoothed data and updated header
    """

    spectral_axis = cube.spectral_axis
    spec_min = spectral_axis.min()
    spec_max = spectral_axis.max()

    new_spectral_axis = np.arange(spec_min.value, spec_max.value, res_t) *au.km/au.s
    interpcube = cube.spectral_interpolate(new_spectral_axis, suppress_smooth_warning=True)

    return interpcube

def get_spatialregrid(cube, header, res_t, size_t=''):

    """Spatially smooth data

    Parameters
    ----------
    cube : spectral cube object
        input cube
    header : hdu header
        header from input cube
        [bit annoying that this has to be given
         should be folded into the function]
    res_t : float
        target *pixel size* (in arcsec)
        *not* resolution - do get_spatialsmooth for this opertaion
    size_t : array (optional)
        target size of regridded map

    Returns
    -------
    interpcube : spectral cube object
        Output regridded spectral cube
    """

    cdelt1 = header['CDELT1']
    cdelt2 = header['CDELT2']
    naxis1 = header['NAXIS1']
    naxis2 = header['NAXIS2']
    crpix1 = header['CRPIX1']
    crpix2 = header['CRPIX2']

    # Assuming square pixels
    res_t = res_t *au.arcsec
    res_c = np.absolute(cdelt2 *3600 *au.arcsec)
    res_fact = (res_t/res_c).value

    if size_t=='':
        print ('[INFO] Regridding to same size as the input map')
        header['CDELT1'] = cdelt1 * res_fact
        header['CDELT2'] = cdelt2 * res_fact
        header['NAXIS1'] = int(naxis1 / res_fact)
        header['NAXIS2'] = int(naxis2 / res_fact)
        header['CRPIX1'] = int(crpix1 / res_fact)
        header['CRPIX2'] = int(crpix2 / res_fact)

    else:
        print ('[INFO] Regridding map to user given size')
        size_t = size_t * 3600. *au.arcsec
        print ('[INFO] Map size of %i x %i arcsec' % (size_t[0].value, size_t[1].value))

        header['CDELT1'] = cdelt1 *res_fact
        header['CDELT2'] = cdelt2 *res_fact

        xh = np.absolute(size_t[0]/(cdelt1 *3600 *au.arcsec))
        yh = np.absolute(size_t[1]/(cdelt2 *3600 *au.arcsec))
        header['NAXIS1'] = int(xh)
        header['NAXIS2'] = int(yh)
        header['CRPIX1'] = math.ceil(int(xh)/2.)
        header['CRPIX2'] = math.ceil(int(yh)/2.)

        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['RADESYS'] = 'ICRS'
        header['EQUINOX'] = 2000.0
        del header['ORIGIN']
        del header['DATE']

    cube_out = cube.reproject(header, order=u'bilinear')

    return cube_out
