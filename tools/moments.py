import numpy as np
import astropy.units as au
from astropy.io import fits
from spectral_cube import SpectralCube

def get_momentmaps(cube, rms, mom_velocity='', velomoms=False):

    """Get moment maps from a given (masked) cube

    Parameters
    ----------
    cube : spectral cube object
        Input cube
    rms : hdu
        rms map from cube
    mom_velocity : list (optional)
        list of velocities over which to deteremine moment maps - e.g. [40,50]

    Returns
    -------
    mom_out : dict
        dictionary of output moment maps - will be one or more of the following
        ['mom0', 'max', 'mom1', 'sigma', 'variance', 'fwhm']
    """

    try:
        velo_min = mom_velocity[0]*au.km / au.s
        velo_max = mom_velocity[1]*au.km / au.s
        cube_slab = cube.spectral_slab(velo_min, velo_max)

    except:
        print('[INFO] No velocity range.')
        cube_slab = cube

    mom0 = cube_slab.moment(order=0, axis=0).hdu
    maximum = cube_slab.max(axis=0).hdu

    if velomoms:
        mom1 = cube_slab.moment(order=1, axis=0).hdu
        sigma = cube_slab.linewidth_sigma().hdu
    else:
        mom1 = None
        sigma = None

    mom_out = dict.fromkeys(['mom0', 'max', 'mom1', 'sigma',
                            'mom0err', 'maxerr', 'mom1err', 'sigmaerr'])
    mom_out['mom0'] = mom0
    mom_out['max'] = maximum
    mom_out['mom1'] = mom1
    mom_out['sigma'] = sigma

    nchan = get_nchancube(cube)
    cdelt3 = np.abs(cube.hdu.header['CDELT3']) *au.km/au.s
    mom0err = get_mom0err(nchan, rms, cdelt3)
    mom0err.data = mom0err.data.value

    mom_out['mom0err'] = mom0err
    mom_out['maxerr'] = rms
    mom_out['mom1err'] = None
    mom_out['sigmaerr'] = None

    mom0s2n = fits.PrimaryHDU(mom0.data/mom0err.data, rms.header)
    maxs2n = fits.PrimaryHDU(maximum.data/rms.data, rms.header)

    mom_out['mom0s2n'] = mom0s2n
    mom_out['maxs2n'] = maxs2n
    mom_out['mom1s2n'] = None
    mom_out['sigmas2n'] = None

    return mom_out

def get_nchancube(cube):
    """Get map of number of channels along each pixel within a cube"""
    nchans = np.nansum(~np.isnan(cube)*1, axis=0)
    return(nchans)

def get_nchanmask(mask):
    """Get map of number of channels along each pixel within a mask"""
    arr1 = mask.include()
    arr2 = arr1 * 1
    nchans = np.sum(arr2, axis=0)
    return(nchans)

def get_nchanbool(data):
    """Gets number of channels witin each LOS"""
    nchans = np.nansum(data*1, axis=0)
    return(nchans)

def get_mom0err(nchan, rms_map, velo_res):
    """Get map of err on mom0
    sigma_W = rms * V_res * N_chan**0.5

    Parameters
    ----------
    nchan : (2d) np.array
        Input map of number of channels
    rms_map : fits hdu object
        Input map of rms (inc header for output hdu)
    velo_res : float
        Input velocity resolution (in kms-1 usally)

    Returns
    -------
    mom0err : fits hdu object
        Output map of err on mom0
    """
    cdelt = np.absolute(velo_res)
    nchan_sqrt = np.sqrt(nchan)

    mom0err = rms_map.data * nchan_sqrt * velo_res
    id_nan = np.where(mom0err==0)
    mom0err[id_nan] = np.nan

    header = rms_map.header
    header['BUNIT'] = 'K km s-1'
    mom0err_hdu = fits.PrimaryHDU(mom0err, header)
    return mom0err_hdu
