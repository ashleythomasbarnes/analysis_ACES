import numpy as np
import astropy.units as au
import astropy.stats as stats
from astropy.io import fits
from spectral_cube import SpectralCube
from astropy.utils.console import ProgressBar

from scipy.special import erf
from scipy.optimize import minimize
# from scipy import stats

def get_rms(cube, rms_velocities, outputfile=''):

    """Get RMS map from a given cube

    Parameters
    ----------
    cube : spectral cube object
        Input cube
    rms_velocities : list of velocities used to return the rms - unitless
        List of rms values - e.g [[20,30],[50,60]]
    outputfile : string (optional)
        Output filename

    Returns
    -------
    rms_map : fits.PrimaryHDU object
        map of rms values
    """

    rms_cubetot = ['']*len(rms_velocities)

    for i,rms_velocity in enumerate(rms_velocities):

        rms_velocity = rms_velocity*(au.km/au.s)

        rms_cube = cube.spectral_slab(rms_velocity[0], rms_velocity[1])
        rms_cubetot[i]=rms_cube.hdu.data

    rms_cubetot = np.vstack(rms_cubetot)
    rms_map = np.nanstd(rms_cubetot, axis = 0)

    header = cube.header
    del header['*3']
    del header['PV*']
    header['WCSAXES'] = 2
    rms_map = fits.PrimaryHDU(data=rms_map, header=header)

    if outputfile != '':
        rms_map.writeto(outputfile, overwrite = True)

    return rms_map


def erf0(x, erftarg):
    """
    Error function

    Parameters
    ----------
    x :
        x-sigma interval of which erf gives us the probability that noise is in range [-x,x]
    erftarg :
        targeted error function value

    Returns
    -------
    errorfunction :
    """
    return abs(erf(x)-erftarg)

def get_rmsrob(cube):

    """ Get ROBUST RMS map from a given cube

    To generate a map of the errors in a data cube by
   	taking the RMS along a pixel and then rejecting anything over
   	a threshold determined by the number of channels in a
   	spectrum, so that there is 25% chance of noise producing a spike
  	above the threshold. (original authours: I. Beslic & J. den Brok)

    Parameters
    ----------
    cube : spectral cube object
        Input cube

    Returns
    -------
    rms_map : fits.PrimaryHDU object
        map of rms values
    """
    # data, header = fits.getdata(file_to_cube, header=True)
    data = cube.hdu.data
    header = cube.hdu.header
    cube_dim = np.shape(data)

    # If cube does not have 3 dimensions, cannot do analysis, so return with error message
    if len(cube_dim)<3:
        raise IndexError("Holly Moly, Cube not 3-D! Try again.")

    if len(cube_dim)==4:
        if cube_dim[0]==1:
            data = data[0]
            cube_dim = np.shape(data)
        else:
            raise IndexError("Do not understand cube dimesnions, please provide cube in shape (z,x,y)")
    channels = cube_dim[0]
    x, y = cube_dim[1],cube_dim[2]
    # Initial guess (taken from errmap_rob.pro)
    x0 = 3
    # The targeted errf is defined like that by E. Rosollowski, we don't knwo why
    #print(channels)
    args = 1-(5-1)/channels

    sig_false = minimize(erf0, x0, args=(args)).x

    cube_map = np.zeros((x, y)) + stats.mad_std(data, axis=0)

    for i in range(x):

        for j in range(y):
            spec = data[:,i,j]

            # spec_neg: part of spectrum, where values are negative => Sure to only include noise
            # The following line is just a fancy way of ignoring the nans. It does the same as
            # spec_neg = spec[spec<0], but doesn't produce the runtime error warning

            spec_neg = spec[np.less(spec, 0., where=~np.isnan(spec))&  ~np.isnan(spec)]

            if len(spec_neg) < 10:
                continue

            concat = np.concatenate((spec_neg, -spec_neg))
            sigma = stats.mad_std(concat[~np.isnan(concat)])
            # simple code: cube_map[i,j] = np.nanstd(spec[spec<sig_false*sigma])
            spec_ = spec[np.less(spec, sig_false*sigma, where=~np.isnan(spec))]
            cube_map[i,j] = stats.mad_std(spec_[~np.isnan(spec_)])

    del header['*3']
    del header['PV1_*']
    header['WCSAXES'] = 2
    rms_map = fits.PrimaryHDU(data=cube_map, header=header)

    return rms_map
