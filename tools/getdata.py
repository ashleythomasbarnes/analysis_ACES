import numpy as np
import astropy.units as au
import astropy.stats as stats
from astropy.io import fits
from spectral_cube import SpectralCube
from astropy.utils.console import ProgressBar

def get_cube(inputfile, x_units=(au.km / au.s), y_units=au.K):

    """Get datacube

    Parameters
    ----------
    inputfile : str
        input filename
    x_units:
        units of x-axis
    y_units:
        units of y-axis

    Returns
    -------
    cube : spectral-cube object
        output cube in xunit and yunit
    """

    cube = SpectralCube.read(inputfile, format = 'fits', unit = y_units, beam_threshold = 0.05)
    cube = cube.with_spectral_unit(x_units, velocity_convention='radio')
    cube.beam_threshold = 0.5
    cube.allow_huge_operations=True

    if cube.unit == "Jy / beam":
        cube = cube.to(au.K)

    return cube
