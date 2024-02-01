###  some simple modules for dealing with nans values in array

import os
import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import interpolate_replace_nans


def replacenans_2d(data, value = ''):

    """Replace nans within array with value - max value of array if no value given (good for saturated values)
    
        Input: 
            data = np.array of data values
            value = value to replace, if no value is given then maximum value within array is determined and used 
        Output: 
            data_out = np.array of data value with no nan values"""

    ids = np.isnan(data)
    data_out = np.copy(data)
    
    if value == '':
        data_out[ids] = np.nanmax(data_out)
        
    else:
        data_out[ids] = value
        
    return data_out


def replacenans_2d_interpolate(data):

    """Interpolate over nans within array using astropy interpolate_replace_nans function
    
        Input: 
            data = np.array of data values
        Output: 
            data_out = np.array of data value with no nan values"""

    
    # We smooth with a Gaussian kernel with x_stddev=1 (and y_stddev=1)
    # It is a 9x9 array
    kernel = Gaussian2DKernel(x_stddev=1)

    # create a "fixed" image with NaNs replaced by interpolated values
    data_out = interpolate_replace_nans(data, kernel)
        
    return data_out