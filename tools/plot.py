from astropy.io import fits, ascii
from astropy.table import Table, join, vstack
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from astropy.wcs import WCS

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["ytick.right"] = 'True'
plt.rcParams["xtick.top"] = 'True'

def get_cmapscaling(data, base=1000):
    """Get log scaling"""
    array = ((data*base)+1)
    data_ = np.log(array)/np.log(base)
    return(data_)

def get_plot(data, hdu=None, scaling=99.9, base=1000):
    """Create simple plot of fits image in matplotlib"""
    try:
        data_ = get_cmapscaling(data.data)
    except:
        data_ = get_cmapscaling(data)

    vmin = np.nanpercentile(data_[data_!=0], 100-scaling)
    vmax = np.nanpercentile(data_[data_!=0], scaling)
#     print(vmin, vmax)
    try:
        wcs = WCS(data.header)
        fig = plt.figure()
        ax = fig.add_subplot(projection=wcs)
    except:
        fig = plt.figure()
        if hdu != None:
            wcs = WCS(hdu.header)
            ax = fig.add_subplot(projection=wcs)
        else:
            ax = fig.add_subplot()
    im = ax.imshow(data_, origin='lower', cmap='magma', vmin=vmin, vmax=vmax)
    return((fig, ax), (vmin, vmax))
