# to create bin density data points for plotting over scatter plot

import numpy as np
from scipy import ndimage

def get_bins_(data, nbins=100, logbins=False):

    """Return bins for pixel density"""

    min=np.nanmin(data)
    max=np.nanmax(data)

    if logbins:
        min = np.nanmin(data[data>0])
        bins = np.logspace(np.log10(min), np.log10(max), nbins+1)
    else:
        bins = np.linspace(min, max, nbins+1)

    return(bins)

def get_bindensity_2d(data1, data2, xbins='', ybins='', nbins=100, logbins=False, sigma=1):

    """Return x, y, and contours for pixel density plot
        Input:
            data1 = input x-data
            data2 = input y-data
            xbins = input x-bins, if not given determined from min/max of x-data
            ybins = input y-bins, if not given determined from min/max of y-data
            nbins = input number of bins, but ignored if bins given: default: 100
            logbins = log-spaced bins, but ignored if bins given: default: False
            sigma = smoothing parameter for contours: default = 1 
        Output:
            binsx_c = x-axis bin points
            binsy_c = y-axis bin points
            counts = contour data points (normalised)
        Notes:
            Plot using e.g.:
                binsx_contour, binsy_contour, contour = get_bindensity_2d(data1, data2, logbins=True, nbins=70)
                levels = np.array([0.25, 0.5, 0.75, 1]) * np.nanmax(contour)
                ax1.contour(binsx_contour, binsy_contour, contour, levels=levels, linewidths=0.75, alpha=1, colors='red')
    """

    data1 = data1.flatten()
    data2 = data2.flatten()

    id_nonnan = np.where(~np.isnan(data1))
    data1 = data1[id_nonnan]
    data2 = data2[id_nonnan]

    id_nonnan = np.where(~np.isnan(data2))
    data1 = data1[id_nonnan]
    data2 = data2[id_nonnan]

    binsx = get_bins_(data1, nbins=nbins, logbins=logbins)
    binsy = get_bins_(data2, nbins=nbins, logbins=logbins)

    counts, _, _ = np.histogram2d(data1, data2, bins=(binsx, binsy))
    counts = np.ma.masked_where(counts < 0.1, counts)

    binsx_c = np.empty(len(binsx) - 1) * np.nan

    for i in range(len(binsx_c)):
        binsx_c[i] = np.nanmean(binsx[i:i+1])

    binsy_c = np.empty(len(binsy) - 1) * np.nan
    for i in range(len(binsy_c)):
        binsy_c[i] = np.nanmean(binsy[i:i+1])

    counts = ndimage.gaussian_filter1d(counts, sigma, 1)

    return binsx_c, binsy_c, counts.T
