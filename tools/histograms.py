import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV

def get_hist(data, bins='', nbins=50, logbins=False, norm=True, cum=False):

    """Get histogram

    Parameters
    ----------
    data : np.array
        input data
    bins : list
        input bin edges for histogram calculaiton; default=''
    nbins : int
        number of bins to determine if bins is not given; defult=50
    logbins : bool
        logarithmically spaced bins if bins is not given
    norm : bool
        normalise such that max is equal to unity; default=True
    cum : bool
        cumulative distorbution; otherwise probability distorbution
    Returns
    -------
    bins : list
        bin edges for histogram calculaiton
    bin_cent : np.array
        bin centres, for easy plotting in matplotlib
    hist : np.array
        histogram data for each bin centre
    """

    data = data.flatten()

    if bins=='':
        vmin=np.nanmin(data)
        vmax=np.nanmax(data)

        bmin = vmin - (np.absolute(vmin)*1)
        bmax = vmax + (np.absolute(vmax)*0.3)

        if logbins:
            min = np.nanmin(data[data>0])
            bins = np.logspace(np.log10(bmin), np.log10(bmax), nbins+1)
        else:
            bins = np.linspace(bmin, bmax, nbins+1)
    else:
        nbins = len(bins)-1

    bins_cent = np.empty([nbins])

    for i in range(nbins):
        bins_cent[i] = np.nanmean([bins[i], bins[i+1]])

    hist = np.histogram(data.flatten(), bins)[0]

    if cum:
        hist = np.cumsum(hist)
    if norm:
        hist = hist/np.nanmax(hist)

    return(bins, bins_cent, hist)


def get_histvalwt(data, bins='', nbins=50, logbins=False, norm=True, cum=False, cum_err=False):


    """Get value weighted histogram of dataset, using bins from another dataset

    Parameters
    ----------
    data : np.array
        input data
    bins : list
        input bin edges for histogram calculaiton; default=''
    nbins : int
        number of bins to determine if bins is not given; defult=50
    logbins : bool
        logarithmically spaced bins if bins is not given
    norm : bool
        normalise such that abs(max) is equal to unity; default=True
    cum : bool
        cumulative distorbution; otherwise probability distorbution
    Returns
    -------
    bins : list
        bin edges for histogram calculaiton
    bin_cent : np.array
        bin centres, for easy plotting in matplotlib
    hist : np.array
        histogram value weighted (i.e. summed) data for each bin centre
        - instead of returning the numbe of points within each bin
            returns the sum of the values within the bins
    """

    data = data.flatten()

    if bins=='':
        vmin=np.nanmin(data)
        vmax=np.nanmax(data)

        bmin = vmin - (np.absolute(vmin)*1)
        bmax = vmax + (np.absolute(vmax)*0.3)

        if logbins:
            min = np.nanmin(data[data>0])
            bins = np.logspace(np.log10(bmin), np.log10(bmax), nbins+1)
        else:
            bins = np.linspace(bmin, bmax, nbins+1)
    else:
        nbins = len(bins)-1

    bins_cent = np.empty([nbins])
    hist = np.empty(nbins)*np.nan

    for i in range(nbins):

        bins_cent[i] = np.nanmean([bins[i], bins[i+1]])
        ids = np.where((data > bins[i]) & (data < bins[i+1]))
        hist[i] = np.nansum(data[ids])

    if cum:
        hist = np.cumsum(hist)
    elif cum_err:
        hist = np.sqrt(np.cumsum(hist**2))
    if norm:
        # hist = hist / np.nanmax(np.absolute(hist))
        hist = hist / np.nanmax(hist)

    return(bins, bins_cent, hist)

def get_histvalwt_fromanother(datax, datay, bins='', nbins=50, logbins=False, norm=True, cum=False, cum_err=False):


    """Get value weighted histogram of dataset, using bins from another dataset

    Parameters
    ----------
    datax : np.array
        input x-axis data (used for binning)
    datay : np.array
        input y-axis data (data to bin)
    bins : list
        input bin edges for histogram calculaiton; default=''
    nbins : int
        number of bins to determine if bins is not given; defult=50
    logbins : bool
        logarithmically spaced bins if bins is not given
    norm : bool
        normalise such that abs(max) is equal to unity; default=True
    cum : bool
        cumulative distorbution; otherwise probability distorbution
    cum_err : cool
        determines error on cumulative distobution as sqrt(sum(data^2))
    Returns
    -------
    bins : list
        bin edges for histogram calculaiton
    bin_cent : np.array
        bin centres, for easy plotting in matplotlib
    hist : np.array
        histogram value weighted (i.e. summed) data for each bin centre
        - instead of returning the numbe of points within each bin
            returns the sum of the values within the bins
    """

    datax = datax.flatten()
    datay = datay.flatten()

    if bins=='':
        min=np.nanmin(datax)
        max=np.nanmax(datax)
        if logbins:
            min = np.nanmin(datax[datax>0])
            bins = np.logspace(np.log10(min), np.log10(max), nbins+1)
        else:
            bins = np.linspace(min, max, nbins+1)
    else:
        nbins = len(bins)-1

    bins_cent = np.empty([nbins])
    hist = np.empty(nbins)*np.nan

    for i in range(nbins):

        bins_cent[i] = np.nanmean([bins[i], bins[i+1]])
        ids = np.where((datax > bins[i]) & (datax < bins[i+1]))
        hist[i] = np.nansum(datay[ids])

    if cum:
        hist = np.cumsum(hist)
    elif cum_err:
        hist = np.sqrt(np.cumsum(hist**2))
    if norm:
        # hist = hist / np.nanmax(np.absolute(hist))
        hist = hist / np.nanmax(hist)

    return(bins, bins_cent, hist)

# import numpy as np
# from sklearn.neighbors import KernelDensity
# from sklearn.model_selection import GridSearchCV

def get_histkde(data, bins='', nbins=50, norm=True, cum=False):

    """Get Kernel Density Estimation histogram
    See: https://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/
    See: https://seaborn.pydata.org/tutorial/distributions.html
    Note: Do not try in log space

    Parameters
    ----------
    data : np.array
        input data
    bins : list
        input bin edges for histogram calculaiton; default=''
    nbins : int
        number of bins to determine if bins is not given; defult=50
    logbins : bool
        logarithmically spaced bins if bins is not given
    norm : bool
        normalise such that max is equal to unity; default=True
    cum : bool
        cumulative distorbution; otherwise probability distorbution

    Returns
    -------
    bins : list
        bin edges for histogram calculaiton
    bin_cent : np.array
        bin centres, for easy plotting in matplotlib
    hist : np.array
        histogram data for each bin centre
    """

    data = data.flatten()
    nsample = len(data)

    if bins=='':
        vmin=np.nanmin(data)
        vmax=np.nanmax(data)

        bmin = vmin - (np.absolute(vmin)*1)
        bmax = vmax + (np.absolute(vmax)*0.3)

        bins = np.linspace(bmin, bmax, nbins+1)
    else:
        nbins = len(bins)-1

    bins_cent = np.empty([nbins])

    for i in range(nbins):
        bins_cent[i] = np.nanmean([bins[i], bins[i+1]])

    data_ = data[:, np.newaxis]
    bins_cent_ = bins_cent[:, np.newaxis]

    #Bandwidth Cross-Validation in Scikit-Learn
    #Fuck knows, see: https://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/
    if nsample <20:
        grid = GridSearchCV(KernelDensity(), {'bandwidth': bins_cent}, cv=nsample) # 20-fold cross-validation
    else:
        grid = GridSearchCV(KernelDensity(), {'bandwidth': bins_cent}, cv=20) # 20-fold cross-validation
    grid.fit(data_)
    bw = grid.best_params_['bandwidth']

    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_)

    log_density_values = kde.score_samples(bins_cent_)
    hist = np.exp(log_density_values)

    if cum:
        hist = np.cumsum(hist)
    if norm:
        hist = hist/np.nanmax(hist)

    return(bins, bins_cent, hist)
