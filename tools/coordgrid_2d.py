import numpy as np
import astropy.wcs
from astropy.io import fits
from astropy import coordinates, units as au
from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft

def get_coordinate_grid_gal(hdu, frame='fk5'):

    """Returns world coordinate grid of an image
        [Image has to be in J2000 at the moment]
        Gives both RA/DEC and GALACTIC
        Input:
            hdu = Input fits hdu in J2000 coordinates
            frame = Input reference frame; default is fk5
        Output:
            (coordinate_grid_l, coordinate_grid_b) = coordinates of each array index in GALACTIC
            (coordinate_grid_ra, coordinate_grid_dec) = coordinates of each array index in RA/DEC
        """

    wcs = astropy.wcs.WCS(hdu)
    shape = hdu.shape

    index_grid_x = np.empty(shape)
    index_grid_y = np.empty(shape)
    coordinate_grid_ra = np.empty(shape)
    coordinate_grid_dec = np.empty(shape)
    coordinate_grid_l = np.empty(shape)
    coordinate_grid_b = np.empty(shape)

    for row in range(shape[1]):
        for col in range(shape[0]):

            index_grid_x[col, row] = row
            index_grid_y[col, row] = col

    coordinate_grid_radec = wcs.all_pix2world(index_grid_x, index_grid_y, 0)
    coordinate_grid_ra, coordinate_grid_dec = coordinate_grid_radec

    lb = coordinates.SkyCoord(coordinate_grid_ra,
                              coordinate_grid_dec,
                              frame=frame,
                              unit=(au.deg, au.deg)).galactic

    coordinate_grid_l = lb.l.deg
    coordinate_grid_b = lb.b.deg

    return (coordinate_grid_l, coordinate_grid_b), (coordinate_grid_ra, coordinate_grid_dec)
