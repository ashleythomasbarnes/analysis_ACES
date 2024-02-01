import reproject
import numpy as np
from astropy.io import fits
from astropy import coordinates, units as au

def regrid_to_template_2d(hdu, hdu_template):

    """Regrids input fits image to given template file
        Input: 
            hdu = Input fits hdu to be regridded
            hdu_template = Input fits hdu containing template for data
        Output: 
            hdu_rg = fits hdu of regridded data 
    """

    header = hdu.header
    data = hdu.data
    
    header_template = hdu_template.header
    data_template = hdu_template.data
    
    data_rg, _ = reproject.reproject_interp(hdu, header_template)

    keys = ['CTYPE1', 
            'CRVAL1', 
            'CDELT1', 
            'CRPIX1', 
            'CROTA1', 
            'CROTA1', 
            'CRVAL2', 
            'CTYPE2', 
            'CDELT2', 
            'CRPIX2', 
            'CROTA2', 
            'EQUINOX', 
            'RADESYS', 
            'PV2_1', 
            'PV2_2', 
            'LONPOLE', 
            'LATPOLE']

    for key in keys:
        if key in header_template.keys():
            header[key] = header_template[key]

        elif key in header.keys():
            del header[key]
            
    hdu_rg = fits.PrimaryHDU(data_rg, header)

    return hdu_rg


##### Not completed
# def j2000_to_gal(hdu):
    
#     """Regrid data from j2000 to Galactic coordinates
#         Input: 
#             hdu = Input fits hdu to be regridded to J2000 coordinates
#         Output: 
#             hdu_gal = Output fits  hdu in coordinates of Galactic
#         """

#     data = hdu.data
#     header = hdu.header
#     header_gal = header.copy()
    
#     wcs = astropy.wcs.WCS(hdu)
    
#     shape = hdu_sm.shape
#     x_mid = shape[0]/2
#     y_mid = shape[1]/2
#     yx_max = np.max(shape)
    
#     radec_mid = wcs.array_index_to_world(x_mid, y_mid)
#     ra_mid = radec_mid.ra
#     dec_mid = radec_mid.dec

#     lb_mid = coordinates.SkyCoord(ra_mid, dec_mid, frame='icrs', unit=(au.deg, au.deg)).galactic
#     l_mid = lb_mid.l.deg
#     b_mid = lb_mid.b.deg
    
#     header_gal['CTYPE1'] = 'GLON-SIN'
#     header_gal['CTYPE2'] = 'GLAT-SIN'

#     header_gal['CRVAL1'] = l_mid
#     header_gal['CRVAL2'] = b_mid
    
#     header_gal['CRPIX1'] = y_mid
#     header_gal['CRPIX2'] = x_mid
    
# #     header_gal['NAXIS1'] = yx_max
# #     header_gal['NAXIS2'] = yx_max

#     data_gal, _ = reproject.reproject_interp(hdu, header_gal)
#     hdu_gal = fits.PrimaryHDU(data_gal, header_gal)

#     return hdu_gal