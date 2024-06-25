#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import healpy as hp
from astropy.io import fits
import healpy as hp
from astropy.table import Table, Column
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.utils.data import get_readable_fileobj
from astropy import units as u
from ligo.skymap import io, plot, postprocess
from matplotlib import pyplot as plt
import glob

def GRB_unnoculted(path_GBM):
    """ Create a new skymap with the earth-oculted portion of the map masked
    
    new map in path.replace('.fit','_unocculted.fit')
    
    haversine and find_greedy_credible_levels from GBM_tools 
    https://fermi.gsfc.nasa.gov/ssc/data/analysis/rmfit/gbm_data_tools/gdt-docs/index.html
    """
    
    def haversine(lon1, lat1, lon2, lat2, deg=True):
        """Calculates the angular separation between two points using the
        haversine equation. If degrees are passed, degrees are returned. else
        the input/output is assumed to be radians.
        lon -> azimuth
        lat -> zenith

        Args:
            lon1 (float): lon/az of first point
            lat1 (float): lat/zen of first point
            lon2 (float): lon/az of second point
            lat2 (float): lat/zen of second point
            deg (bool, optional): True if input/output in degrees.

        Returns:
            float: Angular separation between points
        """
        if deg:
            lon1, lat1, lon2, lat2 = map(np.deg2rad, [lon1, lat1, lon2, lat2])
        d_lat = 0.5 * (lat2 - lat1)
        d_lon = 0.5 * (lon2 - lon1)

        a = np.sin(d_lat) ** 2 + (np.sin(d_lon) ** 2 * np.cos(lat1) * np.cos(lat2))
        alpha = 2. * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))

        if deg:
            alpha = np.rad2deg(alpha)

        return alpha

    def find_greedy_credible_levels(p):
        """Calculate the credible values of a probability array using a greedy
        algorithm.

        Args:
            p (np.array): The probability array

        Returns:    
             np.array: The credible values
        """
        p = np.asarray(p)
        pflat = p.ravel()
        i = np.argsort(pflat)[::-1]
        cs = np.cumsum(pflat[i])
        cls = np.empty_like(pflat)
        cls[i] = cs
        return cls.reshape(p.shape)

    # getting header information
    hdu = fits.open(path_GBM)
    geo_ra = hdu['HEALPIX'].header['GEO_RA']
    geo_dec = hdu['HEALPIX'].header['GEO_DEC']
    geo_location = [geo_ra,geo_dec]

    try:
        geo_radius = hdu['HEALPIX'].header['GEO_RAD']
    except:
        geo_radius = 67.5

    prob, sig= hp.read_map(path_GBM, field=(0, 1), memmap=False,
                                    verbose=False)
    # prob_mask = prob > 0.0
    npix = len(prob)
    nside = hp.npix2nside(npix)

    theta, phi = hp.pix2ang(nside, np.arange(npix))
    ra = np.rad2deg(phi)
    dec = np.rad2deg(np.pi / 2.0 - theta)
    ang = haversine(*geo_location, ra, dec)

    # masking earth
    geo_mask = (ang <= geo_radius)
    new_prob = np.copy(prob)
    new_prob[geo_mask] = 0.0
    # renormalize
    new_prob /= np.sum(new_prob)
    # have to redo the significance
    new_sig = 1.0 - find_greedy_credible_levels(new_prob)

    unocculted_path = path_GBM.replace('.fit','_unocculted.fit')
    # get arrays in proper order, and write the healpix data to disk
    prob_arr = hp.reorder(new_prob, r2n=True)
    sig_arr = hp.reorder(new_sig, r2n=True)
    columns = ['PROBABILITY', 'SIGNIFICANCE']
    hp.write_map(unocculted_path, (prob_arr, sig_arr), nest=True, coord='C',
                 overwrite=True, \
                 column_names=columns,
                 extra_header=hdu[1].header.cards)
    hdulist = fits.open(unocculted_path)
    hdulist[0].header.extend(hdu[0].header.cards)
    hdulist[1].name = 'HEALPIX'
    hdulist[1].header['TTYPE1'] = (
    'PROBABILITY', 'Differential probability per pixel')
    hdulist[1].header['TTYPE2'] = (
    'SIGNIFICANCE', 'Integrated probability')
    hdulist.writeto(unocculted_path, overwrite=True, checksum=True)

# make a map


def pretty_plot(map_path,ra,dec,cand_name='name'):
    
    # loading relevant info
    prob = hp.read_map(map_path)
    ipix = np.argmax(prob)
    npix = len(prob)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, ipix)
    
    # getting the center
    ra_max = np.rad2deg(phi)
    dec_max = np.rad2deg(0.5 * np.pi - theta)-10
    center = str(int(SkyCoord(ra_max,dec_max,unit='deg').ra.hms.h))+'h '+ str(int(SkyCoord(ra_max,dec_max,unit='deg').dec.dms.d))+'d'
    sky_area_per_pix = 4 * 180**2 / np.pi /npix
    area_95 = np.sum(np.cumsum(-np.sort(-prob))<=0.95) * sky_area_per_pix
    area_50 = np.sum(np.cumsum(-np.sort(-prob))<=0.5) * sky_area_per_pix
    print('AREA 50% ',area_50)
    print('AREA 95% ',area_95)
    
    #plotting
    ax = plt.axes(projection='astro globe', center=center)
    ax.imshow_hpx(prob, cmap='cylon')
    cl = postprocess.find_greedy_credible_levels(prob) * 100
    ax.contour_hpx(cl, levels=[50, 95], colors=['black', 'black'], linewidths=[0.5, 0.5])
    ax.grid()
    ax.scatter(ra,dec,marker='*',color='blue',transform=ax.get_transform('world'),zorder=10)
    ax.set_title(fits.open(map_path)[0].header['OBJECT']+' / '+str(fits.open(map_path)[0].header['TRIGTIME'])+' \n'+cand_name)
    plt.savefig(map_path.replace('.fit','_cand.png'),dpi=250)


# In[13]:



def in_out(map_path,ra_obj,dec_obj, top_fraction = 0.95 ):
    
"""
Calculate if a list of objects are inside 'top_fraction' localisation region of a given skymap
"""
    # Read skymap, calculate top pixels
#     top_fraction = 0.95 # limit skymap top 90% region
    skymap = hp.read_map(map_path)
    npix = len(skymap)
    nside = hp.npix2nside(npix)

    # Convert to astropy Table, easier for manipulation
    indices = np.arange(len(skymap))
    tm = Table(data=(indices, skymap), names=('id', 'prob'))
    tm.sort('prob')
    cs = np.cumsum(tm['prob'])
    cs.name='cumsum'
    tm.add_column(cs)

    top_pix = (tm['cumsum'] > 1 - top_fraction)
    tp = Column(data=top_pix, name="top")
    tm.add_column(tp)

    # Cast as a set for easier comparison below
    top_subset = set(tm['id'][tm['top']])

    inside = False
    pix = hp.ang2pix(nside,ra_obj,dec_obj, lonlat=True)
    if pix in top_subset:
        inside = True 
    return inside


# In[14]:


cands = {
    'ZTF20aamsouh':{'ra':158.86648,
                    'dec':31.63394,
                    'maps':['/Users/tahumada/Downloads/glg_healpix_all_bn200211813.fit'],
                   'window':' 2020-02-11 10:03:00 .. 2020-02-12 07:20:00'},
    'ZTF20aakqxsq':{'ra':67.13090,
                    'dec':49.44819,
                    'maps':['/Users/tahumada/Downloads/glg_healpix_all_bn200127758.fit',
                            '/Users/tahumada/Downloads/glg_healpix_all_bn200128153.fit'],
                   'window':'2020-01-28 06:41:00 .. 2020-01-29 07:15:00'}
        }



# In[46]:


for cand in cands:
    for path_GBM in cands[cand]['maps']:
        path_GBM_unnoc = path_GBM.replace('.fit','_unocculted.fit')
        GRB_unnoculted(path_GBM)
        ra,dec = cands[cand]['ra'],cands[cand]['dec']
        print(cand)
        print(path_GBM, '\n inside?',in_out(path_GBM_unnoc,ra,dec))
        pretty_plot(path_GBM_unnoc,ra,dec,cand_name=cand)


# In[ ]:




