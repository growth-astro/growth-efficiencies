import warnings
import glob
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

from astroML.crossmatch import crossmatch_angular
from astropy.coordinates import SkyCoord
import os
from astropy.io import fits
from astropy import wcs

import urllib
import urllib.request

import astropy.units as u
import astropy.coordinates as coord

os.chdir('/app/data/KPED_new/') # change dir


from photometry_utils import sdss_query,panstarrs_query,get_ztf_cand,do_crossmatch


# folder to save output
# files with LCO data; used only *91.fits becuase they have the catalogs

# lco_save = 'science/LCO/lcogtdata-20190425-4/'
# files = glob.glob('./science/LCO/lcogtdata-20190425-4/*91.fits')

lco_save = 'science/LCO/lcogtdata-20190728-8/'
files = glob.glob('./science/LCO/lcogtdata-20190728-8/*91.fits')

# print(files)
# to log into the growth marshal
username='tahumada'
password='123momiaes'

# url for the project; change month or copy paste from growth marshal
# url_cand='http://skipper.caltech.edu:8080/cgi-bin/growth/list_sources_bare.cgi?programidx=1&startdate=2018-09-13+20%3A42%3A02&enddate=2018-09-17++20%3A42%3A02'

# getting all candidates
# name,ra_transient,dec_transient=get_ztf_cand(url_cand, username, password)
# ra_transient,dec_transient = 258.337812514, -9.96204739154

name='ZTF19abjethn'
ra_transient=326.395431 
dec_transient=+20.690590

# to save output
LCO = []
for image in files:

    #get the data image
    hdu = fits.open(image)
    for i in range(len(hdu)):
        try:
            if hdu[i].header['XTENSION'] == 'BINTABLE':
                bintable = i
                print('Catalog in extesion ', bintable)
        except:
            continue
    
    mag_aper=-2.5*np.log10(hdu[bintable].data['FLUXAPER5'])
    mag_err = hdu[bintable].data['FLUXERR5']/hdu[bintable].data['FLUXaper5']
    
    w = wcs.WCS(hdu[0].header)

    x = hdu[bintable].data.field(0) 
    y = hdu[bintable].data.field(1) 
    ra, dec = w.wcs_pix2world(x, y, 1)

    radius_deg = 0.08
    ra_header, dec_header = hdu[0].header['CRVAL1'],hdu[0].header['CRVAL2']
    
    print('Object: ',hdu[0].header['OBJECT'])
    print('RA:' ,ra_header,'DEC:', dec_header)
    
    # load PS1 data
    radius_deg = 0.08
    if filt == 'up':
        radius_deg = 0.18
        ps1_table = sdss_query(ra_header, dec_header, radius_deg)
    else:
        ps1_table = panstarrs_query(ra_header, dec_header, radius_deg,maxsources=300)
    
    
    filt = hdu[0].header['FILTER'][0]+'mag'
    
    # Getting Zero point
    ind_im,match_im = do_crossmatch(ra,dec,ps1_table['ra'],ps1_table['dec'])
    
    ZP_mag = []
    for ii in range(len(ind_im)):
        if ind_im[ii] < len(ps1_table) and np.abs(mag_err[ii])<0.025 :
            ZP_mag.append(-mag_aper[ii]+ps1_table[filt][ind_im[ii]])
    # summary
    print('number of standars used:', len(ZP_mag))
    print('ZP = ',np.round(np.mean(ZP_mag),3),'+-',np.round(np.std(ZP_mag),3))

    # Get Transient and save it in LCO list
    ind_tra,match_tra = do_crossmatch(ra_transient,dec_transient,ra,dec)
#     print(ind_tra,match_tra)
    if match_tra.any():
        print('LCO '+filt+' = '+str(np.round((mag_aper[ind_tra]+np.median(ZP_mag))[0],3))+' +- '+str(np.round(np.std(ZP_mag),2)))
        LCO.append([hdu[0].header['OBJECT'],1,filt,np.round((mag_aper[ind_tra]+np.median(ZP_mag))[0],3),np.round(np.std(ZP_mag),2)])
    else: 
        maglim = max(mag_aper[np.abs(mag_err)<0.5]+np.mean(ZP_mag))
        print('No transient detected ',filt,'lim = ',np.round(maglim,3))
        LCO.append([hdu[0].header['OBJECT'],0,filt,np.round(maglim,3),-99])
    print(' \n')
    

