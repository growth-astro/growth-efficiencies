# code to get standards isolated and in a magnitude range
# Tomas Ahumada June 14 2019

import os, sys
from matplotlib import pyplot as plt
import requests

import numpy as np
import scipy.spatial as spatial
from penquins import Kowalski
from photometry_utils import panstarrs_query

from astropy.table import Table
from astropy.io import fits
from scipy.ndimage.filters import gaussian_filter
import astropy.wcs

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
    from urllib.request import urlretrieve
    from urllib import request
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen
    from urllib import urlretrieve


def get_ps1_filtered_table(k,name,rad_query=3,m_min=17,m_max=19.5,d=17,plot=False,control = False):
    '''
    Input
    k: kowalski pointer to query
    name: ztf ID
    rad : radius in arcminutes, bare in mind fov KPED 4.4' , LOT 12' , GIT 60' , LT 10' , P60 13'
    m_min,m_max: magnitudes to filter
    d : distance evaluated in order to determine isolation, in arcsec
    plot: set True to scatter plot ps1 stars
    control: set True to print sources
    Output:
    filtered table in mag (m_min,m_max) and distace 
    '''
    
    q = {"query_type": "general_search", 
         "query": "db['ZTF_alerts'].find({'objectId': {'$eq': '"+name+"'}})" 
         }
    r = k.query(query=q,timeout=10)
    if len(r['result_data']['query_result']) == 0:
        print ('query with empty results')
        return
    else:
        candidate = r['result_data']['query_result'][0]
        ra,dec = candidate['candidate']['ra'],candidate['candidate']['dec']

    radius_deg = 1/60*rad_query # KPED fov 4.4' , LOT 12' , GIT 60' , LT 10' , P60 13'
    ps1_table = panstarrs_query(ra,dec,radius_deg,minmag = 10,maxmag = 21.5)

    coords = np.asarray([ps1_table['ra'],ps1_table['dec']])
    coords = coords.T
    g_mag = ps1_table['gmag']
    r_mag = ps1_table['rmag']
    i_mag = ps1_table['imag']
    z_mag = ps1_table['zmag']
    
    # checking for distances
    point_tree = spatial.cKDTree(coords)
    ind_isolated = []
    for i in range(len(coords)):
        if len(point_tree.query_ball_point(coords[i], 1/60/60*d))==1:
            ind_isolated.append(i)
    ind_isolated = np.asarray(ind_isolated)
    coords_iso=coords[ind_isolated]
    g_iso = g_mag[ind_isolated]
    r_iso = r_mag[ind_isolated]
    i_iso = i_mag[ind_isolated]

    print('# of stars isolated (%d arcsec): '%d,len(ind_isolated))

    if plot== True:
        plt.scatter(coords.T[0],coords.T[1])
        plt.scatter(coords[ind_isolated].T[0],coords[ind_isolated].T[1],c='k')
        plt.show()
    if control == True:
        print('# ra dec gmag rmag imag')
        for i in range(len(ind_isolated)):
            print(coords_iso[i][0],coords_iso[i][1],np.round(g_iso[i],2),np.round(r_iso[i],2),np.round(i_iso[i],2))

    # checking for magnitude / filtering

    cond_mag = (g_iso < m_max) * (r_iso < m_max) * (g_iso > m_min) * (r_iso > m_min) 
    coords_filt=coords_iso[cond_mag]
    g_filt = g_iso[cond_mag]
    r_filt = r_iso[cond_mag]
    i_filt = i_iso[cond_mag]
    
    if plot== True:

        plt.scatter(coords.T[0],coords.T[1])
        plt.scatter(coords_filt.T[0],coords_filt.T[1],c='k')
        plt.show()

    # filtered and isolated stars 

    print('# of stars bright and isolated (25 arcsec , 17 < mag < 19): ',len(coords_filt))
    f = open(name+'.txt','w')
    print('# ra dec gmag rmag imag')
    f.write('# ra dec gmag rmag imag\n')
    for i in range(len(coords_filt)):
        print(coords_filt[i][0],coords_filt[i][1],np.round(g_filt[i],2),np.round(r_filt[i],2),np.round(i_filt[i],2))
        f.write('%.7f %.7f %.2f %.2f %.2f \n'%(coords_filt[i][0],coords_filt[i][1],g_filt[i],r_filt[i],i_filt[i]))
    f.close()

    
    return ra,dec,coords_filt

def get_fits_image(ra, dec, rad, debug=False):
    '''
    Connects to the PS1 or SkyMapper image service to retrieve the fits file to be used as a bse for the finder chart.
    '''
    #If dec> -30, we have Pan-STARRS
    if dec > -30:
        # Construct URL to download Pan-STARRS image cutout, and save to tmp.fits
    
        # First find the index of images and retrieve the file of the image that we want to use.
        image_index_url = 'http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra={0}&dec={1}&filters=r'.format(ra, dec)
        urlretrieve(image_index_url, '/tmp/ps1_image_index.txt')
        ix = Table.read('/tmp/ps1_image_index.txt', format="ascii")
        f = ix['filename'].data[0]
        
        image_url = "http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?red={0}&format=fits&size={1}&ra={2}&dec={3}".format(f, int(np.round(rad*3600*4, 0)), ra, dec)
        if (debug):
            print ("URL:", image_url)
            print ("Downloading PS1 r-band image...")
            
        #Store the object to a fits file.
        urlretrieve(image_url, '/tmp/tmp.fits')
        
            
    #Otherwise, we have SkyMapper   
    else:    
        url="http://skymappersiap.asvo.nci.org.au/dr1_cutout/query?POS=%.6f,%.6f&SIZE=%.3f&FORMAT=image/fits&INTERSECT=center&RESPONSEFORMAT=CSV"%(ra, dec, rad)
        page = urlopen(url)
        content = page.read()
        f = open("/tmp/skymapper_image_index.csv", "wb")
        f.write(content)
        f.close()
        
        ix = Table.read('/tmp/skymapper_image_index.csv', format="ascii.csv")

        mask = ((ix['band']=='r')|(ix['band']=='g'))
        
        ix = ix[mask]
        ix.sort(keys='exptime')
    
        image_url = ix['get_image'][-1]
        urlretrieve(image_url, '/tmp/tmp.fits')



    #Finally, once we have Pan-STARRS or SkyMapper images, we try to open them.
    #If there has been any problem with that, we will just go to the DSS image service.
    try:
        image = fits.open("/tmp/tmp.fits")
        #If everything went well, it shall be a fits image and opening it shall cause no issue.
        return '/tmp/tmp.fits'

        #If there was an error with the fits, we shall go for the DSS image
    except IOError:
        #One of the services may fail, so we need to account for that and provide a backup DSS image service.
        try:
            image_url = 'http://archive.eso.org/dss/dss/image?ra=%.5f&dec=%.5f&x=%.2f&y=%.2f&Sky-Survey=DSS1&mime-type=download-fits' %                 ((ra), (dec), (rad*60), (rad*60))
            if debug: print ("Downloading DSS image...")
            urlretrieve(image_url, '/tmp/tmp.fits')
        except:
            image_url = 'http://archive.stsci.edu/cgi-bin/dss_search?ra=%.6f&dec=%.6f&generation=DSS2r&equinox=J2000&height=%.4f&width=%.4f&format=FITS' %                 (ra, dec, rad*60, rad*60)
            urlretrieve(image_url, '/tmp/tmp.fits')
            
        #We try one more time to open it. If not successful, we return None as the image filename.
        try:
            fits.open("/tmp/tmp.fits")
        except IOError:
            print ("Your fits image could not be retrieved.")
            return None
    

def plot_references(ra, dec,rad_query,coords_filt,image_file = None,debug = True,plot=False):
    radius_deg = 1/60*2*rad_query
    if image_file is None:
            image_file = get_fits_image(ra, dec, radius_deg*1.1, debug=debug)    
            image = fits.open(image_file)

    # Get pixel coordinates of SN, reference stars in DSS image
    wcs = astropy.wcs.WCS(image[0].header)
    target_pix = wcs.wcs_world2pix([(np.array([ra,dec], np.float_))], 1)

    # Plot finder chart

    #Adjust some of the counts to make easier the plotting.
    image[0].data[image[0].data>30000] = 30000
    image[0].data[np.isnan(image[0].data)] = 0

    plt.figure(figsize=(10,10))
    plt.set_cmap('gray_r')
    smoothedimage = gaussian_filter(image[0].data, 1.3)
    plt.imshow(smoothedimage, origin='lower',vmin=np.percentile(image[0].data.flatten(), 10), vmax=np.percentile(image[0].data.flatten(), 99.0))

    # Mark target
    plt.plot([target_pix[0,0]+15,(target_pix[0,0]+10)],[target_pix[0,1],(target_pix[0,1])], 'g-', lw=2)
    plt.plot([target_pix[0,0],(target_pix[0,0])],[target_pix[0,1]+10,(target_pix[0,1])+15], 'g-', lw=2)
    plt.annotate(name, xy=(target_pix[0,0], target_pix[0,1]),  xycoords='data',xytext=(22,-3), textcoords='offset points')

    # # Mark and label reference stars
    ref_pix = wcs.wcs_world2pix(coords_filt, 1)
    for i,ref1_pix in enumerate(ref_pix):
    # if (len(catalog)>0):
        plt.plot([ref1_pix[0]+20,(ref1_pix[0]+15)],[ref1_pix[1],(ref1_pix[1])], 'b-', lw=2)
        plt.plot([ref1_pix[0],(ref1_pix[0])],[ref1_pix[1]+15,(ref1_pix[1])+20], 'b-', lw=2)
        plt.annotate(str(i), xy=(ref1_pix[0], ref1_pix[1]),  xycoords='data',xytext=(22,-3), textcoords='offset points', color="r")
    plt.xlim(0,len(image[0].data))
    plt.ylim(0,len(image[0].data))
    plt.hlines(50,10,10+60/0.258)
    plt.annotate('1 arcmin', xy=(10+10/0.258, 50),  xycoords='data',xytext=(0,3), textcoords='offset points', color="k")#,fontsize=20)
    plt.savefig(name+'.png',dpi=300)
    if plot:
        plt.show()

plot =  False
control = False
image_file = None
debug = True
rad_query=3

username = 'kped'
password = 'queryitEDdy!'

k = Kowalski(username=username, password=password, verbose=False)

names = [
'ZTF19aarykkb',  
'ZTF19aarzaod', 
'ZTF19aasckwd',
'ZTF19aasckkq',
'ZTF19aasbphu',
'ZTF19aaryxjf',
'ZTF19aarxxwb',
'ZTF19aasdajo',
'ZTF19aasbamy',
'ZTF19aarycuy',
'ZTF19aasbaui',
'ZTF19aasejil',
'ZTF19aascxux',
'ZTF19aashlts',
'ZTF19aasfogv',
#'PGIR19bn',
    ]

for name in names:
    ra,dec,coords_filt = get_ps1_filtered_table(k,name,rad_query=rad_query,m_min=17,m_max=20.5,d=15,plot=False,control = False)
    plot_references(ra, dec,rad_query,coords_filt,image_file = None,debug = True,plot=False)


