
import os
import optparse
import warnings
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import numpy as np
import healpy as hp

from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy.coordinates import Angle
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.time import Time
import ligo.skymap.plot
from ligo.skymap.io import fits
from ligo.skymap import postprocess
from ligo.skymap import distance

import matplotlib
#matplotlib.rc('text', usetex=True)
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
from matplotlib import patches
import matplotlib.pyplot as plt
from matplotlib import rcParams

from astropy.visualization.wcsaxes import SphericalCircle
from astropy.table import Table

import gwemopt.ztf_tiling

__version__ = 1.0

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-s", "--skymap", help="GW skymap.", default='../data/GW190425/LALInference.fits.gz')
    parser.add_option("-f","--fields", help="Observed fields.", default='../data/GW190425/ZTF_fields.dat')
    parser.add_option("-l","--transients", help="Transient list.", default='../data/GW190425/transients.dat')

    parser.add_option("-o", "--output", help="output file",default="../output/GW190425/LALInf.pdf")
    parser.add_option("-i", "--inputDir", help="input directory",default="../input/")

    parser.add_option("-t", "--telescope", help="Telescope.", default ="ZTF")
    parser.add_option("-g", "--gps", help="Event time GPS.", default=1240215503.011549, type=float)

    parser.add_option("--doTransients",  action="store_true", default=False)
    parser.add_option("--doMovie",  action="store_true", default=False)

    parser.add_option("--doToOOnly",  action="store_true", default=False)
    parser.add_option("--doChipGaps",  action="store_true", default=False)

    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely. (Default: False)")

    opts, args = parser.parse_args()

    # show parameters
    if opts.verbose:
        print >> sys.stderr, ""
        print >> sys.stderr, "running gwemopt_run..."
        print >> sys.stderr, "version: %s"%__version__
        print >> sys.stderr, ""
        print >> sys.stderr, "***************** PARAMETERS ********************"
        for o in opts.__dict__.items():
          print >> sys.stderr, o[0]+":"
          print >> sys.stderr, o[1]
        print >> sys.stderr, ""

    return opts

# =============================================================================
#
#                                    MAIN
#
# =============================================================================

warnings.filterwarnings("ignore")

# Parse command line
opts = parse_commandline()
inputDir = opts.inputDir

survey_fields = os.path.join(inputDir,'ZTF_Fields.txt')
fields = np.genfromtxt(survey_fields, comments='%', usecols=range(3), names=['field_id', 'ra', 'dec'])
rows = np.vstack((fields['field_id'], fields['ra'], fields['dec'])).T
fields = Table(rows=rows, names=('field_id', 'ra', 'dec'))

try:
    names = ('field_id', 'filter_id', 'time', 'limmag')
    table = Table.read(opts.fields, format='ascii', names=names, data_start=0)
except:
    names = ('field_id', 'filter_id', 'time', 'limmag', 'exposure_time')
    table = Table.read(opts.fields, format='ascii', names=names, data_start=0)

if opts.doTransients:
    try:
        transients = Table.read(opts.transients, format='ascii', data_start=0)
    except:
        with open(opts.transients) as f:
            content = f.readlines()
        lines = [x.strip() for x in content]  
        names, ras, decs = [], [], []
        for line in lines:
            lineSplit = list(filter(None,line.split(" ")))
            names.append(lineSplit[0])
            ras.append(Angle(lineSplit[2],unit=u.hour).deg)
            decs.append(Angle(lineSplit[4],unit=u.deg).deg)
        transients = Table([names,ras,decs], names=('col1', 'col2', 'col3'))

start_time = Time(opts.gps, format='gps').mjd
table["mjd"] = Time(table["time"], format='isot').mjd - start_time
table.sort('mjd')

if opts.doToOOnly:
    table = table[table['exposure_time']>30]

#center = SkyCoord.from_name('NGC 4993')
center = SkyCoord('14.8h', 3.5, unit=['hourangle','deg'])
center2 = SkyCoord('6h', -33.5, unit=['hourangle','deg'])
center3 = SkyCoord('12.2h', -16.5, unit=['hourangle', 'deg'])

fig = plt.figure(figsize=(8, 8), dpi=100)

ax = plt.axes(
    [0.1, 0.3, 0.9, 0.4],
    projection='astro hours mollweide')
    #center=center)

#ax_inset.compass(0.25, 0.1, 0.15)

skymap, metadata = fits.read_sky_map(opts.skymap, nest=None, distances=True)
distmu = skymap[1]
distsigma = skymap[2]
distnorm = skymap[3]
skymap = skymap[0]

skymap_ring = hp.pixelfunc.reorder(skymap, inp='NESTED', out='RING')

#calculate the moments from distmu, distsigma and distnorm
mom_mean, mom_std, mom_norm = distance.parameters_to_moments(distmu,distsigma)

nside_map = hp.npix2nside(len(skymap))
sort_idx = np.argsort(skymap)[::-1]
csm = np.empty(len(skymap))
csm[sort_idx] = np.cumsum(skymap[sort_idx])

cls = 100 * postprocess.find_greedy_credible_levels(skymap)
cs = ax.contour_hpx(
    (cls, 'ICRS'), nested=metadata['nest'],
    colors='k', linewidths=0.5, levels=[90], zorder=0)
fmt = r'%g\%%' if rcParams['text.usetex'] else '%g%%'
#ax.clabel(cs, fmt=fmt, fontsize=4, inline=True)

ax.imshow_hpx(opts.skymap, cmap='cylon')
ax.grid()

insets = []

spc = 0.325

ra = ax.coords[0]
dec = ax.coords[1]

ra.set_ticks_visible(False)
dec.set_ticks_visible(False)

nside = nside_map
ras, decs = [], []
ipixs, ipixs_cumulative, areas = {}, {}, {}
ipix_count = np.zeros(hp.pixelfunc.nside2npix(nside))

bands = {1: 'g', 2: 'r', 3: 'i', 4: 'z', 5: 'J'}
if opts.telescope == "ZTF":
    filts = ["g","r"]
elif opts.telescope == "Gattini":
    filts = ["J"]

for filt in filts:
    ipixs[filt] = []
    ipixs_cumulative[filt] = []
    areas[filt] = []

s = requests.Session()
retries = Retry(total=5,
                backoff_factor=0.1,
                status_forcelist=[ 500, 502, 503, 504 ])
s.mount('http://', HTTPAdapter(max_retries=retries))

loc = EarthLocation.of_site('Palomar')

airmass, alt = [], []
for row in table:
    field = row['field_id']
    filter_id = row['filter_id']
    r = s.get(f'http://skipper.caltech.edu:8081/telescope/%s/field/{field}/json' % opts.telescope)
    res = r.json()

    idx = np.where(fields['field_id'] == field)[0][0]
    ra, dec = fields['ra'][idx], fields['dec'][idx]

    coords = np.asarray(res['geometry']['coordinates']).T * u.deg
    coords = SkyCoord(np.squeeze(coords).T, unit='deg',
                      location=loc,
                      obstime=Time(row['time'],format='isot'))
    ras.append(coords.ra.deg)
    decs.append(coords.dec.deg)
    ax.plot(coords.ra, coords.dec, color='k', transform=ax.get_transform('world'), zorder=1)

    alt.append(np.mean(coords.altaz.az.deg))
    airmass.append(np.mean(coords.altaz.secz))

    #print(field, filter_id, np.mean(coords.altaz.az.deg), np.mean(coords.altaz.alt.deg), np.mean(coords.altaz.secz))

    if opts.doChipGaps:
        ipix = gwemopt.ztf_tiling.get_quadrant_ipix(nside, ra, dec)
        ipix = list({i for _ in ipix for i in _})
        ipix = np.array(ipix)

    else:
        xyz = []
        for r, d in zip(coords.ra.deg, coords.dec.deg):
            xyz.append(hp.ang2vec(r, d, lonlat=True))
        xyz = xyz[:4]
        try:
            ipix = hp.query_polygon(nside, np.array(xyz))
        except:
            ipix = []

    ipix_count[ipix] = ipix_count[ipix] + 1
    for filt in filts:
        if filt == bands[filter_id]:
            ipixs[filt].append(ipix)
        else:
            ipixs[filt].append([])
        ipix_cumulative = list({i for _ in ipixs[filt] for i in _})
        ipixs_cumulative[filt].append(ipix_cumulative)
        areas[filt].append(hp.nside2pixarea(nside, degrees=True) * len(ipix_cumulative))

plotName = opts.output.replace(".pdf","_numobs.pdf")
hp.mollview(ipix_count,title='')
plt.show()
plt.savefig(plotName,dpi=200)
plt.close('all')

idx = np.where(ipix_count >= 1)[0]
print('Cumulative coverage: 1+:', np.sum(skymap_ring[idx]))
idx = np.where(ipix_count >= 2)[0]
print('Cumulative coverage: 2+:', np.sum(skymap_ring[idx]))

npix = hp.nside2npix(nside)
theta, phi = hp.pix2ang(nside, np.arange(npix))
radecs = SkyCoord(ra=phi*u.rad, dec=(0.5*np.pi - theta)*u.rad)
idx = np.where(np.abs(radecs.galactic.b.deg) <= 10.0)[0]
skymap_ring[idx] = 0.0

idx = np.where(ipix_count >= 1)[0]
print('Cumulative coverage (excluding galaxy +-10): 1+:', np.sum(skymap_ring[idx]))
idx = np.where(ipix_count >= 2)[0]
print('Cumulative coverage (excluding galaxy +-10): 2+:', np.sum(skymap_ring[idx]))

table["airmass"] = airmass
table["altitude"] = alt

table["ras"] = ras
table["decs"] = decs    
for filt in filts:
    table["ipixs_%s" % filt] = ipixs[filt]
    table["ipixs_cumulative_%s" % filt] = ipixs_cumulative[filt]
    table["areas_%s" % filt] = areas[filt]

if opts.doTransients:    
    text_str = ""
    ax.scatter(transients['col2'], transients['col3'], transform=ax.get_transform('world'), color='w',zorder=2, s=70)
    if not opts.telescope == "Gattini":
        for ii, transient in enumerate(transients):
            ipix = hp.ang2pix(nside_map, transient['col2'], transient['col3'], lonlat=True)
            #ax.text(transient['col2']+2.0, transient['col3'], transient['col1'][-3:], transform=ax.get_transform('world'), color='k', fontsize=4, zorder=3)
            if ii == 11:
                ax.text(transient['col2']+3.5, transient['col3']-1.5, ii+1, transform=ax.get_transform('world'), color='k', fontsize=8, zorder=3)
            elif ii >= 9:
                ax.text(transient['col2']+3.5, transient['col3']-1.5, ii+1, transform=ax.get_transform('world'), color='k', fontsize=8, zorder=3)
            else:
                ax.text(transient['col2']+2.5, transient['col3']-1.5, ii+1, transform=ax.get_transform('world'), color='k', fontsize=8, zorder=3)
            text_str = "%s %d: %s,"%(text_str, ii+1, transient['col1'])
            #text_str = "%s %d: %s (%.1f),"%(text_str, ii+1, transient['col1'], 100*csm[ipix])
    print(text_str)

outputDir = "/".join(opts.output.split("/")[:-1])
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)
    
plt.savefig(opts.output,dpi=200,bbox_inches='tight')    
plt.close()    

plotName = opts.output.replace(".pdf","_airmass.pdf")
plt.figure()
plt.scatter(table["airmass"], table["limmag"], s=20, c=table["mjd"])
plt.xlabel('Airmass')
plt.ylabel('Limiting Magnitude [ab]')
cbar = plt.colorbar()
cbar.set_label('')
plt.savefig(plotName,dpi=200,bbox_inches='tight')
plt.close()

if opts.doMovie:
    moviedir = os.path.join(outputDir,'movie')
    if not os.path.isdir(moviedir):
        os.makedirs(moviedir)
    
    for ii in range(len(table)):
        fig = plt.figure(figsize=(8, 12))
    
        gs = fig.add_gridspec(3, 1)
        ax1 = fig.add_subplot(gs[0:2, 0], projection='astro hours mollweide')
        ax2 = fig.add_subplot(gs[2, 0])
        #ax3 = fig.add_subplot(gs[3, 0])
    
        #ax = plt.axes(
        #    [0.1, 0.3, 0.9, 0.4],
        #    projection='astro hours mollweide')
        #    #center=center)
    
        #ax_inset.compass(0.25, 0.1, 0.15)
    
        cs = ax1.contour_hpx(
            (cls, 'ICRS'), nested=metadata['nest'],
            colors='k', linewidths=0.5, levels=[90])
    
        ax1.imshow_hpx(opts.skymap, cmap='cylon')
        ax1.grid()
    
        insets = []
    
        ra = ax1.coords[0]
        dec = ax1.coords[1]
    
        ra.set_ticks_visible(False)
        dec.set_ticks_visible(False)
    
        for jj, row in enumerate(table):
            ras, decs = row['ras'], row['decs']
            if ii == jj:
                if row['filter_id'] == 1:
                    color = 'g'
                elif row['filter_id'] == 2:
                    color = 'r'
                else:
                    color = 'k'
                poly = patches.Polygon(np.vstack((ras, decs)).T, transform=ax1.get_transform('world'), alpha=0.5, color=color) 
                ax1.plot(ras, decs, color=color, transform=ax1.get_transform('world'))
                ax1.add_patch(poly)
            elif ii > jj:
                ax1.plot(ras, decs, color='k', transform=ax1.get_transform('world'))
    
        for filt in filts:
            for jj, row in enumerate(table):
                if ii >= jj:
                    mjd = row['mjd']
                    area = row['areas_%s' % filt]
                    if filt == 'g':
                        color = 'g'
                    elif filt == 'r':
                        color = 'r'
                    else:
                        color = 'k'                
                    ax2.plot(mjd, area, '*', color=color)
        ax2.set_xlim([0,np.max(table['mjd'])])
        max_area = np.max([np.max(table['areas_%s' % filt]) for filt in filts])
        ax2.set_ylim([0,max_area])
    
        ax2.set_xlabel('Time since event [days]')
        ax2.set_ylabel('Sky area [sq. deg.]')
    
        #for jj, row in enumerate(table):
        #    if ii >= jj:
        #        mjd = row['mjd']
        #        airmass = row['airmass']

        #        if row['filter_id'] == 1:
        #            color = 'g'
        #        elif row['filter_id'] == 2:
        #            color = 'r'
        #        else:
        #            color = 'k'
        #        ax3.plot(mjd, airmass, '*', color=color)
        #ax3.set_xlim([0,np.max(table['mjd'])])
        #ax3.set_ylim([0.9*np.min(table['airmass']),1.1*np.max(table['airmass'])])
        #ax3.invert_yaxis()

        #ax3.set_xlabel('Time since event [days]')
        #ax3.set_ylabel('Airmass')

        fig.suptitle('')
    
        plotName = os.path.join(moviedir,'movie-%04d.png'%ii)
        plt.savefig(plotName,dpi=200,bbox_inches='tight')
        plt.close()
    
    output = opts.output.split("/")[-1].replace(".pdf","")
    
    moviefiles = os.path.join(moviedir,"movie-%04d.png")
    filename = os.path.join(moviedir,"%s.mpg" % (output))
    ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
    os.system(ffmpeg_command)
    filename = os.path.join(moviedir,"%s.gif" % (output))
    ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
    os.system(ffmpeg_command)
    rm_command = "rm %s/*.png"%(moviedir)
    os.system(rm_command)
