
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

__version__ = 1.0

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-s", "--skymap", help="GW skymap.", default='../data/GW190425/LALInference.fits.gz')
    parser.add_option("-f","--fields", help="Observed fields.", default='../data/GW190425/ZTF_fields.dat')
    parser.add_option("-o", "--output", help="output file",default="../output/GW190425/LALInf_distance.pdf")

    parser.add_option("-t", "--telescope", help="Telescope.", default ="ZTF")

    parser.add_option("-a", "--apparent", help="Apparent magnitude.", default=20.2, type=float)

    parser.add_option("--doObserved",  action="store_true", default=False)

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

try:
    names = ('field_id', 'filter_id', 'time', 'limmag')
    table = Table.read(opts.fields, format='ascii', names=names, data_start=0)
except:
    names = ('field_id', 'filter_id', 'time', 'limmag', 'exposure_time')
    table = Table.read(opts.fields, format='ascii', names=names, data_start=0)

table = table.group_by('field_id').groups

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

skymap, metadata = fits.read_sky_map(opts.skymap, nest=True, distances=True)
distmu = skymap[1]
distsigma = skymap[2]
distnorm = skymap[3]
skymap = skymap[0]
nside = hp.npix2nside(len(skymap))

s = requests.Session()
retries = Retry(total=5,
                backoff_factor=0.1,
                status_forcelist=[ 500, 502, 503, 504 ])
s.mount('http://', HTTPAdapter(max_retries=retries))

loc = EarthLocation.of_site('Palomar')

#calculate the moments from distmu, distsigma and distnorm
mom_mean, mom_std, mom_norm = distance.parameters_to_moments(distmu,distsigma)
distmod = 5.0*np.log10(mom_mean*1e6)-5.0
if opts.doObserved:
    abs_mean = np.inf*np.ones(mom_mean.shape)
    for row in table:
        field = row['field_id'][0]

        filter_id = row['filter_id'][0]
        limmag = np.max(row['limmag'])
        
        r = s.get(f'http://skipper.caltech.edu:8081/telescope/%s/field/{field}/json' % opts.telescope)
        res = r.json()    

        coords = np.asarray(res['geometry']['coordinates']).T * u.deg
        coords = SkyCoord(np.squeeze(coords).T, unit='deg',
                          location=loc,
                          obstime=Time(row['time'],format='isot'))

        xyz = []
        for r, d in zip(coords.ra.deg, coords.dec.deg):
            xyz.append(hp.ang2vec(r, d, lonlat=True))
        xyz = xyz[:4]
        try:
            ipix = hp.query_polygon(nside, np.array(xyz), nest=True)
        except:
            ipix = []
        abs_mean[ipix] = limmag - distmod[ipix]
else:
    abs_mean = opts.apparent - distmod

abs_mean[~np.isfinite(abs_mean)] = np.nanmax(abs_mean)

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

c = ax.imshow_hpx(abs_mean, cmap='cylon_r', nested=True, vmin=-18, vmax=-15)
ax.grid()

insets = []

spc = 0.325

ra = ax.coords[0]
dec = ax.coords[1]

ra.set_ticks_visible(False)
dec.set_ticks_visible(False)

cbar = plt.colorbar(c)
cbar.set_label('Absolute Magnitude')
cbar.ax.invert_yaxis()

for row in table:
    field = row['field_id'][0]
    filter_id = row['filter_id'][0]
    r = s.get(f'http://skipper.caltech.edu:8081/telescope/%s/field/{field}/json' % opts.telescope)
    res = r.json()

    coords = np.asarray(res['geometry']['coordinates']).T * u.deg
    coords = SkyCoord(np.squeeze(coords).T, unit='deg')
    ax.plot(coords.ra, coords.dec, color='gray', transform=ax.get_transform('world'), zorder=1)

outputDir = "/".join(opts.output.split("/")[:-1])
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)
    
plt.savefig(opts.output,dpi=200,bbox_inches='tight')    
plt.close()    

plotName = opts.output.replace(".pdf","_hist.pdf")
bins = np.arange(-20,-10, 0.25)
idx = np.where(np.isfinite(abs_mean))[0]
hist, bin_edges = np.histogram(abs_mean[idx], bins=bins, weights=skymap[idx], density=True)
bins = (bin_edges[1:]+bin_edges[:-1])/2.0
plt.figure()
plt.step(bins, hist, '-', color='k')
plt.gca().invert_xaxis()
plt.xlabel('Absolute Magnitude')
plt.ylabel('Probability Density')
plt.savefig(plotName,dpi=200,bbox_inches='tight')
plt.close()

