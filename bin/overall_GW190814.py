
import os
import optparse
import warnings
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import numpy as np
import healpy as hp

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.time import Time
from astropy_healpix import uniq_to_level_ipix, level_to_nside
from astropy_healpix.core import boundaries_lonlat
import ligo.skymap.plot
from ligo.skymap.io import fits
from ligo.skymap import postprocess

import matplotlib
#matplotlib.rc('text', usetex=True)
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
from matplotlib import patches
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import patches

from astropy.visualization.wcsaxes import SphericalCircle
from astropy.table import Table

__version__ = 1.0

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-s", "--skymap", help="GW skymap.", default='../data/GW190814/LALInference.v1.fits.gz')
    parser.add_option("-f","--fields", help="Observed fields.", default='../data/GW190814/ZTF_fields.dat')

    parser.add_option("-o", "--output", help="output file",default="../output/overall_GW190814.pdf")

    parser.add_option("-t", "--telescope", help="Telescope.", default ="ZTF")

    parser.add_option("--doMovie",  action="store_true", default=False)

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

center = SkyCoord('1.0h', -24.0, unit=['hourangle','deg'])

fig = plt.figure(figsize=(8, 8), dpi=100)

#ax = plt.axes(
#    [0.1, 0.3, 0.9, 0.4],
#    projection='astro hours mollweide')

fig, axs = plt.subplots(
    1, 3, figsize=(12, 4),
    subplot_kw=dict(projection='astro zoom', center=center, radius=15*u.deg))

for ax in axs:
    for key in ['ra', 'dec']:
        ax.coords[key].set_ticklabel_visible(False)
        ax.coords[key].set_ticks_visible(False)

#ax = plt.axes(
#    [0.1, 0.3, 0.9, 0.4],
#    projection='astro zoom',
#    center=center,
#    radius=15*u.deg)

color2 = 'coral'
color1 = 'cornflowerblue'
color3 = 'palegreen'
color4 = 'darkmagenta'
color5 = 'red'
color_names=[color1,color2,color3,color4,color5]
color_name = color_names[0]

events = []
for skymap in opts.skymap.split(","):
    skymapsplit = skymap.split("/")
    events.append(skymapsplit[-2])

prob, metadata = fits.read_sky_map(opts.skymap)
skymap_moc = fits.read_sky_map(opts.skymap, moc=True)

cls = 100 * postprocess.find_greedy_credible_levels(prob)
kwargs = {'linestyle': '--'}
linewidths = 5

ax = axs[0]
ax.set_title('Probability Sky Map')
ax.imshow_hpx(prob, cmap='cylon')
ax.grid()

ax = axs[1]
ax.grid(color='none', visible=False)
ax.set_title('Multi-resolution HEALPix Tiles')
ax.title.set_rasterized(False)
transform = ax.get_transform('world')
level, ipix = uniq_to_level_ipix(skymap_moc['UNIQ'])
nside = level_to_nside(level)
groups = Table({'nside': nside, 'ipix': ipix}).group_by('nside').groups
for key, group in zip(groups.keys, groups):
    nside = key['nside']
    vertices = np.stack(
        boundaries_lonlat(group['ipix'], 2, nside, order='nested'),
        axis=-1
    ).deg
    for verts in vertices:
        ax.add_patch(plt.Polygon(
            verts,
            edgecolor='black',
            facecolor='none',
            linewidth=0.25 if nside <= 64 else 0.125,
            rasterized=True,
            transform=transform
        ))

ax = axs[2]
ax.set_title('Sky Map Tiling')
cs = ax.contour_hpx(
            (cls, 'ICRS'), nested=False,
            colors=color_name, linewidths=linewidths, levels=[90], zorder=10)
kwargs = {'linestyle': '-'}
cs = ax.contour_hpx(
            (cls, 'ICRS'), nested=False,
            colors=color_name, linewidths=linewidths, levels=[50], zorder=10)

fmt = r'%g\%%' if rcParams['text.usetex'] else '%g%%'
ax.grid()

insets = []

spc = 0.325

ra = ax.coords[0]
dec = ax.coords[1]

ra.set_ticks_visible(False)
dec.set_ticks_visible(False)

nside = 256
ipixs, ipixs_cumulative, areas = {}, {}, {}

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

linewidth=0.5
lines = []
for ii, (fields, color_name) in enumerate(zip(opts.fields.split(","),color_names)):
    names = ('field_id', 'filter_id', 'time', 'limmag')
    table = Table.read(fields, format='ascii', names=names, data_start=0)

    if ii == 0:
        table.add_row([1288, 3, '2019-08-15T03:45:27.000014', 20.0])
        table.add_row([1241, 3, '2019-08-15T03:45:27.000014', 20.0])
        table.add_row([202, 3, '2019-08-15T03:45:27.000014', 20.0])

    field_ids = []
    cnt = 0
    for row in table:
        field = row['field_id']
        if field in field_ids: continue

        filter_id = row['filter_id']
        r = s.get(f'http://skipper.caltech.edu:8081/telescope/%s/field/{field}/json' % opts.telescope)
        res = r.json()

        coords = np.asarray(res['geometry']['coordinates']).T * u.deg
        coords = SkyCoord(np.squeeze(coords).T, unit='deg')

        idx1 = np.where(coords.ra.deg>=180.0)[0]
        idx2 = np.where(coords.ra.deg<180.0)[0]
        if (len(idx1)>0 and len(idx2)>0):
            continue

        idx1 = np.where(np.abs(coords.dec.deg)>=87.0)[0]
        if len(idx1) == 4:
            continue

        if ii == 3:
            ax_inset.plot(coords.ra, coords.dec, color=color_name, transform=ax_inset.get_transform('world'), zorder=1, linewidth=linewidth)
            alpha = 0.075
            path = np.vstack((coords.ra.deg, coords.dec.deg))
            path = matplotlib.path.Path(path.T)
            patch = patches.PathPatch(path, alpha=alpha, color=color_name, fill=True, zorder=3, edgecolor=color_name, transform=ax_inset.get_transform('world'))
            ax_inset.add_patch(patch)
        else:
            if ii == 1:
                if coords.ra[0].deg < 16.0*360/24.0: continue
            if cnt == 0:
                cs = ax.plot(coords.ra, coords.dec, color=color_name, transform=ax.get_transform('world'), zorder=1, label=events[ii], linewidth=linewidth)
                lines.append(cs[0])
            else:
                cs = ax.plot(coords.ra, coords.dec, color=color_name, transform=ax.get_transform('world'), zorder=1, linewidth=linewidth)
            alpha = 0.075
            path = np.vstack((coords.ra.deg, coords.dec.deg))
            path = matplotlib.path.Path(path.T)
            patch = patches.PathPatch(path, alpha=alpha, color=color_name, fill=True, zorder=3, edgecolor=color_name, transform=ax.get_transform('world'))
            ax.add_patch(patch)

        field_ids.append(field)
        cnt = cnt + 1
 
ax.legend(lines, events[:-1], loc='upper left', fontsize='xx-small', bbox_to_anchor=(0,1.16), frameon=False).set_zorder(100)

outputDir = "/".join(opts.output.split("/")[:-1])
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)
    
plt.savefig(opts.output,dpi=200,bbox_inches='tight')    
plt.close()    

