
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

    parser.add_option("-s", "--skymap", help="GW skymap.", default='../data/GW190425/LALInference.fits.gz,../data/GW190426/LALInference1.fits.gz,../data/GW190901/LALInference.v2.fits.gz,../data/GW190814/LALInference.v1.fits.gz')
    parser.add_option("-f","--fields", help="Observed fields.", default='../data/GW190425/ZTF_fields.dat,../data/GW190426/ZTF_fields.dat,../data/GW190901/ZTF_fields.dat,../data/GW190814/ZTF_fields.dat')

    parser.add_option("-o", "--output", help="output file",default="../output/overall.pdf")

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

center = SkyCoord('18.0h', 30.5, unit=['hourangle','deg'])
center2 = SkyCoord('1.0h', -24.0, unit=['hourangle','deg'])

fig = plt.figure(figsize=(8, 8), dpi=100)

#ax = plt.axes(
#    [0.1, 0.3, 0.9, 0.4],
#    projection='astro hours mollweide')
ax = plt.axes(
    [0.1, 0.3, 0.9, 0.4],
    projection='astro globe',
    center=center)

ax_inset = plt.axes(
    [0.63, 0.65, 0.1, 0.1],
    projection='astro zoom',
    center=center2,
    radius=15*u.deg)

for key in ['ra', 'dec']:
    ax_inset.coords[key].set_ticklabel_visible(False)
    ax_inset.coords[key].set_ticks_visible(False)
ax.mark_inset_axes(ax_inset)
ax.connect_inset_axes(ax_inset, 'upper left')
ax.connect_inset_axes(ax_inset, 'lower left')

color2 = 'coral'
color1 = 'cornflowerblue'
color3 = 'palegreen'
color4 = 'darkmagenta'
color5 = 'red'
color_names=[color1,color2,color3,color4,color5]

events = []
for skymap in opts.skymap.split(","):
    skymapsplit = skymap.split("/")
    events.append(skymapsplit[-2])

lines = []
for ii, (skymap, color_name) in enumerate(zip(opts.skymap.split(","),color_names)):
    skymap, metadata = fits.read_sky_map(skymap, nest=None)
    cls = 100 * postprocess.find_greedy_credible_levels(skymap)
    kwargs = {'linestyle': '--'}

    linewidths = 1.5
    if ii == 3:
        cs = ax_inset.contour_hpx(
            (cls, 'ICRS'), nested=metadata['nest'],
            colors=color_name, linewidths=linewidths, levels=[90], zorder=10)
        kwargs = {'linestyle': '-'}
        cs = ax_inset.contour_hpx(
            (cls, 'ICRS'), nested=metadata['nest'],
            colors=color_name, linewidths=linewidths, levels=[50], zorder=10)
        ax_inset.set_title(events[ii], fontsize='xx-small')
    else:
        cs = ax.contour_hpx(
            (cls, 'ICRS'), nested=metadata['nest'],
            colors=color_name, linewidths=linewidths, levels=[90], 
            zorder=10, linestyle='--')
        lines.append(cs)
        cs = ax.contour_hpx(
            (cls, 'ICRS'), nested=metadata['nest'],
            colors=color_name, linewidths=linewidths, levels=[50],
            zorder=10, linestyle='--')

    fmt = r'%g\%%' if rcParams['text.usetex'] else '%g%%'
    #ax.clabel(cs, fmt=fmt, fontsize=4, inline=True)

    #ax.imshow_hpx(opts.skymap, cmap='cylon')
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

    if ii == 3:
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

