
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
from astropy.table import Table, unique

__version__ = 1.0

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-s", "--skymap", help="GW skymap.", default='../data/GW190426/LALInference1.fits.gz,../data/bn191017391/glg_healpix_all_bn191017391.fit')
    parser.add_option("-f","--fields", help="Observed fields.", default='../data/serendipitous/2019-10-05-2019-10-11_normal.txt')

    parser.add_option("-o", "--output", help="output file",default="../output/serendipitous/overall.pdf")

    parser.add_option("-t", "--telescope", help="Telescope.", default ="ZTF")

    parser.add_option("--time", help="start time.", default="2019-10-05T00:00:00")

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

def rotate_map(hmap, rot_theta, rot_phi):
    """
    Take hmap (a healpix map array) and return another healpix map array 
    which is ordered such that it has been rotated in (theta, phi) by the 
    amounts given.
    """
    nside = hp.npix2nside(len(hmap))

    # Get theta, phi for non-rotated map
    t,p = hp.pix2ang(nside, np.arange(hp.nside2npix(nside))) #theta, phi

    # Define a rotator
    r = hp.Rotator(deg=False, rot=[rot_phi,rot_theta])

    # Get theta, phi under rotated co-ordinates
    trot, prot = r(t,p)

    # Interpolate map onto these co-ordinates
    rot_map = hp.get_interp_val(hmap, trot, prot)

    return rot_map

# =============================================================================
#
#                                    MAIN
#
# =============================================================================

warnings.filterwarnings("ignore")

# Parse command line
opts = parse_commandline()

center = SkyCoord('6.0h', 30.5, unit=['hourangle','deg'])
center2 = SkyCoord('1.0h', -24.0, unit=['hourangle','deg'])

fig = plt.figure(figsize=(8, 8), dpi=100)

#ax = plt.axes(
#    [0.1, 0.3, 0.9, 0.4],
#    projection='astro hours mollweide')
ax = plt.axes(
    [0.1, 0.3, 0.9, 0.4],
    projection='astro globe',
    center=center)

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
    skymap, metadata = fits.read_sky_map(skymap, nest=True)
    nside = hp.get_nside(skymap)
    npix = hp.nside2npix(nside)
    skymap = skymap[hp.ring2nest(nside, np.arange(npix))]
    if ii == 0:
        skymap = rotate_map(skymap, np.deg2rad(0), np.deg2rad(120))
    else:
        skymap = rotate_map(skymap, np.deg2rad(0), np.deg2rad(-55))
    skymap = skymap[hp.nest2ring(nside, np.arange(npix))]

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
color_names=[color3,color4,color5]
for ii, (fields, color_name) in enumerate(zip(opts.fields.split(","),color_names)):
    color_name=color3
    table = Table.read(fields, format='ascii.fixed_width',
                       data_start=2, data_end=-1)

    survey_start = Time(opts.time, format='isot', scale='utc').jd
    dt = 2.0
    survey_end = survey_start + dt

    table =table[table["jd"]>=survey_start]
    table =table[table["jd"]<=survey_end]

    names = ('field', 'fid')
    renames = ('field_id', 'filter_id')
    table.rename_columns(names, renames)
    table = unique(table, keys=['jd', 'field_id'])
    table = table.group_by('field_id').groups

    if ii == 3:
        table.add_row([1288, 3, '2019-08-15T03:45:27.000014', 20.0])
        table.add_row([1241, 3, '2019-08-15T03:45:27.000014', 20.0])
        table.add_row([202, 3, '2019-08-15T03:45:27.000014', 20.0])

    field_ids = []
    cnts = [0,0,0]
    labels = {}
    labels[0] = "1 Observation"
    labels[1] = "2 Observations"
    labels[2] = "3+ Observations"

    for row in table:
        field = row['field_id'][0]
        jj = len(row) - 1
        if jj > 2:
            jj = 2
        color_name = color_names[jj]
        if field in field_ids: continue

        filter_id = row['filter_id']
        r = s.get(f'http://skipper.caltech.edu:8081/telescope/%s/field/{field}/json' % opts.telescope)
        res = r.json()

        coords = np.asarray(res['geometry']['coordinates']).T * u.deg
        coords = SkyCoord(np.squeeze(coords).T, unit='deg')

        #idx1 = np.where(coords.ra.deg>=180.0)[0]
        #idx2 = np.where(coords.ra.deg<180.0)[0]
        #if (len(idx1)>0 and len(idx2)>0):
        #    continue

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
            if cnts[jj] == 0:
                cs = ax.plot(coords.ra, coords.dec, color=color_name, transform=ax.get_transform('world'), zorder=1, label=labels[jj], linewidth=linewidth)
                lines.append(cs[0])
                cnts[jj] = 1
            else:
                cs = ax.plot(coords.ra, coords.dec, color=color_name, transform=ax.get_transform('world'), zorder=1, linewidth=linewidth)
            alpha = 0.075
            path = np.vstack((coords.ra.deg, coords.dec.deg))
            path = matplotlib.path.Path(path.T)
            patch = patches.PathPatch(path, alpha=alpha, color=color_name, fill=True, zorder=3, edgecolor=color_name, transform=ax.get_transform('world'))
            ax.add_patch(patch)

        field_ids.append(field)
 
ax.legend(loc='upper left', fontsize='xx-small', bbox_to_anchor=(0,1.16), frameon=False).set_zorder(100)

outputDir = "/".join(opts.output.split("/")[:-1])
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)
    
plt.savefig(opts.output,dpi=200,bbox_inches='tight')    
plt.close()    

