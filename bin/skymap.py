
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
from matplotlib.lines import Line2D

from astropy.visualization.wcsaxes import SphericalCircle
from astropy.table import Table

__version__ = 1.0

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-s", "--skymap", help="GW skymap.", default='../data/GW190425/bayestar.fits.gz,../data/GW190425/LALInference.fits.gz')
    parser.add_option("-o", "--output", help="output file",default="../output/GW190425/skymap.pdf")

    parser.add_option("-r", "--ra", default='18.0h')
    parser.add_option("-d", "--dec", default=30.5, type=float)

    parser.add_option("-p", "--projection", default='astro globe')

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

center = SkyCoord(opts.ra, opts.dec, unit=['hourangle','deg'])
center2 = SkyCoord('1.0h', -24.0, unit=['hourangle','deg'])

fig = plt.figure(figsize=(8, 8), dpi=100)

if opts.projection == "astro zoom":
    ax = plt.axes(
        [0.1, 0.3, 0.9, 0.4],
        projection=opts.projection,
        center=center,
        radius=15*u.deg)
elif opts.projection == "astro globe":
    ax = plt.axes(
        [0.1, 0.3, 0.9, 0.4],
        projection=opts.projection,
        center=center)
else:
    ax = plt.axes(
        [0.1, 0.3, 0.9, 0.4],
        projection=opts.projection)

color2 = 'coral'
color1 = 'cornflowerblue'
color3 = 'palegreen'
color4 = 'darkmagenta'
color_names=[color1,color2,color3,color4]

skymapnames = []
eventname = []
for skymap in opts.skymap.split(","):
    skymapsplit = skymap.split("/")
    skymapname = skymapsplit[-1].replace(".fits.gz","")
    skymapname = skymapname.capitalize()
    skymapnames.append(skymapname)
    eventname = skymapsplit[-2]

legend_elements = []
for ii, (skymap, color_name) in enumerate(zip(opts.skymap.split(","),color_names)):
    skymap, metadata = fits.read_sky_map(skymap, nest=None)
    cls = 100 * postprocess.find_greedy_credible_levels(skymap)
    kwargs = {'linestyle': '--'}

    linewidths = 1.5
    if ii == 3:
        cs = ax_inset.contour_hpx(
            (cls, 'ICRS'), nested=metadata['nest'],
            colors=color_name, linewidths=linewidths, levels=[99], zorder=10)
        kwargs = {'linestyle': '-'}
        cs = ax_inset.contour_hpx(
            (cls, 'ICRS'), nested=metadata['nest'],
            colors=color_name, linewidths=linewidths, levels=[90], zorder=10)
        ax_inset.set_title(skymaps[ii], fontsize='x-small')
    else:
        cs = ax.contour_hpx(
            (cls, 'ICRS'), nested=metadata['nest'],
            colors=color_name, linewidths=linewidths, levels=[90], 
            zorder=10, linestyle='--')
        cs = ax.contour_hpx(
            (cls, 'ICRS'), nested=metadata['nest'],
            colors=color_name, linewidths=linewidths, levels=[50],
            zorder=10, linestyle='--')

    fmt = r'%g\%%' if rcParams['text.usetex'] else '%g%%'
    #ax.clabel(cs, fmt=fmt, fontsize=4, inline=True)

    #ax.imshow_hpx(opts.skymap, cmap='cylon')

    legend_elements.append(Line2D([0], [0], color=color_name, linewidth=linewidths, label=skymapnames[ii]))

ax.grid()

insets = []

spc = 0.325

ra = ax.coords[0]
dec = ax.coords[1]

ra.set_ticks_visible(False)
dec.set_ticks_visible(False)

ax.legend(handles=legend_elements, loc='upper left', fontsize='x-small').set_zorder(100)
ax.set_title(eventname)

outputDir = "/".join(opts.output.split("/")[:-1])
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)
    
plt.savefig(opts.output,dpi=200,bbox_inches='tight')    
plt.close()    

