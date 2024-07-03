
import os
import optparse
import warnings
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import numpy as np
import healpy as hp
import pandas as pd

from astropy import coordinates
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.time import Time, TimeDelta
from astropy.table import Table
import ligo.skymap.plot
from ligo.skymap.io import fits
from ligo.skymap import postprocess
import ligo.skymap.distance as ligodist

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
import pyvo.dal

import gwemopt.gracedb, gwemopt.ztf_tiling

__version__ = 1.0

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-f", "--observations", help="observation file.", default='../data/serendipitous/2019-04-25-2019-10-03_subdata.txt')
    parser.add_option("-o", "--outputDir", help="output directory.", default='../output/hist')

    parser.add_option("--doHists",  action="store_true", default=False)

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

def get_ztf_quadrants():
    """Calculate ZTF quadrant footprints as offsets from the telescope
    boresight."""
    quad_prob = gwemopt.ztf_tiling.QuadProb(0, 0)
    ztf_tile = gwemopt.ztf_tiling.ZTFtile(0, 0)
    quad_cents_ra, quad_cents_dec = ztf_tile.quadrant_centers()
    offsets = np.asarray([
        quad_prob.getWCS(
            quad_cents_ra[quadrant_id],
            quad_cents_dec[quadrant_id]
        ).calc_footprint(axes=quad_prob.quadrant_size)
        for quadrant_id in range(64)])
    return np.transpose(offsets, (2, 0, 1))



# =============================================================================
#
#                                    MAIN
#
# =============================================================================

warnings.filterwarnings("ignore")

# Parse command line
opts = parse_commandline()

outputDir = opts.outputDir
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)

all_obs_table = pd.read_csv(opts.observations, delimiter = '|')
all_obs_table.columns = [c.replace(' ', '') for c in all_obs_table.columns]
use_obs = all_obs_table['status'] == 1
obstable = all_obs_table[use_obs]
obstable['ccdid'] = np.floor(obstable["rcid"]/4)

obstable = obstable[['field', 'rcid', 'fid', 'expid', 'jd', 'exptime', 'sciinpseeing', 'diffmaglim', 'scimaglim', 'programid']]
names = ('field', 'rcid', 'diffmaglim', 'scimaglim', 'fid', 'programid')
renames = ('fieldid', 'chid', 'limMag', 'limMag1', 'filterid', 'progid')
obstable = obstable.rename(columns = {'field': 'fieldid', 'rcid': 'chid', 'diffmaglim': 'limMag', 'fid': 'filterid', 'programid': 'progid'})

obstable = Table.from_pandas(obstable)

data = {}
for exptime in list(set(obstable["exptime"])):
    data[exptime] = {}
    for filterid in list(set(obstable["filterid"])):

        idx = obstable["exptime"] == exptime
        obstable2 = obstable[idx]

        idx = np.where(obstable2["filterid"] == filterid)[0]
        obstable_slice = obstable2[idx]

        data[exptime][filterid] = obstable_slice["limMag"]

        try:
            bins = np.linspace(np.min(obstable_slice["limMag"]),np.max(obstable["limMag"]),50)
            hist1, bin_edges = np.histogram(obstable_slice["limMag"], bins=bins)
            hist1 = hist1 / float(np.sum(hist1))
            bins1 = (bins[1:] + bins[:-1])/2.0
        except:
            bins1 = [0, 0]
            hist1 = [0, 0]
    
        if opts.doHists: 
            bins = np.linspace(np.min(obstable_slice["limMag"]),np.max(obstable["limMag"]),50)
            hist2, bin_edges = np.histogram(obstable_slice["limMag"], bins=bins)
            hist2 = hist2 / float(np.sum(hist2))
            bins2 = (bins[1:] + bins[:-1])/2.0
            
            neg2 = np.percentile(obstable_slice["limMag"],2.5)
            neg1 = np.percentile(obstable_slice["limMag"],16)
            med = np.percentile(obstable_slice["limMag"],50)
            pos1 = np.percentile(obstable_slice["limMag"],84)
            pos2 = np.percentile(obstable_slice["limMag"],97.5)
        
            print(exptime, neg2, neg1, med, pos1, pos2)
            
            plt.figure(figsize=(12,10))
            plt.plot(bins1, hist1,'g')
            plt.plot(bins2, hist2,'r')
            plt.xlabel('Limiting Magnitude')
            plt.ylabel('Probability Density Function')
            plt.plot([neg2,neg2],[0,1],'k--')
            plt.plot([neg1,neg1],[0,1],'k--')
            plt.plot([med,med],[0,1],'k--')
            plt.plot([pos1,pos1],[0,1],'k--')
            plt.plot([pos2,pos2],[0,1],'k--')
            plt.xlim([18.0,22.0])
            plt.ylim([0,0.2])
            plt.show()
            plotName = os.path.join(outputDir,'limmag_%d.pdf' % exptime)
            plt.savefig(plotName,dpi=200)
            plt.close('all')

labels = []
tt = np.linspace(1, 320, 100)

keys = list(data.keys())
fig = plt.figure(figsize=(24,18))
ax = fig.add_subplot(111)
for ii, key in enumerate(keys):
    data_out = data[key]
    filter_ids = list(data_out.keys())
    for jj, filter_id in enumerate(filter_ids):
        if filter_id == 1:
            color = 'b'
            label = 'g'
            dt = -5
        elif filter_id == 2:
            color = 'springgreen'
            label = 'r'
            dt = 0
        elif filter_id == 3:
            color = 'red'
            label = 'i'
            dt = 5

        if len(data_out[filter_id]) == 0:
            continue
        if key == 20: continue

        parts = plt.violinplot(data_out[filter_id],[key+dt],widths=5,
                               quantiles=[0.1,0.9], showextrema=False)

        #for partname in ('cbars','cmins','cmaxes'):
        #    vp = parts[partname]
        #    vp.set_edgecolor(color)
        #    vp.set_linewidth(3)
        #for partname in ['cquantiles']:
        #    vp = parts[partname]
        #    vp.set_edgecolor(color)
        #    vp.set_linewidth(3)
        for pc in parts['bodies']:
            pc.set_facecolor(color)
            pc.set_edgecolor(color)

        perc50 = np.percentile(data_out[filter_id],50)
        if not filter_id in labels:
            plt.plot([key+dt-2,key+dt+2],[perc50,perc50],'-',color=color,label=label,linewidth=3)
            labels.append(filter_id)
        else:
            plt.plot([key+dt-2,key+dt+2],[perc50,perc50],'-',color=color,linewidth=3)

        perc90 = np.percentile(data_out[filter_id],90)
        plt.plot([key+dt-2,key+dt+2],[perc90,perc90],'-',color=color,linewidth=3)

        perc10 = np.percentile(data_out[filter_id],10)
        plt.plot([key+dt-2,key+dt+2],[perc10,perc10],'-',color=color,linewidth=3)

        percbad = float(len(np.where(data_out[filter_id] <= (perc50-2))[0]))/float(len(data_out[filter_id]))

        if key == 30:
            mag = -2.5*np.log10(np.sqrt(30/tt))
            plt.plot(tt, mag+perc50, '--', color=color,linewidth=3)

        print(key, label, perc50, percbad)

plt.xlabel('Exposure Time [s]',fontsize=48)
plt.ylabel('Limiting Magnitude',fontsize=48)
plt.grid()
plt.xlim([10, 320])
plt.ylim([17, 22])
plt.legend(loc=4, prop={'size':36})
plt.yticks(fontsize=48)
plt.xticks(fontsize=48)
plotName = os.path.join(outputDir,'limmag.pdf')
plt.savefig(plotName)
plt.close()
 
