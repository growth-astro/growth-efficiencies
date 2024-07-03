#!/usr/bin/env python

import os
import optparse
import warnings
from astropy.table import Table
from astropy.coordinates import Distance
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline as Spline2d
from scipy.interpolate import interpolate as interp

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time, TimeDelta
from ligo.skymap.io import fits
from ligo.skymap.distance import parameters_to_marginal_moments
#import pysqlite3
from astropy import coordinates
from astropy.cosmology import Planck15
import simsurvey
import sncosmo
#from ampel.ztf.archive.ArchiveDB import ArchiveDB
#import simsurvey_tools as sst
import datetime
import pickle
import pyvo.dal

import healpy as hp
import pandas as pd

import matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rc('font', **font)
import matplotlib.cm as cm

import gwemopt.ztf_tiling

__version__ = 1.0

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-s", "--skymap", help="GW skymap.", default='../data/GW190425/LALInference.fits.gz')
    parser.add_option("-t", "--observation_type", help="observation type.", default='obsfile')
    parser.add_option("-f", "--observations", help="observation file.", default='../data/serendipitous/2019-04-25-2019-10-03_subdata.txt')
    parser.add_option("-o", "--outputDir", help="output file",default="../output/GW190425/simsurvey/")
    parser.add_option("-i", "--inputDir", help="input directory",default="../input/")

    parser.add_option("--start_time", help="start time.", default='2019-04-25T08:18:25')
    parser.add_option("-d", "--dt", help="length of time to check.", default=3.0, type=float)

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

# Parse command line
opts = parse_commandline()

def clean_df(df, zp=30):

    df = df[df.status!=0]

    #df['skynoise'] = [10 ** (-0.4 * (i - zp)) / 5. for i in df['diffmaglim']]   # df['limMag']]
    df["filter"]  = [rename_filter(i) for i in df['filterid']] #FIXME ONCE WE HAVE FID DATA
    #df = df[df["fieldid"]<1896]
    #df = df[df["limMag"]>19] #not sure abut this; it will reduce LC points on bright SNe
    print('Survey pointings for all programs:')
    print(len(df))
    #df = df[df["progid"]==1] #MSIP only
    print('Survey pointings for MSIP programs:')
    print(len(df[df["progid"]==1]))
    #survey_start = Time("2018-06-01 00:00:00.000").jd #Survey start
    #survey_end = Time("2018-12-31 00:00:00.000").jd
    #df =df[df["jd"]>survey_start]
    #df =df[df["jd"]<survey_end]
    return df

def rename_filter(fid):
    if fid == 1:
        return "ztfg"
    if fid == 2:
        return "ztfr"
    if fid==3:
        return "ztfi"

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

def compute_ccdid_from_channel(chid):
    #this maps eeach channel to CCDID 0 -15 to match simulation definition
    if chid in list(range(0, 4)):
        return 0
    elif chid in list(range(4, 8)):
        return 1
    elif chid in list(range(8, 12)):
        return 2
    elif chid in list(range(12, 16)):
        return 3
    elif chid in list(range(16, 20)):
        return 4
    elif chid in list(range(20, 24)):
        return 5
    elif chid in list(range(24, 28)):
        return 6
    elif chid in list(range(28, 32)):
        return 7
    elif chid in list(range(32, 36)):
        return 8
    elif chid in list(range(36, 40)):
        return 9
    elif chid in list(range(40, 44)):
        return 10
    elif chid in list(range(44, 48)):
        return 11
    elif chid in list(range(48, 52)):
        return 12
    elif chid in list(range(52, 56)):
        return 13
    elif chid in list(range(56, 60)):
        return 14
    elif chid in list(range(60, 64)):
        return 15

# =============================================================================
#
#                                    MAIN
#
# =============================================================================

warnings.filterwarnings("ignore")

client = pyvo.dal.TAPService('https://irsa.ipac.caltech.edu/TAP')

#set up directory
inputDir = opts.inputDir
outputDir = opts.outputDir
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)

survey_fields = os.path.join(inputDir,'ZTF_Fields.txt')

pcklFile = os.path.join(outputDir,'df_sim.pkl')
start_time = Time(opts.start_time, format='isot', scale='utc')
end_time = start_time + TimeDelta(opts.dt * u.day)

if not os.path.isfile(pcklFile):
    if opts.observation_type == "TAP":
        obstable = client.search("""
        SELECT field,rcid,fid,expid,obsjd,exptime,seeing,airmass,maglimit,ipac_gid
        FROM ztf.ztf_current_meta_sci WHERE (obsjd BETWEEN {0} AND {1})
        """.format(start_time.jd, end_time.jd)).to_table()

        obstable['ccdid']=np.floor(obstable["rcid"]/4)
        names = ('obsjd', 'field', 'rcid', 'maglimit', 'fid', 'ipac_gid')
        renames = ('jd', 'fieldid', 'chid', 'limMag', 'filterid', 'progid')
        obstable.rename_columns(names, renames)
    elif opts.observation_type == "obsfile":
        obstable = Table.read(opts.observations, format='ascii.fixed_width',
                              data_start=2, data_end=-1)
        obstable['ccdid']=np.floor(obstable["rcid"]/4).astype(int)
        #names = ('field', 'rcid', 'scimaglim', 'fid', 'programid', 'subtractionstatus')
        names = ('field', 'rcid', 'diffmaglim', 'fid', 'programid', 'subtractionstatus')
        renames = ('fieldid', 'chid', 'limMag', 'filterid', 'progid', 'status')
        obstable.rename_columns(names, renames)
    else:
        print("observation type not understood")
        exit(0)

    obstable = obstable.to_pandas()
    obstable =obstable[obstable["jd"]>start_time.jd]
    obstable =obstable[obstable["jd"]<end_time.jd]

    obstable.to_pickle(pcklFile)

df = pd.read_pickle(pcklFile)
df = clean_df(df)

skymap, metadata = fits.read_sky_map(opts.skymap, nest=False, distances=True)
map_struct = {}
map_struct["prob"] = skymap[0]
map_struct["distmu"] = skymap[1]
map_struct["distsigma"] = skymap[2]
map_struct["distnorm"] = skymap[3]

nside = hp.get_nside(skymap[0])
npix = hp.nside2npix(nside)

sort_idx = np.argsort(map_struct["prob"])[::-1]
csm = np.empty(len(map_struct["prob"]))
csm[sort_idx] = np.cumsum(map_struct["prob"][sort_idx])
#map_struct["prob"][csm>0.9] = 0.0

fields = np.genfromtxt(survey_fields, comments='%', usecols=range(3), names=['field_id', 'ra', 'dec'])
rows = np.vstack((fields['field_id'], fields['ra'], fields['dec'])).T
fields = Table(rows=rows, names=('field_id', 'ra', 'dec'))

quadrant_coords = get_ztf_quadrants()
skyoffset_frames = coordinates.SkyCoord(fields['ra'], fields['dec'], unit=u.deg).skyoffset_frame()

quadrant_coords_icrs = coordinates.SkyCoord(
    *np.tile(quadrant_coords[:, np.newaxis, ...], (len(fields), 1, 1)), unit=u.deg,
    frame=skyoffset_frames[:, np.newaxis, np.newaxis]).transform_to(coordinates.ICRS)
quadrant_xyz = np.moveaxis(quadrant_coords_icrs.cartesian.xyz.value, 0, -1)

field_ids = np.unique(df.fieldid)
probs = {}
fieldlist_1, fieldlist_2 = [], []
probstot_1, probstot_2 = [], []
ipixs = {}
ccds = []
for field_id, xyz in zip(fields['field_id'], quadrant_xyz):
    if not field_id in probs:
        probs[field_id] = []
    if not field_id in ipixs:
        ipixs[field_id] = []
    for ii, xyz in enumerate(xyz):
        ipix = hp.query_polygon(nside, xyz)
        probs[field_id].append(np.sum(map_struct["prob"][ipix]))
        ipixs[field_id].append(ipix)
    probs[field_id] = np.array(probs[field_id])
    if field_id < 1000:
        fieldlist_1.append(field_id)
        probstot_1.append(np.sum(probs[field_id]))
    else:
        fieldlist_2.append(field_id)
        probstot_2.append(np.sum(probs[field_id]))

fieldlist_1, probstot_1 = np.array(fieldlist_1), np.array(probstot_1)
fieldlist_2, probstot_2 = np.array(fieldlist_2), np.array(probstot_2)

sort_idx = np.argsort(probstot_1)[::-1]
csm_1 = np.empty(len(probstot_1))
csm_1[sort_idx] = np.cumsum(probstot_1[sort_idx])

sort_idx = np.argsort(probstot_2)[::-1]
csm_2 = np.empty(len(probstot_2))
csm_2[sort_idx] = np.cumsum(probstot_2[sort_idx])

fieldlist_1 = fieldlist_1[csm_1<0.8]
fieldlist_2 = fieldlist_2[csm_2<0.8]
print(fieldlist_1, fieldlist_2)

filename = os.path.join(outputDir, 'fields.dat')
fid = open(filename, 'w')

df_group = df.groupby('jd')
rows_to_concat, limmags, ccdids_all = [], [], []
for ii, (g, data) in enumerate(df_group):
    if np.mod(ii, 100) == 0:
        print('Process group %d/%d' % (ii, len(df_group)))
    df_sub = df_group.get_group(g)
    ccdids = df_sub["ccdid"].tolist()
    limmag_median = np.median(df_sub["limMag"].tolist())

    fieldid = df_sub["fieldid"].tolist()[0]
    filterid = df_sub["filterid"].tolist()[0]
    obstime = Time(df_sub["jd"].tolist()[0], format='jd').isot
    exposuretime = df_sub["exptime"].tolist()[0]

    if (np.isin(fieldid, fieldlist_1) or np.isin(fieldid, fieldlist_2)):
        fid.write('%d %d %s %.5f %d\n' % (fieldid, filterid,
                                          obstime, limmag_median,
                                          exposuretime))
fid.close()

