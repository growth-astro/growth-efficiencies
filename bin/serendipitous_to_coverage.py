
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

    parser.add_option("-i", "--inputDir", help="input directory",default="../input/")
    parser.add_option("-g", "--graceid", help="gracedb id.", default='S190412m')
    parser.add_option("-s", "--skymap", help="GW skymap.", default='../data/GW190425/LALInference.fits.gz')
    parser.add_option("-t", "--observation_type", help="observation type.", default='obsfile')
    parser.add_option("-f", "--observations", help="observation file.", default='../data/serendipitous/2019-04-25-2019-10-03_subdata.txt')
    parser.add_option("-o","--outfields", help="Output fields file.", default='../data/GW190425/ZTF_fields_MSIP.dat')

    parser.add_option("--start_time", help="start time.", default='2019-04-25T08:18:25')
    parser.add_option("-d", "--dt", help="length of time to check.", default=3.0, type=float)

    parser.add_option("--doGraceDB",  action="store_true", default=False)
    parser.add_option("--doMSIPOnly",  action="store_true", default=False)

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
start_time = opts.start_time

if opts.doGraceDB:
    params = {}
    params["event"] = opts.graceid
    params["outputDir"] = "../data/%s"%opts.graceid
    if not os.path.isdir(params["outputDir"]):
        os.makedirs(params["outputDir"])
    skymapfile, eventinfo = gwemopt.gracedb.get_event(params)
    start_time = Time(eventinfo["gpstime"], format="gps")
else:
    skymapfile = opts.skymap
    start_time = Time(opts.start_time, format='isot', scale='utc')

inputDir = opts.inputDir
survey_fields = os.path.join(inputDir,'ZTF_Fields.txt')

client = pyvo.dal.TAPService('https://irsa.ipac.caltech.edu/TAP')

end_time = start_time + TimeDelta(opts.dt * u.day)

if opts.observation_type == "TAP":
    obstable = client.search("""
    SELECT field,rcid,fid,expid,obsjd,exptime,seeing,airmass,maglimit,ipac_gid
    FROM ztf.ztf_current_meta_sci WHERE (obsjd BETWEEN {0} AND {1})
    """.format(start_time.jd, end_time.jd)).to_table()

    obstable['ccdid']=np.floor(obstable["rcid"]/4)
    names = ('obsjd', 'field', 'rcid', 'maglimit', 'fid', 'ipac_gid')
    renames = ('jd', 'fieldid', 'chid', 'limMag', 'filterid', 'progid')
    obstable.rename_columns(names, renames)
    obstable = obstable.to_pandas()
    obstable =obstable[obstable["jd"]>start_time.jd]
    obstable =obstable[obstable["jd"]<end_time.jd]

elif opts.observation_type == "obsfile":
    obstable = Table.read(opts.observations, format='ascii.fixed_width',
                          data_start=2, data_end=-1)
    obstable['ccdid']=np.floor(obstable["rcid"]/4)
    names = ('field', 'rcid', 'scimaglim', 'fid', 'programid')
    renames = ('fieldid', 'chid', 'limMag', 'filterid', 'progid')
    obstable.rename_columns(names, renames)
    obstable = obstable.to_pandas()
    obstable =obstable[obstable["jd"]>start_time.jd]
    obstable =obstable[obstable["jd"]<end_time.jd]

elif opts.observation_type == "allobsfile":
    all_obs_table = pd.read_csv(opts.observations, delimiter = '|')
    all_obs_table.columns = [c.replace(' ', '') for c in all_obs_table.columns]
    use_obs = (all_obs_table['jd'] > start_time.jd) & (all_obs_table['jd'] < end_time.jd) & (all_obs_table['status'] == 1)
    obstable = all_obs_table[use_obs]
    obstable['ccdid'] = np.floor(obstable["rcid"]/4)

    obstable = obstable[['field', 'rcid', 'fid', 'expid', 'jd', 'exptime', 'sciinpseeing', 'diffmaglim', 'programid']]
    names = ('field', 'rcid', 'diffmaglim', 'fid', 'programid')
    renames = ('fieldid', 'chid', 'limMag', 'filterid', 'progid')
    obstable = obstable.rename(columns = {'field': 'fieldid', 'rcid': 'chid', 'diffmaglim': 'limMag', 'fid': 'filterid', 'programid': 'progid'})
else:
    print("observation type not understood")
    exit(0)

if opts.doMSIPOnly:
    obstable =obstable[obstable["progid"] == 1]
df = obstable

skymap, metadata = fits.read_sky_map(skymapfile, nest=False, distances=True)
map_struct = {}
map_struct["prob"] = skymap[0]
map_struct["distmu"] = skymap[1]
map_struct["distsigma"] = skymap[2]
map_struct["distnorm"] = skymap[3]

nside = 256
map_struct["prob"], map_struct["distmu"], map_struct["distsigma"], map_struct["distnorm"] = ligodist.ud_grade(map_struct["prob"], map_struct["distmu"], map_struct["distsigma"], nside)

#unit = ""
#plotName = opts.outfields.replace(".dat",".pdf")
#hp.mollview(map_struct["prob"],title='',unit=unit,min=np.percentile(map_struct["prob"],1),max=np.percentile(map_struct["prob"],99))
#plt.show()
#plt.savefig(plotName,dpi=200)
#plt.close('all')

nside = hp.get_nside(map_struct["prob"])
npix = hp.nside2npix(nside)

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

field_prob_1, field_ids_1 = [], []
field_prob_2, field_ids_2 = [], []

ipixs = {}
ipixs_all = np.empty((0,1))
totprob = 0.0
coverage = np.zeros(map_struct["prob"].shape)
for field_id, xyz in zip(fields['field_id'], quadrant_xyz):
    if not field_id in probs:
        probs[field_id] = []
    if not field_id in ipixs:
        ipixs[field_id] = []
    for ii, xyz in enumerate(xyz):
        ipix = hp.query_polygon(nside, xyz)
        probs[field_id].append(np.sum(map_struct["prob"][ipix]))
        ipixs[field_id].append(ipix)
        coverage[ipix] = coverage[ipix] + 1
        if field_id in field_ids:
            ipixs_all = np.append(ipixs_all,ipix)
    probs[field_id] = np.array(probs[field_id])
    if field_id < 1000:
        field_ids_1.append(field_id)
        field_prob_1.append(np.sum(probs[field_id]))
    else:
        field_ids_2.append(field_id)
        field_prob_2.append(np.sum(probs[field_id]))

ipix = np.unique(ipixs_all).astype(int)
prob = np.sum(map_struct["prob"][ipix])

fid = open(opts.outfields.replace(".dat",".prob"),'w')
fid.write('%.5e' % prob)
fid.close()
print('Event: %s, Integrated probability: %.5f' % (opts.graceid, prob))

field_ids_1, field_ids_2 = np.array(field_ids_1), np.array(field_ids_2)
field_prob_1 = np.array(field_prob_1)
field_prob_2 = np.array(field_prob_2)

idx = np.argsort(field_prob_1)[::-1]
field_ids_1, field_prob_1 = field_ids_1[idx], field_prob_1[idx]
idx = np.argsort(field_prob_2)[::-1]
field_ids_2, field_prob_2 = field_ids_2[idx], field_prob_2[idx]

cumsum = np.cumsum(field_prob_1)
idx2 = np.argmin(np.abs(cumsum-0.90))
field_ids_1 = field_ids_1[:idx2]
field_prob_1 = field_prob_1[:idx2]

cumsum = np.cumsum(field_prob_2)
idx2 = np.argmin(np.abs(cumsum-0.90))
field_ids_2 = field_ids_2[:idx2]
field_prob_2 = field_prob_2[:idx2]

field_ids = np.hstack((field_ids_1,field_ids_2))

jds = []
fid = open(opts.outfields,'w')
cnt = 0
for ii, observation in df.iterrows():
    if observation.fieldid not in field_ids: continue
    if observation.limMag < 0: continue
    if observation.jd not in jds:
        tt = Time(observation.jd, format='jd').isot
        fid.write('%d %d %s %.5f\n' % (observation.fieldid, observation.filterid, tt, observation.limMag))
        jds.append(observation.jd)
fid.close()
