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
import random

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

import ligo.skymap.plot

import matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rc('font', **font)
import matplotlib.cm as cm
import h5py

import gwemopt.ztf_tiling

__version__ = 1.0

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-t", "--observation_type", help="observation type.", default='obsfile')
    parser.add_option("-f", "--observations", help="observation file.", default='../data/serendipitous/ztf_obsfile_status1_v2_0501_0901.csv')
    parser.add_option("-o", "--outputDir", help="output file",default="../output/GW190425/simsurvey/")
    parser.add_option("-i", "--inputDir", help="input directory",default="../input/")
    parser.add_option("-x", "--sfdDir", help="sfd directory",default="../input/sfd98/")
    parser.add_option("-p", "--possisDir", help="possis directory",default="../data/possis/")
    parser.add_option("--kasenDir", help="possis directory",default="../data/kasen/")
    parser.add_option("-b", "--bsfile", help="Bright stars loss file",default="../input/bright_star_loss.csv")
    parser.add_option("-n", "--ntransient", help="int number", default=10000, type=int)
    parser.add_option("-r", "--rate", help="kilonovae rate", default=1e-6, type=float)
    parser.add_option("-z", "--zeropoint", help="zeropoint", default=30, type=float)
    parser.add_option("--doRotate",  action="store_true", default=False)
    parser.add_option("--theta", help="theta rotation.", default=0.0, type=float)
    parser.add_option("--phi", help="phi rotation.", default=0.0, type=float)
    parser.add_option("--peakmag", help="Limit for peak magnitude.", default=21.0, type=float)
    parser.add_option("--deltadecay", help="Maximum relative error for delta mag", default=.9, type=float)

    parser.add_option("--start_time", help="start time.", default='2010-07-25T00:00:00')
    parser.add_option("--end_time", help="start time.", default='2029-08-25T00:00:00')

    parser.add_option("--thetadist", help="(sine, uniform)",default="uniform")
    parser.add_option("--thetai", help="Initial theta", default=0, type=float)
    parser.add_option("--thetaf", help="Final theta", default=90, type=float)

    parser.add_option("--mag", help="mag.", default=-16, type=float)
    parser.add_option("--dmag", help="dmag.", default=0.0, type=float)

    parser.add_option("--dmin", help="minimum allowed luminosity distance in Mpc", default=1.0, type=float)
    parser.add_option("--dmax", help="maximum allowed luminosity distance in Mpc", default=200.0, type=float)

    parser.add_option("-m", "--modelType", help="(Bulla, Bulladynwind, kasen, Tophat, afterglow)",default="Tophat")
    parser.add_option("--mdyn", help="Mass dynamical ejecta.", default=0.02, type=float)
    parser.add_option("--mwin", help="Mass disk wind.", default=0.13, type=float)

    parser.add_option("--NPoints", help="How many points for detection?", default=2, type=int)
    parser.add_option("--threshold", help="Required SNR", default=5., type=float)
    parser.add_option("--minseparation", help="Minimun data points time separation in days", default=0.0104, type=float)
    parser.add_option("--maxseparation", help="Maximum data points time separation in days", default=7., type=float)
    parser.add_option("--gal", help="galactic latitude limit in degrees", default=0., type=float)
    parser.add_option("--doBS",  action="store_true", default=False)

    parser.add_option("--E0", help="E0.", default=1e53, type=float)

    parser.add_option("--mej", help="mej", default=0.04, type=float)
    parser.add_option("--opening_angle", help="phi", default=30, type=int)
    parser.add_option("--temp", help="T", default=5000, type=int)

    parser.add_option("--doSimSurvey",  action="store_true", default=False)
    parser.add_option("--doFilter",  action="store_true", default=False)
    parser.add_option("--host",  action="store_true", default=False)

    parser.add_option("--doPlots",  action="store_true", default=False)
    parser.add_option("--useRate", action="store_true", default=False)
    parser.add_option("--sourcenoise", action="store_true", default=False)

    parser.add_option("--pickleFile", help="file to store lightcurves", default="sim.pkl", type=str)
    parser.add_option("--pickleFileFilter", help="file to store filtered lightcurves", default="sim_filter.pkl", type=str)

    parser.add_option("--doAllQuadrants",  action="store_true", default=False)

    parser.add_option("--seed", help="fix the random number generator seed", type=int, default=None)

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

def load_ztf_bands(bandpass_dir="./"):
    bands = {
        'ztfg' : 'ztfg_eff.txt',
        'ztfr' : 'ztfr_eff.txt',
        'ztfi' : 'ztfi_eff.txt',
    }

    for bandname in bands.keys() :
        fname = bands[bandname]
        b = np.loadtxt(os.path.join(bandpass_dir,fname))
        band = sncosmo.Bandpass(b[:,0], b[:,1], name=bandname)
        sncosmo.registry.register(band, force=True)

def clean_df(df, zp=30):

    df = df[df.status!=0]

    df['skynoise'] = [10 ** (-0.4 * (i - zp)) / 5. for i in df['limMag']]
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

def load_pointings():
    file = np.load("data/ZTF_survey/FullAlertQuery_20190408.npy")
    alert_jd = []
    alert_limMag = []
    alert_chid = []
    alert_ccdid = []
    alert_progid = []
    alert_fieldid = []
    alert_filterid = []

    for i in range(len(file)):
        alert_jd.append(file[i][0])
        alert_limMag.append(file[i][1])
        alert_chid.append(int(file[i][2]))
        alert_ccdid.append(int(compute_ccdid_from_channel(file[i][2])))
        alert_progid.append(file[i][3])
        alert_fieldid.append(int(file[i][4]))
        alert_filterid.append(int(file[i][5]))

    alert_data = {}
    alert_data["jd"] =alert_jd
    alert_data["limMag"] =alert_limMag
    alert_data["chid"] =alert_chid # dont need
    alert_data["ccdid"] =alert_ccdid # dont need
    alert_data["progid"] =alert_progid
    alert_data["fieldid"] =alert_fieldid
    alert_data["filterid"] =alert_filterid
    #need rcid
    #jdstart and end of ref
    #FWHM

    return pd.DataFrame(alert_data)

def load_fields_ccd(fields_file, ccd_file, mwebv=False, galactic=False, num_segs=64):
    fields_raw = np.genfromtxt(fields_file, comments='%')

    fields = {'field_id': np.array(fields_raw[:,0], dtype=int),
              'ra': fields_raw[:,1],
              'dec': fields_raw[:,2]}

    if mwebv:
        fields['mwebv'] = fields_raw[:,3]
    if galactic:
        fields['l'] = fields_raw[:,4]
        fields['b'] = fields_raw[:,5]


    ccd_corners = np.genfromtxt(ccd_file, skip_header=True)
    ccds = [ccd_corners[4*k:4*k+4, :2] for k in range(num_segs)]
    #ccds = []
    #for ccd_corner in ccd_corners:
    #    ccd_list = []
    #    for k in range(64):
    #        ccd_list.append(ccd_corner)
    #    ccd_list = np.array(ccd_list)[:,:2]
    #    ccds.append(ccd_list)
    #ccds = [ccd_corners[np.array([0,1,2,3])+4*k, :2] for k in range(64)]
    #ccds = [ccd_corners[np.array([0,1,2,3])+4*k, :2] for k in range(64)]

    return fields, ccds

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

def mattia_template(dynwind=False, dataDir=".",mej=0.04,mdyn=0.02, mwin=0.13, phi=30, temp=5000):
    if dynwind:
        l = dataDir+'nsns_nph1.0e+06_mejdyn'+'{:.3f}'.format(mdyn)+'_mejwind'+'{:.3f}'.format(mwin)+'_phi'+'{:.0f}'.format(phi)+'.txt'
    else:
        l = dataDir+'nph1.0e+06_mej'+'{:.3f}'.format(mej)+'_phi'+'{:.0f}'.format(phi)+'.txt'
    f = open(l)
    lines = f.readlines()

    nobs = float(lines[0])
    nwave = float(lines[1])
    line3 = (lines[2]).split(' ')
    ntime = float(line3[0])
    t_i = float(line3[1])
    t_f = float(line3[2])

    cos_theta = np.linspace(0, 1, nobs)  # 11 viewing angles
    phase = np.linspace(t_i, t_f, ntime)  # epochs

    file_ = np.genfromtxt(l, skip_header=3)

    wave = file_[0:int(nwave),0]
    flux = []
    for i in range(int(nobs)):
        flux.append(file_[i*int(nwave):i*int(nwave)+int(nwave),1:])
    flux = np.array(flux).T

    phase = np.linspace(t_i, t_f, len(flux.T[0][0]))  # epochs

    return phase, wave, cos_theta, flux


def kasen_(dataDir="."):
    bluename = dataDir+'knova_d1_n10_m0.025_vk0.30_Xlan1e-4.5.h5'
    bluefin    = h5py.File(bluename,'r')
    redname = dataDir+'knova_d1_n10_m0.040_vk0.10_Xlan1e-2.0.h5'
    redfin    = h5py.File(redname,'r')

    fin = bluefin # or redname, it is the same

    # array of time in seconds
    times = np.array(fin['time'])
    # covert time to days
    phase = times/3600.0/24.0

    # frequency in Hz
    nu    = np.array(fin['nu'],dtype='d')

    # specific luminosity (ergs/s/Hz)
    # this is a 2D array, Lnu[times][nu]
    Lnu_blue   = np.array(bluefin['Lnu'],dtype='d')
    Lnu_red    = np.array(redfin['Lnu'],dtype='d')
    Lnu = Lnu_blue + Lnu_red

    # if you want thing in Flambda (ergs/s/Angstrom)
    c    = 2.99e10
    lam  = c/nu*1e8
    Llam = Lnu*nu**2.0/c/1e8

    dl = 10 * u.pc
    Flam = Llam / (4*np.pi*dl.to(u.cm).value**2)

    wave = np.flipud(lam)
    flux = Flam[:,::-1]

    return phase, wave, flux

def easyint(x,y,xref):
    ir = (xref>=min(x))&(xref<=max(x))
    yint = interp.interp1d(x[np.argsort(x)],y[np.argsort(x)])(xref[ir])
    #yout = np.zeros(len(xref),dmodel=float)
    yout = np.zeros(len(xref),)
    yup = y[-1]
    ylow = y[0]
    yout[ir] = yint
    yout[xref<np.min(x)] = ylow
    yout[xref>np.max(x)] = yup
    return yout

def tophat_template(bandpass_dir, passband='ztfg', mag0=0.0, dmag=0.0):

    bands = {
        'ztfg' : 'ztfg_eff.txt',
        'ztfr' : 'ztfr_eff.txt',
        'ztfi' : 'ztfi_eff.txt',
    }

    nlambda, lambda_min, lambda_max = 100, 1000, 25000
    wave = np.linspace(lambda_min, lambda_max, nlambda)

    ntime, t_i, t_f = 30, 0., 15.25
    phase = np.linspace(t_i, t_f, ntime)  # epochs

    magdiff = mag0 + phase*dmag
    F_Lxlambda2 = 10**(-(magdiff+2.406)/2.5)

    waves, F_Lxlambda2s = np.meshgrid(wave, F_Lxlambda2)
    flux = F_Lxlambda2s/(waves)**2

    return phase, wave, flux

def afterglow_template(jetType = 0, specType = 0, ksiN = 1.0, dL = 3.09e19,
                       theta_v = 0.0, E0 = 1e53,
                       theta_c = np.pi/4.0, theta_w = np.pi/4.0,
                       n = 1, p = 2.5, epsilon_E = 1e-1, epsilon_B = 1e-2,
                       b = 0):

    from afterglowpy import fluxDensity

    nlambda, lambda_min, lambda_max = 100, 1000, 25000
    wave = np.linspace(lambda_min, lambda_max, nlambda)
    nu = 3e8/(wave*1e-10)

    ntime, t_i, t_f = 30, 0.25, 15.25
    phase = np.linspace(t_i, t_f, ntime)  # epochs

    Y = np.array([theta_v, E0, theta_c, theta_w, b, 0, 0, 0,
                  n, p, epsilon_E, epsilon_B, ksiN, dL])

    flux = []
    for t in phase:
        t = t*np.ones(nu.shape)
        mJy = fluxDensity(t, nu, jetType, specType, *Y)
        flux.append(mJy * 1e-3 * 2.99792458E-05 / (wave**2))

    return phase, wave, np.array(flux)

class AngularTimeSeriesSource(sncosmo.Source):
    """A single-component spectral time series model.
        The spectral flux density of this model is given by
        .. math::
        F(t, \lambda) = A \\times M(t, \lambda)
        where _M_ is the flux defined on a grid in phase and wavelength
        and _A_ (amplitude) is the single free parameter of the model. The
        amplitude _A_ is a simple unitless scaling factor applied to
        whatever flux values are used to initialize the
        ``TimeSeriesSource``. Therefore, the _A_ parameter has no
        intrinsic meaning. It can only be interpreted in conjunction with
        the model values. Thus, it is meaningless to compare the _A_
        parameter between two different ``TimeSeriesSource`` instances with
        different model data.
        Parameters
        ----------
        phase : `~numpy.ndarray`
        Phases in days.
        wave : `~numpy.ndarray`
        Wavelengths in Angstroms.
        cos_theta: `~numpy.ndarray`
        Cosine of
        flux : `~numpy.ndarray`
        Model spectral flux density in erg / s / cm^2 / Angstrom.
        Must have shape ``(num_phases, num_wave, num_cos_theta)``.
        zero_before : bool, optional
        If True, flux at phases before minimum phase will be zeroed. The
        default is False, in which case the flux at such phases will be equal
        to the flux at the minimum phase (``flux[0, :]`` in the input array).
        name : str, optional
        Name of the model. Default is `None`.
        version : str, optional
        Version of the model. Default is `None`.
        """

    _param_names = ['amplitude', 'theta']
    param_names_latex = ['A', r'\theta']

    def __init__(self, phase, wave, cos_theta, flux, zero_before=True, zero_after=True, name=None,
                 version=None):
        self.name = name
        self.version = version
        self._phase = phase
        self._wave = wave
        self._cos_theta = cos_theta
        self._flux_array = flux
        self._parameters = np.array([1., 0.])
        self._current_theta = 0.
        self._zero_before = zero_before
        self._zero_after = zero_after
        self._set_theta()

    def _set_theta(self):
        logflux_ = np.zeros(self._flux_array.shape[:2])
        for k in range(len(self._phase)):
            adding = 1e-10 # Here we add 1e-10 to avoid problems with null values
            f_tmp = Spline2d(self._wave, self._cos_theta, np.log(self._flux_array[k]+adding),
                             kx=1, ky=1)
            logflux_[k] = f_tmp(self._wave, np.cos(self._parameters[1]*np.pi/180)).T

        self._model_flux = Spline2d(self._phase, self._wave, logflux_, kx=1, ky=1)


        self._current_theta = self._parameters[1]

    def _flux(self, phase, wave):
        if self._current_theta != self._parameters[1]:
            self._set_theta()
        f = self._parameters[0] * (np.exp(self._model_flux(phase, wave)))
        if self._zero_before:
            mask = np.atleast_1d(phase) < self.minphase()
            f[mask, :] = 0.
        if self._zero_after:
            mask = np.atleast_1d(phase) > self.maxphase()
            f[mask, :] = 0.
        return f

def model_(phase=0, wave=0, cos_theta=0, flux=0, viewingangle_dependece=True):
    '''phase, wave, cos_theta and flux are 1D arrays with the same length'''
    if viewingangle_dependece:
        source = AngularTimeSeriesSource(phase, wave, cos_theta, flux)
    else:
        source = sncosmo.TimeSeriesSource(phase, wave, flux)
    dust = sncosmo.CCM89Dust()
    return sncosmo.Model(source=source,effects=[dust], effect_names=['MW'], effect_frames=['obs'])

def random_parameters_(redshifts, model,
                       r_v=2., ebv_rate=0.11,
                       **kwargs):
    # Amplitude
    amp = []
    for z in redshifts:
        amp.append(10**(-0.4*Planck15.distmod(z).value))

    if opts.thetadist=='sine':
        distribution = np.arccos(np.random.random(len(redshifts))) / np.pi * 180 # distribution of angles. Sine distribution of theta (viewing angle) ($\^{circ}$)
    elif opts.thetadist=='uniform':
        distribution = np.random.uniform(opts.thetai, opts.thetaf,size=len(redshifts)) # distribution of angles. Sine distribution of theta (viewing angle) ($\^{circ}$)
    if opts.host:
        return {
            'amplitude': np.array(amp),
            'theta': distribution,
            'hostr_v': r_v * np.ones(len(redshifts)),
            'hostebv': np.random.exponential(ebv_rate, len(redshifts))
            }
    else:
        return {
            'amplitude': np.array(amp),
            'theta': distribution
            }

def random_parameters_notheta_(redshifts, model,
                               r_v=2., ebv_rate=0.11,
                               **kwargs):
    # Amplitude
    amp = []
    for z in redshifts:
        amp.append(10**(-0.4*Planck15.distmod(z).value))

    if opts.host:
        return {
            'amplitude': np.array(amp),
            'hostr_v': r_v * np.ones(len(redshifts)),
            'hostebv': np.random.exponential(ebv_rate, len(redshifts))
            }
    else:
        return {
            'amplitude': np.array(amp)
            }

# =============================================================================
#
#                                    MAIN
#
# =============================================================================

warnings.filterwarnings("ignore")

if opts.seed is not None:
    np.random.seed(int(opts.seed))

#set up directory
inputDir = opts.inputDir
outputDir = opts.outputDir
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)

sfd98_dir = opts.sfdDir
survey_fields = os.path.join(inputDir,'ZTF_Fields.txt')
ccd_corners = os.path.join(inputDir,'ztf_rc_corners.txt')

start_time = Time(opts.start_time, format='isot', scale='utc')
end_time = Time(opts.end_time, format='isot', scale='utc')

'''if not os.path.isfile(pcklFile):
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
        obstable['ccdid']=np.floor(obstable["rcid"]/4.0).astype(int)
        names = ('field', 'rcid', 'scimaglim', 'fid', 'programid', 'status')
        renames = ('fieldid', 'chid', 'limMag', 'filterid', 'progid', 'status')
        obstable.rename_columns(names, renames)

    elif opts.observation_type == "obsfilecomma":
        obstable = Table.read(opts.observations, format='ascii.csv',
                              data_start=2, data_end=-1)
        if "rcid" in list(obstable.columns):
            obstable['ccdid']=np.floor(obstable["rcid"]/4.0).astype(int)
            names = ('field', 'rcid', 'scimaglim', 'fid', 'programid', 'status')
            renames = ('fieldid', 'chid', 'limMag', 'filterid', 'progid', 'status')
            obstable.rename_columns(names, renames)
        else:
            obstable['ccdid']=np.floor(obstable["chid"]/4.0).astype(int)
    else:
        print("observation type not understood")
        exit(0)

    obstable = obstable.to_pandas()
    obstable =obstable[obstable["jd"]>start_time.jd]
    obstable =obstable[obstable["jd"]<end_time.jd]

    obstable.to_pickle(pcklFile)'''

'''df = pd.read_pickle(pcklFile)
df = clean_df(df, zp=opts.zeropoint)
df = df.astype({'fieldid': 'int64'})'''

df = pd.read_csv(opts.observations)
df = clean_df(df, zp=opts.zeropoint)
df =df[df["jd"]>start_time.jd]
df =df[df["jd"]<end_time.jd]

if opts.doAllQuadrants:
    df_group = df.groupby('jd')
    rows_to_concat, limmags, ccdids_all = [], [], []
    for ii, (g, data) in enumerate(df_group):
        if np.mod(ii, 100) == 0:
            print('Process group %d/%d' % (ii, len(df_group)))
        df_sub = df_group.get_group(g)
        ccdids = df_sub["ccdid"].tolist()
        limmag_median = np.median(df_sub["limMag"].tolist())
        if len(ccdids) == 64: continue
        quadrantIDs = np.arange(64)
        missing_quadrants = np.setdiff1d(quadrantIDs, ccdids)
        for missing_quadrant in missing_quadrants:
            row = df_sub.iloc[0]
            limmags.append(limmag_median)
            ccdids_all.append(missing_quadrant)
            rows_to_concat.append(row)
    df_new = pd.concat(rows_to_concat, ignore_index=True, axis=1)
    df_new = df_new.T
    df_new["limMag"] = limmags
    df_new["chid"] = ccdids_all
    df = pd.concat([df, df_new], ignore_index=True, axis=0)
    df = df.astype({"fieldid": int})

if opts.doPlots:
    plotName = os.path.join(outputDir,'limmag_'+opts.pickleFile[:-4]+'.pdf')
    plt.figure()
    bands = list(set(df['filter']))
    for i in range(len(bands)):
        if bands[i]=='ztfg':
            color='blue'
        elif bands[i]=='ztfr':
            color='red'
        elif bands[i]=='ztfi':
            color='orange'
        else:
            color='black'
        mask = df['filter']==bands[i]
        plt.scatter(df.jd[mask]-np.min(df.jd),df.limMag[mask], marker='.',color=color, label=bands[i], alpha=.5)
    plt.gca().invert_yaxis()
    plt.legend(loc=0)
    plt.xlabel('Time [days]')
    plt.ylabel('Magnitude limit')
    plt.tight_layout()
    plt.savefig(plotName)
    plt.close()

if opts.doSimSurvey:
    # Load in the field and ccd corners, as well as filter files
    #sst.load_ztf_filters()
    fields, ccds = load_fields_ccd(survey_fields,ccd_corners)
    limMags = []
    lims = np.array(df['limMag'])
    for field_id in fields['field_id']:
        idx = np.where(df['fieldid'] == field_id)[0]
        if len(idx) == 0:
            limMags.append(0.0)
        else:
            limMags.append(np.median(lims[idx]))
    fields['limMags'] = np.array(limMags)

    load_ztf_bands(bandpass_dir=inputDir)

    zp = np.ones(len(df))*opts.zeropoint

    print('Loading survey')
    plan = simsurvey.SurveyPlan(time=df['jd'],
                                band=df['filter'],
                                obs_field=df['fieldid'].astype(int),
                                obs_ccd=df['chid'],
                                skynoise=df['skynoise'],
                                zp=zp,
                                comment=df['progid'],
                                fields={k: v for k, v in fields.items()
                                        if k in ['ra', 'dec', 'field_id']},
                                ccds=ccds)

    mjd_range = (plan.pointings['time'].min()-20, plan.pointings['time'].max()+20)

    distance_lower = Distance(opts.dmin * u.Mpc)
    distance_upper = Distance(opts.dmax * u.Mpc)
    z_min = distance_lower.z
    z_max = distance_upper.z

    if opts.modelType == "Bulla":
        mej = opts.mej
        phi = opts.opening_angle
        temp = opts.temp
        phase, wave, cos_theta, flux = mattia_template(dataDir=opts.possisDir,mej=mej, phi=phi, temp=temp)
        if phi==0 or phi==90:
            model = model_(phase=phase, wave=wave, cos_theta=cos_theta, flux=flux, viewingangle_dependece=False) # If viewingangle_dependece=False,it won't take into account cos_theta
            transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_notheta_)
        else:
            model = model_(phase=phase, wave=wave, cos_theta=cos_theta, flux=flux, viewingangle_dependece=True) # If viewingangle_dependece=False,it won't take into account cos_theta
            transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_)
    elif opts.modelType == "Bulladynwind":
        mdyn = opts.mdyn
        mwin = opts.mwin
        phi = opts.opening_angle
        temp = opts.temp
        phase, wave, cos_theta, flux = mattia_template(dynwind=True, dataDir=opts.possisDir,mdyn=mdyn, mwin=mwin, phi=phi, temp=temp)
        if phi==0 or phi==90:
            model = model_(phase=phase, wave=wave, cos_theta=cos_theta, flux=flux, viewingangle_dependece=False) # If viewingangle_dependece=False,it won't take into account cos_theta
            transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_notheta_)
        else:
            model = model_(phase=phase, wave=wave, cos_theta=cos_theta, flux=flux, viewingangle_dependece=True) # If viewingangle_dependece=False,it won't take into account cos_theta
            transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_)
    elif opts.modelType == "kasen":
        phase, wave, flux = kasen_(dataDir=opts.kasenDir)
        model = model_(phase=phase, wave=wave, flux=flux, viewingangle_dependece=False) # If viewingangle_dependece=False,it won't take into account cos_theta
        transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_notheta_)
    elif opts.modelType == "afterglow":
        E0 = opts.E0
        n = 1.0
        theta_v = 0.0
        theta_c = np.pi/4.0
        theta_w = np.pi/4.0
        p = 2.5
        epsilon_B = 1e-2
        epsilon_E = 1e-1

        jetType = 0
        specType = 0
        ksiN = 1.0
        dL = 3.09e19

        phase, wave, flux = afterglow_template(jetType = jetType,
                                               specType = specType,
                                               ksiN = ksiN,
                                               dL = dL,
                                               theta_v = theta_v, E0 = E0,
                                               theta_c = theta_c,
                                               theta_w = theta_w,
                                               n = n, p = p,
                                               epsilon_E = epsilon_E,
                                               epsilon_B = epsilon_B)

        model = model_(phase=phase, wave=wave, flux=flux, viewingangle_dependece=False) # If viewingangle_dependece=False,it won't take into account cos_theta
        transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_notheta_)
    elif opts.modelType == "Tophat":
        phase, wave, flux = tophat_template(bandpass_dir=inputDir, mag0=opts.mag, dmag=opts.dmag)
        model = model_(phase=phase, wave=wave, flux=flux, viewingangle_dependece=False) # If viewingangle_dependece=False,it won't take into account cos_theta
        transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_notheta_)


    ntransient = opts.ntransient # if ntransient==None, it will simulate LCs following the rate
    if opts.useRate==True:
        ntransient = None
    dec_range = (-30,90)
    ra_range = (0, 360)

    # path to galactic extinction data
    sfd98_dir = opts.sfdDir

    rate = opts.rate # Even if we fix a number of events we need to give a rate.
    tr = simsurvey.get_transient_generator([z_min, z_max],
                                          ntransient=ntransient,
                                          ratefunc=lambda z: rate,
                                          dec_range=dec_range,
                                          ra_range=ra_range,
                                          mjd_range=(mjd_range[0],
                                                     mjd_range[1]),
                                          transientprop=transientprop,
                                          sfd98_dir=sfd98_dir
                                          )

    if opts.doPlots:
        plotName = os.path.join(outputDir,'field_hist_'+opts.pickleFile[:-4]+'.pdf')
        plt.figure()
        _,_,_ = plt.hist(df["fieldid"])
        plt.xlabel('fieldid')
        plt.ylabel('N')
        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()

    survey = simsurvey.SimulSurvey(generator=tr, plan=plan, sourcenoise=opts.sourcenoise,
                                   n_det=opts.NPoints, threshold=opts.threshold) # By default 2 detections required (n_det=2)
    if survey._properties['instruments'] is None:
        print('There is no hope... setting results to 0.')
        filename = os.path.join(outputDir, 'KNe.dat')
        fid = open(filename, 'w')
        fid.write('%.5e %.5e %.5e' % (0.0, 0.0, 0.0))
        fid.close()

        exit(0)

    lcs = survey.get_lightcurves(progress_bar=True, notebook=True)

    if lcs.meta_full == None:
        print("No lightcurves... returning.")
        exit(0)

    pcklfilename = opts.pickleFile

    pcklFile = os.path.join(outputDir,pcklfilename)
    pickle.dump(lcs,open(pcklFile, "wb" ) )


    def program1(lc):
        return lc[lc['comment']==1]
    def program2(lc):
        return lc[lc['comment']==2]
    def program3(lc):
        return lc[lc['comment']==3]

    print()
    print('From simsurvey:', len(lcs.lcs))
    print()
    try:
        print('In programid==1: ', len(lcs.filter(filterfunc=program1).filter(n_det=opts.NPoints, threshold=opts.threshold).lcs))
    except:
        pass
    try:
        print('In programid==2: ', len(lcs.filter(filterfunc=program2).filter(n_det=opts.NPoints, threshold=opts.threshold).lcs))
    except:
        pass
    try:
        print('In programid==3: ', len(lcs.filter(filterfunc=program3).filter(n_det=opts.NPoints, threshold=opts.threshold).lcs))
    except:
        pass
    print()
    ntransient = len(lcs.meta_notobserved)+len(lcs.meta_full)
    print('Number of created kNe:', ntransient)
    print('Number of created kNe falling in the covered area:', len(lcs.meta_full['z']))

    if opts.doPlots:
        # plot N vs redshift
        plotName = os.path.join(outputDir,'redshift_'+opts.pickleFile[:-4]+'.pdf')
        plt.figure()
        range_bins=(0,lcs.meta_full['z'].max())
        plt.hist(lcs.meta_full['z'],range=range_bins, label='created')
        plt.legend(loc=0)
        plt.xlabel('Redshift')
        plt.ylabel('N')
        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()

        # N vs detection phase
        plotName = os.path.join(outputDir,'detection_phase_'+opts.pickleFile[:-4]+'.pdf')
        plt.figure()
        plt.hist(lcs.stats['p_det'])
        plt.xlabel('Detection date (d)')
        plt.ylabel('N')
        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()

        # N vs last 5 sigma detection
        plotName = os.path.join(outputDir,'last5s_'+opts.pickleFile[:-4]+'.pdf')
        plt.figure()
        plt.hist(lcs.stats['p_last'])
        plt.xlabel('Last date (d)')
        plt.ylabel('N')
        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()

        cmap = cm.summer

        plotName = os.path.join(outputDir,'radec_'+opts.pickleFile[:-4]+'.pdf')
        plt.figure()
        ax = plt.axes(
            [0.05, 0.05, 0.9, 0.9],
            projection='geo degrees mollweide'
            )
        ax.grid()

        ax.scatter(lcs.meta_notobserved['ra'], lcs.meta_notobserved['dec'],  transform=ax.get_transform('world'), marker='*',color='grey', label='not_observed', alpha=0.7)
        ax.scatter(lcs.meta_rejected['ra'], lcs.meta_rejected['dec'],  transform=ax.get_transform('world'), marker='*',color='blue', label='rejected', alpha=0.7)
        ax.scatter(lcs.meta['ra'], lcs.meta['dec'],  transform=ax.get_transform('world'), marker='*', color='red', label='detected', alpha=0.7)

        ax.legend(loc=0)
        ax.set_ylabel('DEC (deg)')
        ax.set_xlabel('RA (deg)')

        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()

if opts.doFilter:

    pcklfilename = opts.pickleFile

    fields, ccds = load_fields_ccd(survey_fields,ccd_corners)
    pcklFile = os.path.join(outputDir,pcklfilename)
    lcs = pickle.load(open(pcklFile, "rb" ) )

    def bright_stars_loss(lc, bst):
        f, c, idx, w = [],[],[],[]
        for i in range(len(lc.lcs)):
            try:
                loss = float(bst[(bst.field==np.min(lc.lcs[i]['field']))&(bst.ccd==np.min(lc.lcs[i]['ccd']))].loss.values)
            except:
                loss = 0
            w.append(loss)
            f.append(np.min(lc.lcs[i]['field']))
            c.append(np.min(lc.lcs[i]['ccd']))
            idx.append(lc.meta['idx_orig'][i])
        w = np.array(w)
        f = np.array(f)
        c = np.array(c)
        ndet = 0
        for i in np.unique(np.array([f,c]).T, axis=0):
            try:
                loss = float(bst[(bst.field==i[0])&(bst.ccd==i[1])].loss.values)
            except:
                loss = 0
            ndet += len(f[(f==i[0])&(c==i[1])]) * (1 - loss)
        ndet = round(ndet)
        idx_orig = []
        idx_ = np.array(idx)
        w_ = np.array(w)
        while len(np.unique(idx_orig))<(len(lc.lcs)-ndet):
            size = (len(lc.lcs)-ndet) - len(np.unique(idx_orig))
            n = random.choices(idx_, weights=w_, k=size)
            n = np.unique(n)
            wh = np.in1d(idx_, n)
            idx_ = idx_[~wh]
            w_ = w_[~wh]
            idx_orig = idx_orig + list(n)
        return ndet, idx_orig

    if opts.doBS:
        bst = pd.read_csv(opts.bsfile, names = ['field', 'ccd', 'loss'])# Z distritbution: it should be all for 20 mag
        ndet, idx_orig = bright_stars_loss(lcs, bst)
    else:
        idx_orig = []

    def dm(m1, m2, m1err, m2err):
        dm = m2 - m1
        dmerr = np.sqrt(m1err**2 + m2err**2)
        return dmerr/dm

    def filterfunc(lc):
        lc_ = lc[lc['time']>=np.argmax(lc['time'])]
        time = lc_['time']
        mag  = 30 - 2.5 * np.log10(lc_['flux'])
        magerr = 2.5 / np.log(10) * lc_['fluxerr']/lc_['flux']
        bands = lc_['band']

        go = False
        for b in np.unique(bands):
            mag_    = mag[bands==b]
            magerr_ = magerr[bands==b]
            for i in range(len(mag_)):
                if dm(mag_[0], mag_[i], magerr_[0], magerr_[i])<opts.deltadecay:
                    go = True

        time = lc[lc['flux']/lc['fluxerr']>=opts.threshold]['time']
        dt_days = (time[-1]-time[0])
        separation_min = opts.minseparation #days
        separation_max = opts.maxseparation #days
        mask_dt = dt_days>=separation_min and dt_days<=separation_max

        c = SkyCoord(ra=lc.meta['ra'], dec=lc.meta['dec'], unit=(u.degree, u.degree))
        gal_lat_limit = opts.gal    # galactic latitude filter (+- in degrees)
        mask_gal = c.galactic.b.degree>=gal_lat_limit or c.galactic.b.degree <= -gal_lat_limit

        if mask_dt and mask_gal and (lc.meta['idx_orig'] not in idx_orig) and go and np.nanmin(mag)<=opts.peakmag:
            return lc


    lcs_filter = lcs.filter(filterfunc=filterfunc, n_det=opts.NPoints, threshold=opts.threshold)

    pcklfilenamefilter = opts.pickleFileFilter
    pcklFile_filter = os.path.join(outputDir,pcklfilenamefilter)
    pickle.dump(lcs_filter,open(pcklFile_filter, "wb" ) )

    if lcs_filter.lcs is None:
        print('There is no hope... setting results to 0.')
        filename = os.path.join(outputDir, 'KNe.dat')
        fid = open(filename, 'w')
        fid.write('%.5e %.5e %.5e' % (0.0, 0.0, 0.0))
        fid.close()

        exit(0)

    # Search for peak and count filtered kNe
    peak_ = []
    kN_n=0
    kN_1=0
    kN_2=0
    for i in range(len(lcs_filter.lcs)):
        lc = lcs_filter.lcs[i]
        f = lc[lc['flux']/lc['fluxerr']>=opts.threshold]['flux']
        if len(f)>=1:
            kN_1 += 1
        if len(f)>=2:
            kN_2 += 1
        if len(f)>=opts.NPoints:
            kN_n += 1
        if len(f) and (f[0]<np.max(f) and np.max(f)<f[-1]):
            peak = True
        else:
            peak = False
        peak_.append(peak)

    def program1(lc):
        return lc[lc['comment']==1]
    def program2(lc):
        return lc[lc['comment']==2]
    def program3(lc):
        return lc[lc['comment']==3]

    print()
    print('Before filtering kNe:', len(lcs.lcs))
    print('Filtered kNe:', len(lcs_filter.lcs))
    print('kNe showing a maximum:', sum(peak_))
    print()
    try:
        print('In programid==1: ', len(lcs_filter.filter(filterfunc=program1).filter(filterfunc=filterfunc, n_det=opts.NPoints, threshold=opts.threshold).lcs))
    except:
        pass

    try:
        print('In programid==2: ', len(lcs_filter.filter(filterfunc=program2).filter(filterfunc=filterfunc, n_det=opts.NPoints, threshold=opts.threshold).lcs))
    except:
        pass

    try:
        print('In programid==3: ', len(lcs_filter.filter(filterfunc=program3).filter(filterfunc=filterfunc, n_det=opts.NPoints, threshold=opts.threshold).lcs))
    except:
        pass
    print()
    ntransient = len(lcs.meta_notobserved)+len(lcs.meta_full)
    print('Number of created kNe:', ntransient)
    print('Number of created kNe falling in the covered area:', len(lcs.meta_full['z']))

    '''print('Number of detected over all created:', len(lcs.lcs)/ntransient, ' (efficiency?)')
    print('Number of filtered over all created:', len(lcs_filter.lcs)/ntransient, ' (efficiency?)')

    print('Number of detected over created in the observed fields:', len(lcs.lcs)/len(lcs.meta_full['z']), ' (efficiency)')
    print('Number of filtered over created in the observed fields:', len(lcs_filter.lcs)/len(lcs.meta_full['z']), ' (efficiency)')'''


    filename = os.path.join(outputDir, 'KNe.dat')
    fid = open(filename, 'w')
    fid.write('%.5e %.5e %.5e' % (kN_1/ntransient,kN_2/ntransient,kN_n/ntransient))
    fid.close()

    if opts.doPlots:
        jd_, band_, f_, ferr_, z_, p_det_, t0_, tlast5s_, ra_, dec_, th_ = [],[],[],[],[],[],[],[],[],[],[]
        for i in range(len(lcs_filter.lcs)):
            lc = lcs_filter.lcs[i]
            jd_.append(lc['time'])
            band_.append(lc['band'])
            f_.append(lc['flux'])
            ferr_.append(lc['fluxerr'])
            z_.append(lcs_filter.meta['z'][i])
            p_det_.append(lcs_filter.stats['p_det'][i])
            t0_.append(lcs_filter.meta['t0'][i])
            tlast5s_.append(lc['time'][lc['flux']/lc['fluxerr']>=opts.threshold][-1])
            ra_.append(lcs_filter.meta['ra'][i])
            dec_.append(lcs_filter.meta['dec'][i])
            if (opts.modelType=='Bulla' or opts.modelType=='Bulladynwind') and phi!=0 and phi!=90:
                th_.append(lcs_filter.meta['theta'][i])

        x_max = 0
        mag_limit = 22
        y_max = mag_limit
        plotName = os.path.join(outputDir,'lcs.pdf')
        plt.figure()
        filters = ['ztfg', 'ztfr', 'ztfi']
        colors = ['green', 'red', 'orange']
        for filt, color in zip(filters,colors):
            for i in range(len(jd_)):
                mask = band_[i] == filt
                x = jd_[i][mask]-t0_[i]
                y = opts.zeropoint-2.5*np.log10(f_[i][mask])
                yerr = ferr_[i][mask] / f_[i][mask]
                if len(x)>0 and np.max(x)>x_max:
                    x_max = np.max(x)
                if len(y)>0 and np.max(y)<y_max:
                    y_max = np.max(y)
                if i == 0:
                    plt.errorbar(x, y, yerr=yerr, marker='o',color=color, label=filt, alpha=0.5)
                else:
                    plt.errorbar(x, y, yerr=yerr, marker='o',color=color, alpha=0.5)
        plt.xlabel('$t-t_{merger}$')
        plt.ylabel('Magnitude')
        # Add LC fit here if needed
        plt.gca().invert_yaxis()
        plt.xlim(-.5,x_max+.5)
        plt.ylim(mag_limit+0.5,y_max-.5)
        plt.legend(loc=0)
        plt.tight_layout()
        plt.savefig(plotName)
        plt.show()
        plt.close()

        # plot N vs redshift
        plotName = os.path.join(outputDir,'redshift.pdf')
        plt.figure()
        range_bins=(0,lcs.meta_full['z'].max())
        plt.hist(lcs.meta_full['z'],range=range_bins, label='created')
        plt.hist(z_,range=range_bins, label='filtered_kNe')
        plt.legend(loc=0)
        plt.xlabel('Redshift')
        plt.ylabel('N')
        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()

        if lcs.meta is not None:
            # plot jd_detection vs survey time
            plotName = os.path.join(outputDir,'detection_date.pdf')
            jd_bins=np.linspace(0,end_time.jd-start_time.jd,round(end_time.jd-start_time.jd))
            plt.hist(lcs.stats['p_det']+lcs.meta['t0']-start_time.jd,bins=jd_bins)
            plt.hist(np.array(p_det_)+np.array(t0_)-start_time.jd,bins=jd_bins, label='filtered_kNe')
            plt.legend(loc=0)
            plt.xlabel('Detection date (d)')
            plt.ylabel('N')
            plt.tight_layout()
            plt.savefig(plotName)
            plt.close()

        # N vs detection phase
        plotName = os.path.join(outputDir,'detection_phase.pdf')
        plt.figure()
        plt.hist(lcs.stats['p_det'])
        plt.xlabel('Detection date (d)')
        plt.ylabel('N')
        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()

        # N vs last 5 sigma detection
        plotName = os.path.join(outputDir,'last5s.pdf')
        plt.figure()
        plt.hist(np.array(tlast5s_) - np.array(t0_))
        plt.xlabel('Last 5sigma detection (d)')
        plt.ylabel('N')
        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()

        cmap = cm.summer

        plotName = os.path.join(outputDir,'radec.pdf')
        plt.figure()
        idx = np.where(fields['limMags'] > 0)[0]
        sc = plt.scatter(fields['ra'][idx], fields['dec'][idx], marker='s', c=fields['limMags'][idx], alpha=.5, cmap=cmap)
        idx = np.where(fields['limMags'] == 0)[0]
        plt.scatter(fields['ra'][idx], fields['dec'][idx], marker='s', color='orange', alpha=.5, label='Pointings')
        plt.scatter(tr.ra, tr.dec, label='created', alpha=.5, color='blue')
        if lcs.meta is not None:
            plt.scatter(lcs.meta['ra'], lcs.meta['dec'],marker='*', label='detected', alpha=.5, color='green')
        plt.scatter(ra_, dec_, marker='*', label='filtered', alpha=.5, color='red')
        cbar = plt.colorbar(sc)
        cbar.set_label('Limiting Magnitudes')
        plt.legend(loc=0)
        plt.xlabel('RA')
        plt.ylabel('DEC')
        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()

        if (opts.modelType=='Bulla' or opts.modelType=='Bulladynwind') and phi!=0 and phi!=90:
            # 2D histogram of detections per redshift and viewing angle.
            plotName = os.path.join(outputDir,'2dhist.pdf')
            plt.figure()
            binsz = 20
            binsth = 9
            zmax = lcs.meta_full['z'].max()
            H, yedges, xedges, patches = plt.hist2d(z_, th_,bins=[binsz, binsth], range=[[0,zmax],[0,90]])
            plt.imshow(H,aspect='auto',interpolation='none', origin='low', extent=[yedges.min(), yedges.max(), xedges.min(), xedges.max()])
            plt.colorbar(label='N')
            plt.xlabel('z')
            plt.ylabel('viewing angle ($^\circ$)')
            plt.savefig(plotName)
            plt.close()

            # Hist viewing angle
            plotName = os.path.join(outputDir,'viewingangle.pdf')
            plt.figure()
            range_bins=(0,lcs.meta_full['theta'].max())
            plt.hist(lcs.meta_full['theta'],range=range_bins, label='created')
            plt.hist(th_,range=range_bins, label='filtered_kNe')
            plt.legend(loc=0)
            plt.xlabel('viewing angle ($^\circ$)')
            plt.ylabel('N')
            plt.tight_layout()
            plt.savefig(plotName)
            plt.close()
