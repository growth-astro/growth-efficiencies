#!/usr/bin/env python
from __future__ import division
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

#import ligo.skymap.plot

from astropy.coordinates import SkyCoord
from astropy.io import ascii
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

#matplotlib.rcParams['font.family']="sans-serif"
#matplotlib.rcParams['font.sans-serif']='Comic Sans MS'

#font = {'family' : 'normal',
 #       'weight' : 'normal',
  #      'size'   : 16}
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
#matplotlib.rc('font', **font)

#font = {'family' : 'normal',
 #       'weight' : 'normal',
  #      'size'   : 16}
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
#matplotlib.rc('font', **font)


from matplotlib.lines import Line2D
import matplotlib.cm as cm

from gwemopt.chipgaps import ztf
import h5py

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
    parser.add_option("-x", "--sfdDir", help="sfd directory",default="../input/sfd98/")
    parser.add_option("-p", "--possisDir", help="possis directory",default="../data/possis/")
    parser.add_option("-a", "--afterglowDir", help="aftterglow directory",default="../data/possis/")
    parser.add_option("-b", "--bsfile", help="Bright stars loss file",default="../input/bright_star_loss.csv")
    parser.add_option("-n", "--ntransient", help="int number", default=10000, type=int)
    parser.add_option("-r", "--rate", help="kilonovae rate", default=5e-7, type=float)
    parser.add_option("-z", "--zeropoint", help="zeropoint", default=30, type=float)
    parser.add_option("--doRotate",  action="store_true", default=False)
    parser.add_option("--theta", help="theta rotation.", default=0.0, type=float)
    parser.add_option("--phi", help="phi rotation.", default=0.0, type=float)

    parser.add_option("--start_time", help="start time.", default='2019-04-25T08:18:25')
    parser.add_option("-d", "--dt", help="length of time to check.", default=3.0, type=float)

    parser.add_option("--thetadist", help="(sine, uniform)",default="uniform")
    parser.add_option("--thetai", help="Initial theta", default=0, type=float)
    parser.add_option("--thetaf", help="Final theta", default=90, type=float)

    parser.add_option("--mag", help="mag.", default=-16, type=float)
    parser.add_option("--dmag", help="dmag.", default=0.0, type=float)

    parser.add_option("--dmin", help="minimum allowed luminosity distance in Mpc", default=0.1, type=float)
    parser.add_option("--dmax", help="maximum allowed luminosity distance in Mpc", default=1e10, type=float)

    parser.add_option("--doColor", help="Tophat template with default color evolution for GW170817", action='store_true', default=False)
    parser.add_option("--rmag_corr", help="r-band peak correction (compared to g=-16.92)", default=0.07, type=float)
    parser.add_option("--imag_corr", help="i-band peak correction (compared to g=-16.92)", default=0.47, type=float)
    parser.add_option("--diff_gi", help="reddening rate correction in i-band, from the g-band decay rate, (mags/day)", default=-0.83)
    parser.add_option("--diff_gr", help="reddening rate correction in r-band, from the g-band decay rate, (mags/day)", default=-0.49)

    parser.add_option("-m", "--modelType", help="(Bulla, Bulladynwind, BulladynwindNSBH, Tophat, afterglow, kasen)",default="Tophat")

    parser.add_option("--mdyn", help="Mass dynamical ejecta.", default=0.02, type=float)
    parser.add_option("--mwin", help="Mass disk wind.", default=0.09, type=float)

    parser.add_option("--NPoints", help="How many points for detection?", default=2, type=int)
    parser.add_option("--threshold", help="Required SNR", default=5., type=float)
    parser.add_option("--separation", help="Data points time separation in minutes", default=15., type=float)
    parser.add_option("--gal", help="galactic latitude limit in degrees", default=10., type=float)
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
    parser.add_option("--doField", help="report KN efficiency in a single field after full simulation", action="store_true", default=False)
    parser.add_option("--doSingleField", help="only simulate KNe in a single field", action="store_true", default=False)
    parser.add_option("--field", help="field ID to consider", default=800)

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
    df['skynoise'] = [10 ** (-0.4 * (i - zp)) / 5. for i in df['limMag']]  # use scimaglim for deeper limmags
    df["filter"]  = [rename_filter(i) for i in df['filterid']] #FIXME ONCE WE HAVE FID DATA
    print('Survey pointings for all programs:')
    print(len(df))
    print('Survey pointings for MSIP programs:')
    print(len(df[df["progid"]==1]))
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

def mattia_template(dynwind=False, dataDir=".",mej=0.04,mdyn=0.02, mwin=0.13, phi=30, temp=5000, nsns=False, nsbh=False):
    if dynwind:
        if nsns:
            mtype = 'nsns'
        elif nsbh:
            mtype = 'nsbh'
        l = dataDir+mtype+'_nph1.0e+06_mejdyn'+'{:.3f}'.format(mdyn)+'_mejwind'+'{:.3f}'.format(mwin)+'_phi'+'{:.0f}'.format(phi)+'.txt'
    else:
        l = dataDir+'nph1.0e+06_mej'+'{:.3f}'.format(mej)+'_phi'+'{:.0f}'.format(phi)+'.txt'
    f = open(l)
    lines = f.readlines()

    nobs = int(lines[0])
    nwave = float(lines[1])
    line3 = (lines[2]).split(' ')
    ntime = int(line3[0])
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
    bluename = dataDir+'kilonova_models/knova_d1_n10_m0.025_vk0.30_Xlan1e-4.5.h5'
    bluefin    = h5py.File(bluename,'r')
    redname = dataDir+'kilonova_models/knova_d1_n10_m0.040_vk0.10_Xlan1e-2.0.h5'
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

def tophatgw(bandpass_dir, gm=opts.mag, g_dm = opts.dmag, rmag_corr =0.07, imag_corr=0.47, diff_gr=-0.49, diff_gi=-0.83):
    """A model with mag_corr and diff set to zero will product a kilonova with no color or color-evolution."""
    Wave = [[], [], []]
    Flux = [[], [], []]
    bands = {
        'ztfg' : 'ztfg_eff.txt',
        'ztfr' : 'ztfr_eff.txt',
        'ztfi' : 'ztfi_eff.txt'
    }
    bounds = {
        'ztfg': [3676, 5496],
        'ztfr': [5497, 6867],
        'ztfi': [6868, 9022]
    }
    for i, filt in enumerate(bands.keys()):
        f = ascii.read(bandpass_dir+bands[filt])
        mask = (np.array(f["Wavelength (A)"]) > bounds[filt][0]) & (np.array(f["Wavelength (A)"]) < bounds[filt][1])
        Wave[i] = f["Wavelength (A)"][mask]
        ntime, t_i, t_f = 30, 0., 15.25
        phase = np.linspace(t_i, t_f, ntime)  # epochs
        magdiff = gm + phase*g_dm
        if filt == 'ztfg':
            F_Lxlambda2  = 10**(-(magdiff+2.406)/2.5)
        elif filt == 'ztfr': 
            magdiff2 = magdiff + diff_gr*phase + rmag_corr
            F_Lxlambda2  = 10**(-(magdiff2+2.406)/2.5)
        elif filt == 'ztfi': 
            magdiff2 = magdiff + diff_gi*phase + imag_corr
            F_Lxlambda2  = 10**(-(magdiff2+2.406)/2.5)
        waves, F_Lxlambda2s = np.meshgrid(Wave[i], F_Lxlambda2)
        Flux[i] = F_Lxlambda2s/(waves)**2
    Waveb = np.arange(1000, bounds['ztfg'][0])
    Wavea = np.arange(bounds['ztfi'][1]+1, 25000)
    wave = np.concatenate((Waveb, Wave[0], Wave[1], Wave[2], Wavea))
    waves, Fluxb = np.meshgrid(Waveb, phase*0)
    waves, Fluxa = np.meshgrid(Wavea, phase*0)
    flux = np.concatenate((Fluxb, Flux[0], Flux[1], Flux[2], Fluxa), axis=1)
    return phase, wave, flux


def afterglow_template(jetType = 0, specType = 0, ksiN = 1.0, dL = 3.09e19,
                       theta_v = 0.0, E0 = 1e53,
                       theta_c = np.pi/4.0, theta_w = np.pi/4.0,
                       n = 3.7e-3, p = 2.43, epsilon_E = 1e-1, epsilon_B = 1e-2,
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
        mJy = fluxDensity(t*86400.0, nu, jetType, specType, *Y)
        flux.append(mJy * 1e-3 * 2.99792458E-05 / (wave**2))

    return phase, wave, np.cos(theta_v), np.array(flux)

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

class TimeSeriesSource_theta(sncosmo.Source):

    _param_names = ['amplitude', 'theta']
    param_names_latex = ['A', r'\theta']

    def __init__(self, phase, wave, cos_theta, flux, zero_before=False, name=None,
                 version=None):
        self.name = name
        self.version = version
        self._phase = phase
        self._wave = wave
        self._parameters = np.array([1.])
        self._model_flux = Spline2d(phase, wave, flux, kx=2, ky=2)
        self._zero_before = zero_before
        self._cos_theta = cos_theta

    def _flux(self, phase, wave):
        f = self._parameters[0] * self._model_flux(phase, wave)
        if self._zero_before:
            mask = np.atleast_1d(phase) < self.minphase()
            f[mask, :] = 0.
        return f


def model_(phase=0, wave=0, cos_theta=0, flux=0, viewingangle_dependece=True, afterglow=False):
    '''phase, wave, cos_theta and flux are 1D arrays with the same length'''
    if viewingangle_dependece and afterglow==False:
        source = AngularTimeSeriesSource(phase, wave, cos_theta, flux)
    elif viewingangle_dependece and afterglow:
        source = TimeSeriesSource_theta(phase, wave, cos_theta, flux)
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

client = pyvo.dal.TAPService('https://irsa.ipac.caltech.edu/TAP')

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
        if opts.observations[-4:]=='.csv':
            obstable = pd.read_csv(opts.observations)
            obstable['ccdid']=np.floor(obstable['rcid']/4).astype(int)
        else:
            obstable = Table.read(opts.observations, format='ascii.fixed_width',
                                data_start=2, data_end=-1)
            obstable['ccdid']=np.floor(obstable['rcid']/4).astype(int)
            names = ('field', 'rcid', 'scimaglim', 'fid', 'programid', 'status')
            renames = ('fieldid', 'chid', 'limMag', 'filterid', 'progid', 'status')
            obstable.rename_columns(names, renames)
    else:
        print("observation type not understood")
        exit(0)
    try:
        obstable = obstable.to_pandas()
        print(obstable.columns)
        obstable =obstable[obstable["jd"]>start_time.jd]
        obstable =obstable[obstable["jd"]<end_time.jd]
    except AttributeError:
        obstable =obstable[obstable["jd"]>start_time.jd]
        obstable =obstable[obstable["jd"]<end_time.jd]

    obstable.to_pickle(pcklFile)

df = pd.read_pickle(pcklFile)
df = clean_df(df, zp=opts.zeropoint)

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
    '''plotName = os.path.join(outputDir,'limmag.pdf')
    plt.figure()
    plt.scatter(df.jd-np.min(df.jd),df.limMag)
    plt.gca().invert_yaxis()
    plt.xlabel('Time [days]')
    plt.ylabel('Magnitude limit')
    plt.tight_layout()
    plt.savefig(plotName)
    plt.close()'''

    plotName = os.path.join(outputDir,'limmag_b.pdf')
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

if opts.modelType=='afterglow':
    skymap, metadata = fits.read_sky_map(opts.skymap, nest=False, distances=False)
    map_struct = {}
    map_struct["prob"] = skymap
    #map_struct["distmu"] = np.ones(len(skymap))*250
    #map_struct["distsigma"] = np.ones(len(skymap))*1.5e5
    nside = hp.get_nside(skymap)
else:
    skymap, metadata = fits.read_sky_map(opts.skymap, nest=False, distances=True)
    map_struct = {}
    map_struct["prob"] = skymap[0]
    map_struct["distmu"] = skymap[1]
    map_struct["distsigma"] = skymap[2]
    map_struct["distnorm"] = skymap[3]
    nside = hp.get_nside(skymap[0])

npix = hp.nside2npix(nside)
if opts.doRotate:
    for key in map_struct.keys():
        map_struct[key] = rotate_map(map_struct[key], np.deg2rad(opts.theta), np.deg2rad(opts.phi))

sort_idx = np.argsort(map_struct["prob"])[::-1]
csm = np.empty(len(map_struct["prob"]))
csm[sort_idx] = np.cumsum(map_struct["prob"][sort_idx])
#map_struct["prob"][csm>0.9] = 0.0

fields = np.genfromtxt(survey_fields, comments='%', usecols=range(3), names=['field_id', 'ra', 'dec'])
rows = np.vstack((fields['field_id'], fields['ra'], fields['dec'])).T
fields = Table(rows=rows, names=('field_id', 'ra', 'dec'))

quadrant_coords = ztf.get_ztf_quadrants()
skyoffset_frames = coordinates.SkyCoord(fields['ra'], fields['dec'], unit=u.deg).skyoffset_frame()

quadrant_coords_icrs = coordinates.SkyCoord(
    *np.tile(quadrant_coords[:, np.newaxis, ...], (len(fields), 1, 1)), unit=u.deg,
    frame=skyoffset_frames[:, np.newaxis, np.newaxis]).transform_to(coordinates.ICRS)
quadrant_xyz = np.moveaxis(quadrant_coords_icrs.cartesian.xyz.value, 0, -1)

field_ids = np.unique(df.fieldid)
probs = {}
probstot_1, probstot_2 = [], []
ipixs = {}
ccds = []
for field_id, xyz in zip(fields['field_id'], quadrant_xyz):
    if not field_id in field_ids: continue
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
        probstot_1.append(np.sum(probs[field_id]))
    else:
        probstot_2.append(np.sum(probs[field_id]))


coverage = np.zeros(map_struct["prob"].shape)
for field_id in df.fieldid:
    ipix_list = np.array(ipixs[field_id], dtype="object")
    for ipix in ipix_list:
        coverage[ipix] = coverage[ipix] + 1
        if opts.doSingleField:
            if field_id == int(opts.field): 
                map_struct["prob"][ipix] = 1.0
            else:
                map_struct["prob"][ipix] = 1e-20
if opts.doSingleField: 
    map_struct["prob"] = map_struct["prob"]/np.sum(map_struct["prob"])
    df = df[df['fieldid'] == opts.field]

idx0 = np.where(coverage>0)[0]
idx1 = np.where(coverage>1)[0]
idx2 = np.where(coverage>2)[0]

filename = os.path.join(outputDir, 'prob.dat')
fid = open(filename, 'w')
fid.write('%.5f %.5f %.5f' % (np.sum(map_struct["prob"][idx0]),np.sum(map_struct["prob"][idx1]),np.sum(map_struct["prob"][idx2])))
fid.close()

print('Integrated probabilities')
print('1 exposure: %.5f' % np.sum(map_struct["prob"][idx0]))
print('2 exposures: %.5f' % np.sum(map_struct["prob"][idx1]))
print('3+ exposures: %.5f' % np.sum(map_struct["prob"][idx2]))

probarray = []
ipixarray = []
for ii, (fieldid, chid) in enumerate(zip(df.fieldid, df.chid.astype(int))):
    probarray.append(probs[fieldid][chid])
    ipixarray.append(ipixs[fieldid][chid])
df['prob'] = probarray
df['ipix'] = ipixarray

if opts.doSimSurvey:
    # Load in the field and ccd corners, as well as filter files
    #sst.load_ztf_filters()
    fields, ccds = load_fields_ccd(survey_fields,ccd_corners)
    limMags = []
#    limMagErr = []
    lims = np.array(df['limMag'])
    for field_id in fields['field_id']:
        idx = np.where(df['fieldid'] == field_id)[0]
        if len(idx) == 0:
            limMags.append(0.0)
        else:
            limMags.append(np.median(lims[idx]))
    fields['limMags'] = np.array(limMags)

#    for group in df.groupby('expid'):
#        limMagErr.append( np.std(group[:][1]['limMag'])*np.ones(len(group[:][1]['limMag'])) )

    load_ztf_bands(bandpass_dir=inputDir)

    zp = np.ones(len(df))*opts.zeropoint

    print('Loading survey')
    plan = simsurvey.SurveyPlan(time=df['jd'],
                                band=df['filter'],
                                obs_field=df['fieldid'].astype(int),
                                obs_ccd=df['chid'].astype(int),
                                skynoise=df['skynoise'],
                                zp=zp,
#                                limmag_err = np.concatenate(limMagErr),
                                comment=df['progid'],
                                fields={k: v for k, v in fields.items()
                                        if k in ['ra', 'dec', 'field_id']},
                                ccds=ccds)
    #mjd_range = (plan.pointings['time'].min()-10, plan.pointings['time'].max()+10)
    mjd_range = (start_time.jd, start_time.jd)

    if opts.modelType == "afterglow":
        distance_lower = Distance(opts.dmin * u.Mpc)
        distance_upper = Distance(opts.dmax * u.Mpc)
        z_min = distance_lower.z
        z_max = distance_upper.z
    else:
        distmean, diststd = parameters_to_marginal_moments(map_struct["prob"],
                                                           map_struct["distmu"],
                                                           map_struct["distsigma"])

        distance = Distance(distmean * u.Mpc)
        try:
            distance_lower = Distance((distmean - 5*diststd) * u.Mpc)
        except ValueError:
            distance_lower = Distance(opts.dmin * u.Mpc)
        distance_upper = Distance((distmean + 5*diststd) * u.Mpc)
        if distance_upper.value>opts.dmax:
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
        phase, wave, cos_theta, flux = mattia_template(dynwind=True, dataDir=opts.possisDir,mdyn=mdyn, mwin=mwin, phi=phi, temp=temp, nsns=True)
        if phi==0 or phi==90:
            model = model_(phase=phase, wave=wave, cos_theta=cos_theta, flux=flux, viewingangle_dependece=False) # If viewingangle_dependece=False,it won't take into account cos_theta
            transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_notheta_)
        else:
            model = model_(phase=phase, wave=wave, cos_theta=cos_theta, flux=flux, viewingangle_dependece=True) # If viewingangle_dependece=False,it won't take into account cos_theta
            transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_)
    elif opts.modelType == "BulladynwindNSBH":
        mdyn = opts.mdyn
        mwin = opts.mwin
        phi = opts.opening_angle
        temp = opts.temp
        phase, wave, cos_theta, flux = mattia_template(dynwind=True, dataDir=opts.possisDir,mdyn=mdyn, mwin=mwin, phi=phi, temp=temp, nsbh=True)
        if phi==0 or phi==90:
            model = model_(phase=phase, wave=wave, cos_theta=cos_theta, flux=flux, viewingangle_dependece=False) # If viewingangle_dependece=False,it won't take into account cos_theta
            transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_notheta_)
        else:
            model = model_(phase=phase, wave=wave, cos_theta=cos_theta, flux=flux, viewingangle_dependece=True) # If viewingangle_dependece=False,it won't take into account cos_theta
            transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_)
    elif opts.modelType == "afterglow":
        def template(E0 = opts.E0, theta_v = (np.arange(opts.thetai,opts.thetaf+1, 1) * np.pi/180) [::-1], n = 1.0, theta_c = np.pi/4.0, theta_w = np.pi/4.0, p = 2.5, epsilon_B = 1e-2, epsilon_E = 1e-1, jetType = 0, specType = 0, ksiN = 1.0, dL = 3.09e19):
            f = []
            for t in theta_v:
                phase, wave, cos_theta, flux = afterglow_template(jetType = jetType,
                                                       specType = specType,
                                                       ksiN = ksiN,
                                                       dL = dL,
                                                       theta_v = t,
                                                       E0 = E0,
                                                       theta_c = theta_c,
                                                       theta_w = theta_w,
                                                       n = n, p = p,
                                                       epsilon_E = epsilon_E,
                                                       epsilon_B = epsilon_B)
                f.append(flux)
            flux = np.transpose(np.array(f), (1, 2, 0))
            cos_theta = np.cos(theta_v)
            return phase, wave, cos_theta, flux

        filename = opts.afterglowDir+'{:.2e}'.format(opts.E0)+'_'+'{:.0f}'.format(opts.thetai)+'_'+'{:.0f}'.format(opts.thetaf)+'_afterglow.pkl'
        if not os.path.exists(filename):
            phase, wave, cos_theta, flux = template()
            cache = {'phase':phase, 'wave': wave, 'cos_theta': cos_theta, 'flux': flux}
            with open(filename,'wb') as f:
                pickle.dump(cache,f)
        else:
            with open(filename,'rb') as f:
                cache = pickle.load(f)

        model = model_(viewingangle_dependece=True, **cache) # If viewingangle_dependece=False,it won't take into account cos_theta
        #model = model_(phase=phase, wave=wave, flux=flux, viewingangle_dependece=True, afterglow=True) # If viewingangle_dependece=False,it won't take into account cos_theta
        transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_)
    elif opts.modelType == "Tophat":
        if opts.doColor:
            phase, wave, flux = tophatgw(bandpass_dir=inputDir, rmag_corr=opts.rmag_corr, imag_corr=opts.imag_corr, diff_gr = opts.diff_gr, diff_gi = opts.diff_gi)
        else:
            phase, wave, flux = tophat_template(bandpass_dir=inputDir, mag0=opts.mag, dmag=opts.dmag)
        model = model_(phase=phase, wave=wave, flux=flux, viewingangle_dependece=False) # If viewingangle_dependece=False,it won't take into account cos_theta
        transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_notheta_)
    elif opts.modelType == "kasen":
        phase, wave, flux = kasen_(dataDir='/Users/anasaguescarracedo/Dropbox/PhD/github_repositories/Kasen_Kilonova_Models_2017/')
        model = model_(phase=phase, wave=wave, flux=flux, viewingangle_dependece=False) # If viewingangle_dependece=False,it won't take into account cos_theta
        transientprop = dict(lcmodel=model, lcsimul_func=random_parameters_notheta_)


    ntransient = opts.ntransient # if ntransient==None, it will simulate LCs following the rate
    if opts.useRate==True:
        ntransient = None
    dec_range = (-90,90)
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
                                          sfd98_dir=sfd98_dir,
                                          skymap=map_struct,
                                          apply_mwebv=True
                                          ) # When using skymap, zrange, dec_range and ra_range are ignored. 


    if opts.doPlots:
        plotName = os.path.join(outputDir,'radec_.pdf')
        plt.figure()
        plt.scatter(tr.ra, tr.dec, alpha=.5)
        plt.xlabel('RA [degrees]')
        plt.ylabel('Declination [degrees]')
        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()

        plotName = os.path.join(outputDir,'field_hist.pdf')
        plt.figure()
        _,_,_ = plt.hist(df["fieldid"])
        plt.xlabel('fieldid')
        plt.ylabel('N')
        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()

    phase_range=(0,opts.dt)
    survey = simsurvey.SimulSurvey(generator=tr, plan=plan, sourcenoise=opts.sourcenoise,
            phase_range=phase_range, n_det=opts.NPoints, threshold=opts.threshold) # By default 2 detections required (n_det=2)
    if survey._properties['instruments'] is None:
        print('There is no hope... setting results to 0.')
        filename = os.path.join(outputDir, 'KNe.dat')
        fid = open(filename, 'w')
        fid.write('%.5e %.5e %.5e' % (0.0, 0.0, 0.0))
        fid.close()
        exit(0)



    lcs = survey.get_lightcurves(progress_bar=True, notebook=True)

    if len(lcs.meta_full) < 1:
        print("No lightcurves... returning.")
        exit(0)

    if opts.doField:
        nKNe = 0
        for lc in lcs:
            infield = False
            tmp = simsurvey.utils.skybins.SurveyField(fields['ra'][fields['field_id']==int(opts.field)],fields['dec'][fields['field_id']==int(opts.field)])
            infield = tmp.coord_in_field(ra=lc.meta['ra'],dec=lc.meta['dec'])['field']
            if infield:
                nKNe += 1
        print("Number of KNe within field %d: %d" % (int(opts.field), nKNe))



    pcklfilename = opts.pickleFile

    pcklFile = os.path.join(outputDir,pcklfilename)
    pickle.dump(lcs,open(pcklFile, "wb" ) )

    if opts.doPlots:
        # plot N vs redshift
        plotName = os.path.join(outputDir,'redshift_.pdf')
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

        '''cmap = cm.summer

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
        plt.close()'''


if opts.doFilter:

    pcklfilename = opts.pickleFile

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

    def filterfunc(lc):
        time = lc[lc['flux']/lc['fluxerr']>=opts.threshold]['time']
        dt_minutes = (time[-1]-time[0])*60*24
        separation_min = opts.separation #minutes
        mask_dt = dt_minutes>=separation_min

        c = SkyCoord(ra=lc.meta['ra'], dec=lc.meta['dec'], unit=(u.degree, u.degree))
        gal_lat_limit = opts.gal    # galactic latitude filter (+- in degrees)
        mask_gal = c.galactic.b.degree>=gal_lat_limit or c.galactic.b.degree <= -gal_lat_limit

        if mask_dt and mask_gal and (lc.meta['idx_orig'] not in idx_orig):
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


    # Fields observed twice:
    fields2 = {}
    fields2['field_id'] = []
    fields2['ra'] = []
    fields2['dec'] = []
    for f in np.unique(df.fieldid):
        if len(np.unique(df[df.fieldid==f].expid))>=2:
            fields2['field_id'].append(f)
            fields2['ra'].append(fields['ra'][fields['field_id']==f][0])
            fields2['dec'].append(fields['dec'][fields['field_id']==f][0])

    fields2['field_id'] = np.array(fields2['field_id'])
    fields2['ra'] = np.array(fields2['ra'])
    fields2['dec'] = np.array(fields2['dec'])

    # Created in the fields observed twice:
    infields2 = 0
    for lc in lcs.meta_full:
        in_ = False
        for f2 in fields2['field_id']:
            tmp = simsurvey.utils.skybins.SurveyField(fields2['ra'][fields2['field_id']==f2],fields2['dec'][fields2['field_id']==f2])
            in_ = tmp.coord_in_field(ra=lc['ra'],dec=lc['dec'])['field']
            if in_:
                break
        if in_:
            infields2 += 1

    print()
    print('From simsurvey:', len(lcs.lcs))
    print('Passing kNe:', len(lcs_filter.lcs))
    print('kNe showing a maximum:', sum(peak_))
    print()
    ntransient = len(lcs.meta_notobserved)+len(lcs.meta_full)
    print('Number of created kNe:', ntransient)
    print('Number of created kNe falling in the covered area:', len(lcs.meta_full['z']))
    print('Number of created kNe falling in the area covered twice:', infields2)

    print('Number of detected over all created:', len(lcs.lcs)/ntransient, ' (efficiency?)')
    print('Number of filtered over all created:', len(lcs_filter.lcs)/ntransient, ' (efficiency?)')

    print('Number of detected over created in the observed fields:', len(lcs.lcs)/len(lcs.meta_full['z']), ' (efficiency)')
    print('Number of filtered over created in the observed fields:', len(lcs_filter.lcs)/len(lcs.meta_full['z']), ' (efficiency)')

    print('Number of detected over created in the fields observed twice:', len(lcs.lcs)/infields2, ' (efficiency)')
    print('Number of filtered over created in the fields observed twice:', len(lcs_filter.lcs)/infields2, ' (efficiency)')


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
        plotName = os.path.join(outputDir,'lcs_'+opts.pickleFile[:-4]+'.pdf')
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
        plt.close()

        # plot N vs redshift
        plotName = os.path.join(outputDir,'redshift_'+opts.pickleFile[:-4]+'.pdf')
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
            plotName = os.path.join(outputDir,'detection_date_'+opts.pickleFile[:-4]+'.pdf')
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
        plt.hist(np.array(tlast5s_) - np.array(t0_))
        plt.xlabel('Last 5sigma detection (d)')
        plt.ylabel('N')
        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()

        cmap = cm.summer

        plotName = os.path.join(outputDir,'radec_'+opts.pickleFile[:-4]+'.pdf')
        plt.figure()
        idx = np.where(fields['limMags'] > 0)[0]
        sc = plt.scatter(fields['ra'][idx], fields['dec'][idx], marker='s', c=fields['limMags'][idx], alpha=.5, cmap=cmap, label='Pointings')
        idx = np.where(fields['limMags'] == 0)[0]
        #lt.scatter(fields['ra'][idx], fields['dec'][idx], marker='s', color='orange', alpha=.5, label='Pointings')
        plt.scatter(tr.ra, tr.dec, label='created', alpha=.5, color='blue')
        if lcs.meta is not None:
            plt.scatter(lcs.meta['ra'], lcs.meta['dec'],marker='*', label='detected', alpha=.5, color='green')
        plt.scatter(ra_, dec_, marker='*', label='filtered', alpha=.5, color='red')
        cbar = plt.colorbar(sc)
        cbar.set_label('Limiting Magnitudes')
        plt.legend(loc=0)
        plt.xlabel('RA [degrees]')
        plt.ylabel('DEC [degrees]')
        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()

        if opts.doColor:
            plotName = os.path.join(outputDir,'radec_redshift_'+opts.pickleFile[:-4]+'.pdf')
            plt.figure()
            idx = np.where(fields['limMags'] > 0)[0]
            plt.scatter(fields['ra'][idx], fields['dec'][idx], marker='s', c=fields['limMags'][idx], alpha=.5, cmap=cmap, label='Pointings')
            idx = np.where(fields['limMags'] == 0)[0]
            #lt.scatter(fields['ra'][idx], fields['dec'][idx], marker='s', color='orange', alpha=.5, label='Pointings')
            # plt.scatter(tr.ra, tr.dec, alpha=.5, c=tr.z, cmap=cmap, label='created')
            if lcs.meta is not None:
                sc = plt.scatter(lcs.meta['ra'], lcs.meta['dec'],marker='*', label='detected', alpha=.5, c=lcs.meta['z'], cmap=cmap)
            # plt.scatter(ra_, dec_, marker='*', label='filtered', alpha=.5, color='red')
            cbar = plt.colorbar(sc)
            cbar.set_label('Redshift')
            plt.legend(loc=0)
            plt.xlabel('RA [degrees]')
            plt.ylabel('DEC [degrees]')
            plt.tight_layout()
            plt.savefig(plotName)
            plt.close()


            plotName = os.path.join(outputDir,'radec_passbands_'+opts.pickleFile[:-4]+'.pdf')
            gband = [lc['band'] == 'ztfg' for lc in lcs.lcs]
            rband = [lc['band'] == 'ztfr' for lc in lcs.lcs]
            iband = [lc['band'] == 'ztfi' for lc in lcs.lcs]
            for i in range(len(lcs.lcs)):
                if len(lcs[i][gband[i]]) > 0: plt.scatter(lcs.meta['ra'][i], lcs.meta['dec'][i], marker='o', color='g', alpha=0.5, label=None)
                if len(lcs[i][rband[i]]) > 0: plt.scatter(lcs.meta['ra'][i], lcs.meta['dec'][i], marker='o', color='r', alpha=0.5, label=None)
                if len(lcs[i][iband[i]]) > 0: plt.scatter(lcs.meta['ra'][i], lcs.meta['dec'][i], marker='o', color='orange', alpha=0.5, label=None)
            gline = Line2D([],[],color='g', marker='o',label='g-band')
            rline = Line2D([],[],color='r', marker='o',label='r-band')
            iline = Line2D([],[],color='orange', marker='o', label='iband')
            plt.legend(loc=0, handles=[gline, rline, iline])
            plt.xlabel('RA [degrees]')
            plt.ylabel('DEC [degrees]')
            plt.tight_layout()
            plt.savefig(plotName)
            plt.close()

            plotName = os.path.join(outputDir,'hist_passbands_'+opts.pickleFile[:-4]+'.pdf')
            gband = [lc['band'] == 'ztfg' for lc in lcs.lcs]
            rband = [lc['band'] == 'ztfr' for lc in lcs.lcs]
            iband = [lc['band'] == 'ztfi' for lc in lcs.lcs]
            gobs = []
            robs = []
            iobs = []
            for i in range(len(lcs.lcs)):
                gobs.append(len(lcs[i][gband[i]]))
                robs.append(len(lcs[i][rband[i]]))
                iobs.append(len(lcs[i][iband[i]]))
            plt.hist(gobs, bins=range(7), histtype='step', color='g', label='g-band', linewidth=2)
            plt.hist(robs, bins=range(7), histtype='step', color='r', label='r-band', linewidth=2)
            plt.hist(iobs, bins=range(7), histtype='step', color='orange', label='i-band', linewidth=2)
            plt.legend(loc=0)
            plt.xlim(0, 7)
            plt.xlabel('Number of observations')
            plt.ylabel('Number of Recovered KNe')
            plt.tight_layout()
            plt.savefig(plotName)
            plt.close()

        if (opts.modelType=='Bulla' or opts.modelType=='Bulladynwind') and phi!=0 and phi!=90:
            # 2D histogram of detections per redshift and viewing angle.
            plotName = os.path.join(outputDir,'2dhist_'+opts.pickleFile[:-4]+'.pdf')
            plt.figure()
            binsz = 20
            binsth = 9
            zmax = lcs.meta_full['z'].max()
            H, yedges, xedges, patches = plt.hist2d(z_, th_,bins=[binsz, binsth], range=[[0,zmax],[0,90]])
            plt.imshow(H,aspect='auto',interpolation='none', origin='low', extent=[yedges.min(), yedges.max(), xedges.min(), xedges.max()])
            plt.colorbar(label='N detections')
            plt.xlabel('z')
            plt.ylabel('viewing angle ($^\circ$)')
            plt.savefig(plotName)
            plt.close()

            # Hist viewing angle
            plotName = os.path.join(outputDir,'viewingangle_'+opts.pickleFile[:-4]+'.pdf')
            plt.figure()
            range_bins=(0,lcs.meta_full['theta'].max())
            plt.hist(lcs.meta_full['theta'],range=range_bins, label='created')
            plt.hist(th_,range=range_bins, label='filtered_kNe')
            plt.legend(loc=0)
            plt.xlabel('viewing angle ($^\circ$)')
            plt.ylabel('N detections')
            plt.tight_layout()
            plt.savefig(plotName)
            plt.close()
