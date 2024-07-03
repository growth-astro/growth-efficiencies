"""Code for running Yuhan Yao's ForcePhotZTF code on individual ZTF transients. Need to have ForcePhotZTF installed: https://github.com/yaoyuhan/ForcePhotZTF.
Output fits files and plots saved by default in ForcePhotZTF directory."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.io import fits
import pandas as pd
import sys
import argparse
from multiprocessing import cpu_count

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--name", help="transient ZTF name", type=str, default="ZTF20abwysqy")
    parser.add_argument("-s", "--start_time", help="start date for forced photometry (ISOT)", type=str, default=None)
    parser.add_argument("-e", "--end_time", help="end date for forced photometry (ISOT)", type=str, default=None)
    parser.add_argument("-p", "--path", help="path to the ForcePhotZTF directory", type=str, default="/Users/shreyaanand/ZTF/growth-too-papers/")
    parser.add_argument("-r", "--ra", help="target RA (deg)", type=float, default=6.7856144)
    parser.add_argument("-d", "--dec", help="target Dec (deg)", type=float, default=34.0273218)
    parser.add_argument("-a", "--SNR", help="signal-to-noise ratio for detection", type=float, default=5.5)
    parser.add_argument("--doPlots", help="plot forced photometry", action="store_true", default=False)
    opts = parser.parse_args()
    return opts

opts = parse_commandline()

sys.path.append(opts.path)
from ForcePhotZTF.keypairs import get_keypairs
from ForcePhotZTF.force_lc import download_marshal_lightcurve, get_coo_ZTFtarget, astrometry_spread, \
                                    get_cutout_data, download_images_diffpsf_refdiff
from ForcePhotZTF.force_mcmc import get_forced_phot_mcmc
from ForcePhotZTF.refine_lc import get_recerence_jds, plotlcs, read_ipac_lc, read_mcmc_lc
from ForcePhotZTF.force_maxlike import get_forced_phot_maaxlike


def do_forcephot(targetdir, name, ra, dec, start_jd = None, end_jd = None, r_psf = 3, r_bkg_in = 10, r_bkg_out = 15, SNT=5.5, verbose = False):
    lightcurve_file = targetdir + 'lightcurves/marshal_lc_'+name+'.csv'
    download_images_diffpsf_refdiff(targetdir, ra, dec, start_jd, end_jd)
    get_cutout_data(name, targetdir, ra, dec, r_psf = r_psf, 
                r_bkg_in = r_bkg_in, r_bkg_out = r_bkg_out, verbose = False)
    get_recerence_jds(name, targetdir, only_partnership=False, retain_iband = True,
                  oldsuffix = '_info.fits', newsuffix = '_info_ref.fits', verbose=True)
    get_forced_phot_maaxlike(name, targetdir, ra, dec, SNT = SNT, r_psf = r_psf, r_bkg_in = r_bkg_in, r_bkg_out = r_bkg_out, verbose = verbose)
    return None

ncpu = cpu_count()
print("{0} CPUs".format(ncpu))

name = opts.name
targetdir = opts.path + 'ForcePhotZTF/' + name + '/'
ra, dec = opts.ra, opts.dec

try:
    download_marshal_lightcurve(name, targetdir)
except KeyError:
    pass

if opts.start_time is not None:
    start_jd = Time(opts.start_time, format='isot').jd
    end_jd = Time(opts.end_time, format='isot').jd
else:
    start_jd = None
    end_jd = None

do_forcephot(targetdir, name, ra, dec, start_jd=start_jd, end_jd=end_jd, SNT=opts.SNR)

f = fits.open(targetdir+'force_phot_'+name+'_maxlikelihood_lc.fits')
data = f[1].data
f.close()

if opts.doPlots:
    times = Time(data['jdobs'], format='jd')
    daysago = TimeDelta(Time.now() - times, scale='tt')
    nondet = data['mag'] > 98.
    upperlim = data[nondet]
    dayslim = daysago[nondet]
    data = data[nondet^True] # remove the upper limits from the data
    daysago = daysago[nondet^True]
    gband = data['filter'] == 'g'
    rband = data['filter'] == 'r'
    iband = data['filter'] == 'i'
    glim = upperlim['filter'] == 'g'
    rlim = upperlim['filter'] == 'r'
    ilim = upperlim['filter'] == 'i'
    plt.figure()
    plt.plot(dayslim[glim].value, upperlim[glim]['limmag'], 'v', fillstyle='none', color='g', label='g-band upperlim')
    plt.plot(dayslim[rlim].value, upperlim[rlim]['limmag'], 'v', fillstyle='none', color='r', label='r-band upperlim')
    plt.plot(dayslim[ilim].value, upperlim[ilim]['limmag'], 'v', fillstyle='none', color='goldenrod', label='i-band upperlim')
    plt.errorbar(daysago[rband].value, data[rband]['mag'], yerr=data[rband]['mag_unc'], fmt='.', color='r', label='r-band')
    plt.errorbar(daysago[iband].value, data[iband]['mag'], yerr=data[iband]['mag_unc'], fmt='.', color='goldenrod', label='i-band')
    plt.errorbar(daysago[gband].value, data[gband]['mag'], yerr=data[gband]['mag_unc'], fmt='.', color='g', label='g-band')
    plt.xlabel('days ago')
    plt.ylabel('apparent magnitude')
    plt.grid()
    plt.legend()
    plt.savefig(targetdir+name+'_forcedphot.png')