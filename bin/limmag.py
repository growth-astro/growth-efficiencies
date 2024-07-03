
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
matplotlib.rcParams.update({'font.size': 24})
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
from matplotlib import patches
import matplotlib.pyplot as plt
from matplotlib import rcParams

from astropy.visualization.wcsaxes import SphericalCircle
from astropy.table import Table

__version__ = 1.0

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-f","--fields", help="Observed fields.", default='../data/GW190425/ZTF_fields_bayestar.dat,../data/GW190425/Gattini_fields_bayestar.dat')
    parser.add_option("-o", "--output", help="output file",default="../output/GW190425/limmag.pdf")

    parser.add_option("-t", "--telescopes", help="Telescopes.", default ="ZTF,Gattini")
    parser.add_option("-g", "--gps", help="Event time GPS.", default=1240215503.011549, type=float)

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
telescopes = opts.telescopes.split(",")

start_time = Time(opts.gps, format='gps').mjd

if opts.telescopes == "DECam":
    names = ('time', 'filt', 'limmag')
    tables = {}
    for telescope, field in zip(telescopes,opts.fields.split(",")):
        table = Table.read(field, format='ascii', names=names)
        table["mjd"] = Time(table["time"], format='mjd').mjd - start_time
        table.sort('mjd')

        table["exposure_time"] = 300.0
        table["filter_id"] = table["filt"].copy()
        table["filter_id"][table["filter_id"] == "i"] = 3.0
        table["filter_id"][table["filter_id"] == "z"] = 4.0
        table["filter_id"] = table["filter_id"].astype(float)

        table_group = table.group_by('mjd').groups
        tables[telescope] = table_group

else:
    names = ('field_id', 'filter_id', 'time', 'limmag', 'exposure_time')
    tables = {}
    for telescope, field in zip(telescopes,opts.fields.split(",")):
        table = Table.read(field, format='ascii', names=names)
        table["mjd"] = Time(table["time"], format='isot').mjd - start_time
        table.sort('mjd')
        tables[telescope] = table

if len(telescopes) == 2:
    
    fig, ax = plt.subplots(2, 3, figsize=(16, 12))
    limits = [[0,0.2],[0.8,1.2],[1.8,2.2]]
    
    for ii in range(3):
    
        lim = limits[ii]
        idx = np.where(tables["ZTF"]["filter_id"] == 1)[0]
        ax[0,ii].plot(tables["ZTF"]["mjd"][idx], tables["ZTF"]["limmag"][idx], 'gv', label = 'g-band', markerfacecolor='none')
        idx = np.where(tables["ZTF"]["filter_id"] == 2)[0]
        ax[0,ii].plot(tables["ZTF"]["mjd"][idx], tables["ZTF"]["limmag"][idx], 'rv', label = 'r-band', markerfacecolor='none')
        ax[0,ii].set_xlim(lim)
        ax[0,ii].set_ylim([19,23])
        ax[0,ii].invert_yaxis() 
        ax[0,ii].grid()
    
        plt.setp(ax[0,ii].get_xticklabels(), visible=False)
        if not ii == 0:
            plt.setp(ax[0,ii].get_yticklabels(), visible=False)
    
        if ii == 2:
            ax[0,ii].legend()
    
        ax[1,ii].plot(tables["Gattini"]["mjd"][idx], tables["Gattini"]["limmag"][idx]+0.9, 'kv', label = 'J-band')
        ax[1,ii].set_xlim(lim)
        ax[1,ii].set_ylim([13,16])
        ax[1,ii].invert_yaxis()
        ax[1,ii].grid()
    
        if not ii == 0:
            plt.setp(ax[1,ii].get_yticklabels(), visible=False)
    
        if ii == 2:
            ax[1,ii].legend()

elif opts.telescopes == "DECam":

    fig, ax = plt.subplots(1, 3, figsize=(14, 7))
    #limits = [[0.4,0.8],[1.3,2.0],[2.5,3.0]]  
    limits = [[0.35,0.47],[1.3,1.55],[2.3,2.55]]

    for ii in range(3):

        lim = limits[ii]

        for row in tables["DECam"]:
            limmag = np.median(row['limmag'])
            if row["filter_id"][0] == 3:
                ax[ii].plot(row["mjd"][0], limmag, 'v', color='yellow', label = 'i-band', markerfacecolor='yellow', markersize=12, alpha=0.5)
            else:
                ax[ii].plot(row["mjd"][0], limmag, 'v', color='black', label = 'z-band', markerfacecolor='black', markersize=12, alpha=0.5)

        ax[ii].set_xlim(lim)
        ax[ii].set_ylim([19.5,23.0])
        #ax[ii].set_ylim([20.5,23.5])
        ax[ii].invert_yaxis()
        ax[ii].grid()

        if not ii == 0:
            plt.setp(ax[ii].get_yticklabels(), visible=False)

else:

    fig, ax = plt.subplots(1, 3, figsize=(14, 7))
    #limits = [[0.4,0.8],[1.3,2.0],[2.5,3.0]]  
    limits = [[0.05,0.17],[0.90,1.1],[3.05,3.2]] 

    for ii in range(3):
   
        lim = limits[ii]
        idx = np.where((tables["ZTF"]["filter_id"] == 1) & (tables["ZTF"]["exposure_time"] == 30))[0]
        ax[ii].plot(tables["ZTF"]["mjd"][idx], tables["ZTF"]["limmag"][idx], 'gv', label = 'g-band', markerfacecolor='none', markersize=12, alpha=0.5)
        idx = np.where((tables["ZTF"]["filter_id"] == 2) & (tables["ZTF"]["exposure_time"] == 30))[0]
        ax[ii].plot(tables["ZTF"]["mjd"][idx], tables["ZTF"]["limmag"][idx], 'rv', label = 'r-band', markerfacecolor='none', markersize=12, alpha=0.5)

        idx = np.where((tables["ZTF"]["filter_id"] == 1) & (tables["ZTF"]["exposure_time"] > 30))[0]
        ax[ii].plot(tables["ZTF"]["mjd"][idx], tables["ZTF"]["limmag"][idx], 'gv', label = 'g-band', markersize=12, alpha=0.5)
        idx = np.where((tables["ZTF"]["filter_id"] == 2) & (tables["ZTF"]["exposure_time"] > 30))[0]
        ax[ii].plot(tables["ZTF"]["mjd"][idx], tables["ZTF"]["limmag"][idx], 'rv', label = 'r-band', markersize=12, alpha=0.5)

        ax[ii].set_xlim(lim)
        #ax[ii].set_ylim([18,22])
        ax[ii].set_ylim([19.5,22.5])
        ax[ii].invert_yaxis()
        ax[ii].grid()

        if not ii == 0:
            plt.setp(ax[ii].get_yticklabels(), visible=False)

fig.text(0.5, 0.01, 'Time from trigger [days]', ha='center', fontsize=30)
fig.text(0.015, 0.5, 'Limiting Magnitude [AB]', va='center', rotation='vertical', fontsize=30)

plt.subplots_adjust(wspace=0.15, hspace=0.07)

outputDir = "/".join(opts.output.split("/")[:-1])
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)
    
plt.savefig(opts.output,dpi=200,bbox_inches='tight')    
plt.close()    

