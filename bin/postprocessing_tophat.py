from __future__ import division
import matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('xtick', labelsize=16)
matplotlib.rc('ytick', labelsize=16)
matplotlib.rc('font', **font)
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy.time import Time
from astropy.io import ascii
import pickle
import argparse
import sys
import os
import requests
import glob
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from pandas import read_csv, DataFrame
import healpy as hp
from ligo.skymap import postprocess, distance
from ligo.skymap.io import fits
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outputDir", help="output file",default="../output/")
    parser.add_argument('--doPlotEfficiency', help="plot three efficiency curves for each event", action='store_true', default=False)
    parser.add_argument('--doCondor', help="parse files created by condor", action='store_true', default=False)
    parser.add_argument('--doComposite', help='plot all efficiency curves on the same axes', action='store_true', default=False)
    parser.add_argument('--doFiducial', help='plot fiducial curve along with the other efficiency curves for each event', action='store_true', default=False)
    parser.add_argument('--doLuminosityFunc', help='calculate the luminosity function over all events', action='store_true', default=False)
    parser.add_argument('--doTophatDecay', help='calculate luminosity function with linear model', action='store_true', default=False)
    opts = parser.parse_args()
    return opts

def get_absmag(appmag, extinction, z=None, d=None, redshift=True):
    if redshift==True:
        d = cosmo.luminosity_distance(z).to(u.pc)
    else: 
        d = d.to(u.pc)
    M = float(appmag - (5 * np.log10( (d.to(u.pc) / (10 * u.pc)).value) ) - (3.1*extinction))
    return M


opts = parse_commandline()

#absmag = np.linspace(-20,-10, 50)
absmag = np.linspace(-12, -19, num=50)
event_names = ['GW190425_O4', 'GW200105_O4', 'GW200115_O4', 'S230521k', 'S230528a', 'S230529ay', 'S230615az', 'S230627c', 'S231029k', 'S231113bw']
#event_names = ['GW190425', 'S190426c','GW190814', 'S190901ap', 'S190910d', 'S190910h', 'S190923y', 'S190930t', 'S191205ah', 'S191213g', 'S200105ae', 'S200115j', 'S200213t']
#event_times = ['2019-04-25T08:18:05', '2019-04-26T15:21:55', '2019-08-14T21:10:39', '2019-09-01T23:31:01', '2019-09-10T01:26:19', '2019-09-10T08:29:58', \
#    '2019-09-23T12:55:59', '2019-09-30T14:34:07', '2019-12-05T21:52:08', '2019-12-13T04:34:08', '2020-01-05T16:24:26', '2020-01-15T04:23:09', '2020-02-13T04:10:40']
isbns = np.array([True, False, False, True, False, True, False, False, False, True, False, False, True])
lowFAR = np.array([True, True, True, False, False, False, False, False, False, True, True, True, True])
color = ['blue', 'gold', 'black', 'dodgerblue', 'firebrick', 'c', 'peru', 'saddlebrown', 'goldenrod', 'indigo', 'r', 'orange', 'blueviolet']
apparent_mags = np.array([21.5, 21.5, 21.0, 21.0, 20.3, 20.4, 20.1, 21.1, 17.9, 20.4, 20.2, 20.8, 21.2])
event_prob = np.array([0.2413, 0.5233, 0.8857, 0.5694, 0.3299, 0.3326, 0.3899, 0.5063, 0.0568, 0.2750, 0.5239, 0.2221, 0.7217])
event_prob2 = np.array([0.203, 0.094, 0.783, 0.465, 0.257, 0.236, 0.082, 0.0485, 0.368, 0.118, 0.377, 0.125, 0.472])
#extinction = [0.06, 0.54, 0.00, 0.05, 0.10, 0.09, 0.08, 0.01, 0.06, 0.23, 0.05, 0.24, 0.11]
extinction = [0.03, 0.34, 0.02, 0.03, 0.04, 0.08, 0.09, 0.05, 0.04, 0.30, 0.05, 0.13, 0.19]
event_distances = [156,377,267,241,632,230,438,108,385,201,283,340,201]*u.Mpc

# multiply the E(B-V) by 3.1 and overplot 
M = np.zeros(len(event_names))
absmin = np.zeros(len(event_names))
absmax = np.zeros(len(event_names))
magerr = np.zeros(len(event_names))

knsim = {}

for i, name in enumerate(event_names):
    if opts.doCondor:
        # construct a dictionary for all events for easy access; loop through condor directories
        path = opts.outputDir +name+'/'
        #path = '../output/'+name+'/simsurvey/tophat_decay/'
        dirs = os.listdir(path)
        if 'condor' in dirs: dirs.remove('condor')
        if 'combined' in dirs: dirs.remove('combined')
        if 'efficiency_tophat.pdf' in dirs: dirs.remove('efficiency_tophat.pdf')
        if name+'_efficiencies.txt' in dirs: dirs.remove(name+'_efficiencies.txt')
        mags = np.array([float(str(dir).split('_')[1]) for dir in dirs])
        dmags = np.array([float(str(dir).split('_')[2]) for dir in dirs])
        mask  = dmags == 0.0
        dirs = np.array(dirs)[mask]
        inds = np.argsort(mags[mask])
        dirs = dirs[inds]
        knsim[name]={'det_efficiency': [0]*len(absmag), 'obs_efficiency': [0]*len(absmag), 'filt_efficiency': [0]*len(absmag)}
        det_efficiency_all = np.ones(len(absmag))
        obs_efficiency_all = np.ones(len(absmag))
        filt_efficiency_all = np.ones(len(absmag))
        for j, outputdir in enumerate(dirs):
            f1 = open(path+outputdir+'/sim.pkl', 'rb')
            lcs = pickle.load(f1) 
            f1.close()               
            knsim[name]['det_efficiency'][j] = len(lcs.lcs)/10000
            knsim[name]['obs_efficiency'][j] = len(lcs.lcs)/len(lcs.meta_full['z'])
            try:
                f2 = open(path+outputdir+'/sim_filter.pkl', 'rb')
                filtered_lcs = pickle.load(f2)
                f2.close()
                knsim[name]['filt_efficiency'][j] = len(filtered_lcs.lcs)/10000
            except:
                knsim[name]['filt_efficiency'][j] = 0.0

    elif opts.doComposite or opts.doLuminosityFunc or opts.doPlotEfficiency:
        # construct a dictionary for all events for easy access, loop through local files
        knsim[name]={'det_efficiency': [0]*len(absmag), 'obs_efficiency': [0]*len(absmag), 'filt_efficiency': [0]*len(absmag)}
        outputdir = '../output/' + name + '/simsurvey/tophat/'
        det_efficiency_all = np.ones(len(absmag))
        obs_efficiency_all = np.ones(len(absmag))
        filt_efficiency_all = np.ones(len(absmag))
        for j, mag in enumerate(absmag):
            mag = round(mag, 1)
            f1 = open(outputdir+str(mag)+'_sim.pkl', 'rb')
            lcs = pickle.load(f1) 
            f1.close()               
            knsim[name]['det_efficiency'][j] = len(lcs.lcs)/10000
            knsim[name]['obs_efficiency'][j] = len(lcs.lcs)/len(lcs.meta_full['z'])
            try:
                f2 = open(outputdir+str(mag)+'_sim_filter.pkl', 'rb')
                filtered_lcs = pickle.load(f2)
                f2.close()
                knsim[name]['filt_efficiency'][j] = len(filtered_lcs.lcs)/len(lcs.meta_full['z'])
            except:
                knsim[name]['filt_efficiency'][j] = 0.0

#    path = '../data/'+name+'/'
#    files = os.listdir(path)
#    filename = [fname for fname in files if fname.startswith('LAL')]
#    if len(filename) < 1: 
#        print('no LALInference skymap available')
#        filename = [fname for fname in files if fname.startswith('bay')]

#    apparent = apparent_mags[i]

#    skymap, metadata = fits.read_sky_map(path+filename[0], nest=True, distances=True)
#    distmu = skymap[1]
#    distsigma = skymap[2]
#    distnorm = skymap[3]
#    skymap = skymap[0]
#    nside = hp.npix2nside(len(skymap))

    if opts.doPlotEfficiency:
        plt.figure(i)

        if opts.doFiducial:

            s = requests.Session()
            retries = Retry(total=5,
                            backoff_factor=0.1,
                            status_forcelist=[ 500, 502, 503, 504 ])
            s.mount('http://', HTTPAdapter(max_retries=retries))

            loc = EarthLocation.of_site('Palomar')

            #calculate the moments from distmu, distsigma and distnorm
            mom_mean, mom_std, mom_norm = distance.parameters_to_moments(distmu,distsigma)
            distmod = 5.0*np.log10(mom_mean*1e6)-5.0
            distmod_min = 5.0*np.log10((mom_mean-mom_std)*1e6)-5.0
            distmod_max = 5.0*np.log10((mom_mean+mom_std)*1e6)-5.0

            abs_mean = np.inf*np.ones(mom_mean.shape)
            abs_min = np.inf*np.ones(mom_mean.shape)
            abs_max = np.inf*np.ones(mom_mean.shape)

            if 'ZTF_fields.dat' in os.listdir(path):
                names = ('field_id', 'filter_id', 'time', 'limmag')
                try:
                    table = Table.read(path+'/ZTF_fields.dat', format='ascii', names=names, data_start=0)
                except:
                    names = ('field_id', 'filter_id', 'time', 'limmag', 'exptime')
                    table = Table.read(path+'/ZTF_fields.dat', format='ascii', names=names, data_start=0)
                table = table.group_by('field_id').groups
                for row in table:
                    field = row['field_id'][0]
                    filter_id = row['filter_id'][0]
                    limmag = np.min(row['limmag'])

                    #r = s.get(f'http://skipper.caltech.edu:8081/telescope/%s/field/{field}/json' % ('ZTF'))
                    res = r.json()

                    coords = np.asarray(res['geometry']['coordinates']).T * u.deg
                    coords = SkyCoord(np.squeeze(coords).T, unit='deg',
                                      location=loc,
                                      obstime=Time(row['time'],format='isot'))

                    xyz = []
                    for r, d in zip(coords.ra.deg, coords.dec.deg):
                        xyz.append(hp.ang2vec(r, d, lonlat=True))
                    xyz = xyz[:4]
                    try:
                        ipix = hp.query_polygon(nside, np.array(xyz), nest=True)
                    except:
                        ipix = []
                    abs_mean[ipix] = limmag - distmod[ipix]
                    abs_min[ipix] = limmag - distmod_min[ipix]
                    abs_max[ipix] = limmag - distmod_max[ipix]
            else:
                print('no ZTF.fields.dat file for %s' %(name))
                abs_mean = apparent - distmod
                abs_min = apparent - distmod_min
                abs_max = apparent - distmod_max

            abs_mean[~np.isfinite(abs_mean)] = np.nanmax(abs_mean)
            abs_min[~np.isfinite(abs_min)] = np.nanmax(abs_min)
            abs_max[~np.isfinite(abs_max)] = np.nanmax(abs_max)

            nside_map = hp.npix2nside(len(skymap))
            sort_idx = np.argsort(skymap)[::-1]
            csm = np.empty(len(skymap))
            csm[sort_idx] = np.cumsum(skymap[sort_idx])

            cls = 100 * postprocess.find_greedy_credible_levels(skymap)
            bins = np.arange(-20,-10, 0.2)
            idx = np.where(np.isfinite(abs_mean))[0]
            hist_mean, bin_edges, patches = plt.hist(abs_mean[idx], bins=bins, weights=skymap[idx], density=True, cumulative=-1, histtype='step', color='k')
            hist_min, bin_edges, patches = plt.hist(abs_min[idx], bins=bins, weights=skymap[idx], density=True, cumulative=-1, histtype='step', color='k')
            hist_max, bin_edges, patches = plt.hist(abs_max[idx], bins=bins, weights=skymap[idx], density=True, cumulative=-1, histtype='step', color='k')
            bins = (bin_edges[1:]+bin_edges[:-1])/2.0
            plt.fill_between(bins, hist_max, hist_mean, facecolor='lightgray')
            plt.fill_between(bins, hist_mean, hist_min, facecolor='lightgray')

        figname = opts.outputDir+'/'+name+'/simsurvey/tophat/efficiency_tophat.pdf'
        plt.scatter(absmag, knsim[name]['det_efficiency'], color='m', marker='o')
        plt.scatter(absmag, knsim[name]['obs_efficiency'], color='goldenrod', marker='o')
        plt.scatter(absmag, knsim[name]['filt_efficiency'], color='g', marker='o')
        mline = Line2D([],[],color='m', marker='o',label='detected in 90% C.V.')
        yline = Line2D([],[],color='goldenrod', marker='o', label='detected in obs. area')
        gline = Line2D([],[],color='green', marker='o', label='filtered in obs. area') 
        #kline = Line2D([],[],color='k', marker='o', label='fraction deeper than absolute mag')
        plt.xlim(-10.+0.5,-20.-0.5)
        plt.ylim(0.0, 1.0)
        plt.xlabel('absolute magnitude')
        plt.ylabel('KN detection efficiency')
        plt.legend(handles=[mline, yline, gline])
        plt.grid()
        plt.savefig(figname)

    if opts.doComposite or opts.doLuminosityFunc:
        distmean, diststd = distance.parameters_to_marginal_moments(skymap,
                                                   distmu,
                                                   distsigma)
        distmean = distmean*u.Mpc
        diststd = diststd*u.Mpc
        M[i] = get_absmag(apparent_mags[i], extinction[i], d=event_distances[i], redshift=False)
        absmin[i] = get_absmag(apparent_mags[i], extinction[i], d=distmean+2*diststd, redshift=False)
        absmax[i] = get_absmag(apparent_mags[i], extinction[i], d=distmean-2*diststd, redshift=False)
        print(name, event_distances[i], M[i])
#magerr = np.abs([M-absmax, M-absmin]).reshape(2,len(M))

if opts.doComposite:

    matplotlib.rcParams.update({'xtick.labelsize':14})
    matplotlib.rcParams.update({'ytick.labelsize':14})

    fig, ax = plt.subplots(1,2, sharey=True, figsize=(10,9))
    #fig.tight_layout(pad=2.0)
    for i, name in enumerate(event_names):
        if isbns[i]==True: j = 0
        else: j = 1
        #err = np.array([magerr[0][i], magerr[1][i]]).reshape(2,1)
        ax[j].plot(absmag, knsim[name]['det_efficiency'], color=color[i], label=name)
        #ax[j].errorbar(M[i]-extinction[i], 0.5*event_prob[i], xerr=err, marker='*', color=color[i], label=None)
    # make a plot showing all of the luminosity functions together
    ax[0].set_xlim(-10,-20)
    ax[1].set_xlim(-10,-20)
    ax[0].set_ylim(0.0, 1.0)
    ax[0].set_xlabel('Absolute Magnitude')
    ax[1].set_xlabel('Absolute Magnitude')
    ax[0].set_ylabel('Recovery Efficiency')
    ax[0].legend(fontsize=14)
    ax[1].legend(fontsize=14)
    ax[0].grid()
    ax[1].grid()
    plt.savefig('../output/composite_detected_efficiency.pdf')

    fig, ax = plt.subplots(1,2, sharey=True, figsize=(8,5))
    fig.tight_layout(pad=3.0)
    for i, name in enumerate(event_names):
        if isbns[i]==True: j = 0
        else: j = 1
        #err = np.array([magerr[0][i], magerr[1][i]]).reshape(2,1)
        ax[j].plot(absmag, knsim[name]['obs_efficiency'], color=color[i], label=name)
        #ax[j].errorbar(M[i]-extinction[i], 0.5, xerr=err, marker='*', color=color[i], label=None)
    # make a plot showing all of the luminosity functions together
    ax[0].set_xlim(-10,-20)
    ax[1].set_xlim(-10,-20)
    ax[0].set_ylim(0.0, 1.0)
    ax[0].set_xlabel('absolute magnitude')
    ax[1].set_xlabel('absolute magnitude')
    ax[0].set_ylabel('KN detection efficiency')
    ax[0].legend(fontsize=14)
    ax[1].legend(fontsize=14)
    ax[0].grid()
    ax[1].grid()
    plt.savefig('../output/composite_observed_efficiency.pdf')

    fig, ax = plt.subplots(1,2, sharey=True, figsize=(8,5))
    fig.tight_layout(pad=3.0)
    for i, name in enumerate(event_names):
        if isbns[i]==True: j = 0
        else: j = 1
        #err = np.array([magerr[0][i], magerr[1][i]]).reshape(2,1)
        ax[j].plot(absmag, knsim[name]['filt_efficiency'], color=color[i], label=name)
        #ax[j].errorbar(M[i]-extinction[i], 0.5*event_prob2[i], xerr=err, marker='*', color=color[i], label=None)
    # make a plot showing all of the luminosity functions together
    ax[0].set_xlim(-10,-20)
    ax[1].set_xlim(-10,-20)
    ax[0].set_ylim(0.0, 1.0)
    ax[0].set_xlabel('absolute magnitude')
    ax[1].set_xlabel('absolute magnitude')
    ax[0].set_ylabel('KN detection efficiency')
    ax[0].legend(fontsize=14)
    ax[1].legend(fontsize=14)
    ax[0].grid()
    ax[1].grid()
    plt.savefig('../output/composite_filtered_efficiency.pdf')

if opts.doTophatDecay:
    # efficiency = np.ones(1050)
    # for name in event_names:
    #     path = '../output/' + name + '/simsurvey/tophat_decay/'
    #     filename = name+'_efficiencies.txt'
    #     t = Table(ascii.read(path+filename, delimiter=' '))
    #     names = ['col1', 'col2', 'col3']
    #     renames = ['mag', 'dmag', 'efficiency']
    #     t.rename_columns(names, renames)
    #     efficiency = efficiency*(np.ones(len(t)) - t['efficiency'])
    # efficiency = np.ones(len(t)) - efficiency
#    mags_unique = np.linspace(-12, -19, num=50)
#    dmags_unique = np.arange(-0.5, 1.6, 0.1)
    Z = np.ones((21, 50))
    ones = np.ones(Z.shape)
    for name in event_names:
        path = opts.outputDir + name + '/'
#        path = opts.outputDir + name + '/simsurvey/tophat_decay/'
        outputDirs = glob.glob(os.path.join(path, '*_*_*_*_*'))
        mags, dmags = [], []
        sim = np.empty(0, float)
        for outputDir in outputDirs:
            outsplit = list(filter(None,outputDir.split("/")[-1].split("_")))
            jd, mag, dmag, phi, theta = np.array(outsplit, dtype=float)
            mags.append(mag)
            dmags.append(dmag)
            try:
                f1 = open(outputDir+'/sim.pkl', 'rb')
                lcs = pickle.load(f1)
                f1.close()
                data_out = len(lcs.lcs)/10000
            except:
                print('file not found in', outputDir)
                data_out = 0.0
            sim = np.append(sim, data_out)
        mags_unique = np.unique(mags)
        dmags_unique = np.unique(dmags)
        
        for ii, mag in enumerate(mags_unique):
            for jj, dmag in enumerate(dmags_unique):
                if dmag < 0:
                    mag_zero = mag - 3*dmag
                    diff = np.abs(mags - mag_zero)
                    mindiff = np.amin(diff)
                    idx = np.where((dmag == dmags) & (diff == mindiff))[0]
                    #idx = np.where((mag_zero == mags) & (dmag == dmags))[0]
                    if len(idx) == 0: continue
                    if mag_zero > -10.0: continue
                else:
                    idx = np.where((mag == mags) & (dmag == dmags))[0]
                    if len(idx) == 0: continue

                Z[jj,ii] *= (1 - sim[idx][0])


    X, Y = np.meshgrid(mags_unique, dmags_unique)    
    Z = ones - Z
    Z[Z==0.0] = 1e-3
    imag = mags_unique == -16.5
    jdmag = dmags_unique == 1.0
    imag = mags_unique	== -16.0
    jdmag = dmags_unique == 1.0

    #models = [['SN Ia', -0.200, -0.186, -19.00, -16.61],
              #['SN IIP', -0.157, -0.166, -16.84, -17.03],
              #['SN IIb', -0.169, -0.261, -17.44, -17.32],
              #['SN IIn', 0.105, 0.115, -17.58, -15.32],
              #['SN Ib', -0.142, -0.146, -17.4, -17.25],
              #['SN Ic', -0.120, -0.143, -17.48, -14.92],
    models = [['GW170817', 1.0, 1.0, -16.6, -16.6]]
    # mags_unique = np.unique(t['mag'])
    # dmags_unique = np.unique(t['dmag'])
    # Z = efficiency.reshape(len(mags_unique), len(dmags_unique)).transpose()
    # tophat_lumfunc = np.array(Z[dmags_unique == 0.0])
    # X, Y = np.meshgrid(mags_unique, dmags_unique)

    #diffs = [[0.44,-0.15], [0.4,0.25],
            # [0.47,0.32], #[0.5,-0.25],
            # [0.4,-0.2], [0.45, 0.3],
    diffs = [[0.45, 0.05]]

    plotName = os.path.join(opts.outputDir,'O4a_overall_tophat_decay.pdf')
    fig = plt.figure(figsize=(8,6))
    c = plt.pcolor(X,Y,Z,vmin=1e-3,vmax=1.0,cmap='coolwarm')
    cbar = fig.colorbar(c)
    plt.xlabel('Peak Magnitude')
    plt.ylabel('Decay Rate (mag day$^{-1}$)')
    for model, diff in zip(models,diffs):

        plt.plot(model[3], model[1], marker='*', markerfacecolor='none',
                    markeredgecolor='white', markersize=30)

        diffx, diffy = diff
        x, y = model[3], model[1]

        x1text = x*diffx*2.0
        y1text = y*diffy*12.0

        x1arrow = x*diffx*2.2
        y1arrow = y*diffy*13.0

        plt.text(x1text, y1text, model[0], fontsize=20)
        plt.annotate('', xy=(x,y), xytext=(x1arrow,y1arrow), 
            arrowprops=dict(facecolor='black', arrowstyle='->'),
                    )

    cbar.set_label('Efficiency')
    plt.ylim([0,1])
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.savefig(plotName,bbox_inches='tight')
    plt.close()

if opts.doTophatGW:
    plt.figure()
    for i, name in enumerate(event_names):
        path = '../output/'+name+'/simsurvey/tophat/color/'
        f = open(path+'sim.pkl', 'rb')
        lcs = pickle.load(f)
        f.close()
        try:
            f1 = open(path+'sim.pkl', 'rb')
            filtered_lcs = pickle.load(f1)
            f1.close()
            nfilt = len(filtered_lcs.lcs)/10000
        except:
            nfilt = 0
        ndet = len(lcs.lcs)/10000
        nobs = len(lcs.lcs)/len(lcs.meta_full['z'])

        plt.plot(i+1, ndet, 'r')
        plt.plot(i+1, nobs, 'b')
        plt.plot(i+1, nfilt, 'goldenrod')

    plt.set_xticklabels(event_names)
    plt.savefig('GW170817_empirical_tophat.pdf')

if opts.doLuminosityFunc:
    overall = np.ones((len(M), len(absmag)))
    BNS = np.ones((len(M), len(absmag)))
    NSBH = np.ones((len(M), len(absmag)))
    highsig_bns = np.ones((len(M), len(absmag)))   
    highsig_nsbh = np.ones((len(M), len(absmag)))

    for i, mag in enumerate(M):
        bns = event_prob[np.where((isbns == True) & (M >= mag))[0]]
        nsbh = event_prob[np.where((isbns == False) & (M >= mag))[0]]
        promising_bns = event_prob[np.where((lowFAR == True) & (isbns == True) & (M >= mag))[0]]
        promising_nsbh = event_prob[np.where((lowFAR == True) & (isbns == False) & (M >= mag))[0]]
        probs = event_prob[M >= mag]
 
        for prob in probs:
            overall[0][i] = overall[0][i]*(1-prob)
        for prob in bns: 
            BNS[0][i] = BNS[0][i]*(1-prob)
        for prob in nsbh: 
            NSBH[0][i] = NSBH[0][i]*(1-prob)
        for prob in promising_bns: 
            highsig_bns[0][i] = highsig_bns[0][i]*(1-prob)
        for prob in promising_nsbh: 
            highsig_nsbh[0][i] = highsig_nsbh[0][i]*(1-prob)

    for i, name in enumerate(event_names):
        overall[1] = overall[1]*(np.ones(len(absmag)) - knsim[name]['det_efficiency'])
        if isbns[i]==True:
            BNS[1] = BNS[1]*(np.ones(len(absmag)) - knsim[name]['det_efficiency'])
        elif isbns[i]==False:
            NSBH[1] = NSBH[1]*(np.ones(len(absmag)) - knsim[name]['det_efficiency'])
        if (lowFAR[i]==True) and (isbns[i]==True):
            highsig_bns[1] = highsig_bns[1]*(np.ones(len(absmag)) - knsim[name]['det_efficiency'])
        elif (lowFAR[i]==True) and (isbns[i]==False):
            highsig_nsbh[1] = highsig_nsbh[1]*(np.ones(len(absmag)) - knsim[name]['det_efficiency'])

        # obs_efficiency_all = obs_efficiency_all*(np.ones(len(absmag)) - knsim[name]['obs_efficiency'])
        # filt_efficiency_all = filt_efficiency_all*(np.ones(len(absmag)) - knsim[name]['filt_efficiency'])

    ind = np.argsort(M)

    t = Table([absmag, overall[1], NSBH[1], BNS[1]], names=['M', 'overall', 'NSBH', 'BNS'])
    ascii.write(t, '../output/simsurvey_flat_efficiencies.dat', overwrite=True)

    t2 = Table([M[ind], overall[0][ind], NSBH[0][ind], BNS[0][ind]], names=['M', 'overall', 'NSBH', 'BNS'])
    ascii.write(t2, '../output/median_flat_efficiencies.dat', overwrite=True)

    # composite luminosity function for KNe
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.step(M[ind], overall[0][ind], 'goldenrod', where='mid', linewidth=2.0, label=None)
    ax.step(M[ind], NSBH[0][ind], 'r', where='mid', label=None, linewidth=2.0)
    ax.step(M[ind], BNS[0][ind], 'b', where='mid', label=None, linewidth=2.0)
    # ax.step(absmag, highsig_bns[0], 'c', label=None)
    # ax.step(absmag, highsig_nsbh[0], 'm', label=None)


    ax.plot(absmag, overall[1], 'o', color='goldenrod', label=None)
    ax.plot(absmag, NSBH[1], 'ro', label=None)
    ax.plot(absmag, BNS[1], 'bo', label=None)
    # ax.plot(absmag, highsig_bns[1], 'co', label=None)
    # ax.step(absmag, highsig_nsbh[1], 'mo', label=None)


    gline = Line2D([],[],color='goldenrod', marker='o',label='overall')
    rline = Line2D([],[],color='r', marker='o',label='NSBH')
    bline = Line2D([],[],color='b', marker='o', label='BNS')
    mline = Line2D([],[],color='m', marker='o', label='high significance NSBH')
    cline = Line2D([],[],color='c', marker='o', label='high significance BNS')
    ax.set_yscale('log')
    ax.set_xlim(-9.8, -20.2)
    plt.ylim(1e-2,1.1)
    plt.xlabel('Absolute Magnitude')
    plt.ylabel('Probability of Zero Detections')
    plt.legend(handles=[gline, rline, bline])
    plt.grid()
    plt.savefig('../output/luminosity_function.pdf')

if opts.doEfficiencySum:
    for i, name in enumerate(event_names):
        
