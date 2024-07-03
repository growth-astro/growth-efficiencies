#!/usr/bin/env python
# coding: utf-8

# In[2]:


import matplotlib
#matplotlib.rc('text', usetex=True)
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14, 'legend.fontsize': 10})
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import numpy as np
import os
from astropy.io import ascii
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u
from astropy.table import Table
from astropy.time import Time, TimeDelta
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
from astropy.io import fits
from astropy.io import ascii
from penquins import Kowalski
import pandas as pd
import datetime
import json
import requests
import glob
import os


# In[31]:




k = Kowalski(username='sanand', password='PrettyStr0ngPa$$w0rd')

# k = Kowalski(host='localhost', port=8000, protocol='http', username='admin', password='admin')

connection_ok = k.check_connection()
print(f'Connection OK: {connection_ok}')

collection_ZTF_alerts = 'ZTF_alerts'
collection_ZTF_alerts_aux = 'ZTF_alerts_aux'

def k_query(name,radec=False):
    q = {"query_type": "find",
         "query": {
             "catalog": collection_ZTF_alerts,
             "filter": {"objectId": name},
             "projection": {"_id": 0, "cutoutScience": 0, "cutoutTemplate": 0, "cutoutDifference": 0},
         }
         }
    r = k.query(query=q)
    alerts = r['result_data']['query_result']
    ra,dec = alerts[0]['candidate']['ra'],alerts[0]['candidate']['dec']
    #print(f'Number of alerts for {name}:', len(alerts))
    dfs = []

    for packet in alerts:
        df = pd.DataFrame(packet['candidate'], index=[0])
        dfs.append(df)

    # drop duplicate entries. decide using jd
    dfs = pd.concat(dfs, ignore_index=True, sort=False).drop_duplicates(subset='jd').reset_index(drop=True)

    q = {"query_type": "find_one",
         "query": {
             "catalog": collection_ZTF_alerts_aux,
             "filter": {"_id": name},
             "projection": {"_id": 0},
         }
         }
    r = k.query(query=q)
    alert_aux = r['result_data']['query_result']

    df_prv = pd.DataFrame(alert_aux['prv_candidates'])
    dflc = pd.concat([dfs, df_prv],
                     ignore_index=True,
                     sort=False).drop_duplicates(subset='jd').reset_index(drop=True).sort_values(by=['jd'])

    dflc, lc = assemble_lc(dflc,name)
    
    if radec:
        return ra,dec
    else:
        return dflc, lc
    # display(dflc)
    
def is_star(dflc, match_radius_arcsec=1.5, star_galaxy_threshold=0.4):
    try:
        return (dflc.iloc[-1].distpsnr1 < match_radius_arcsec) & (dflc.iloc[-1].sgscore1 > star_galaxy_threshold)
    except Exception as _e:
        print(_e)
        return False
    
def assemble_lc(dflc, objectId, composite=False, match_radius_arcsec=1.5, star_galaxy_threshold=0.4):
    # mjds:
    dflc['mjd'] = dflc.jd - 2400000.5

    dflc['datetime'] = dflc['mjd'].apply(lambda x: Time(x, format='mjd').datetime)
    
    # strings:
    dflc['dt'] = dflc['datetime'].apply(lambda x: x.strftime('%Y-%m-%d %H:%M:%S'))

    dflc.sort_values(by=['mjd'], inplace=True)

    # fractional days ago
    dflc['days_ago'] = dflc['datetime'].apply(lambda x:
                                              (datetime.datetime.utcnow() - x).total_seconds() / 86400.)

    if is_star(dflc, match_radius_arcsec=match_radius_arcsec, star_galaxy_threshold=star_galaxy_threshold):
        # print('It is a star!')
        # variable object/star? take into account flux in ref images:
        lc = []

        # fix old alerts:
        dflc.replace('None', np.nan, inplace=True)

        # prior to 2018-11-12, non-detections don't have field and rcid in the alert packet,
        # which makes inferring upper limits more difficult
        # fix using pdiffimfilename:
        w = dflc.rcid.isnull()
        if np.sum(w):
            dflc.loc[w, 'rcid'] = dflc.loc[w, 'pdiffimfilename'].apply(lambda x:
                                                      ccd_quad_2_rc(ccd=int(os.path.basename(x).split('_')[4][1:]),
                                                                    quad=int(os.path.basename(x).split('_')[6][1:])))
            dflc.loc[w, 'field'] = dflc.loc[w, 'pdiffimfilename'].apply(lambda x:
                                                                        int(os.path.basename(x).split('_')[2][1:]))

        grp = dflc.groupby(['fid', 'field', 'rcid'])
        impute_magnr = grp['magnr'].agg(lambda x: np.median(x[np.isfinite(x)]))
        impute_sigmagnr = grp['sigmagnr'].agg(lambda x: np.median(x[np.isfinite(x)]))

        for idx, grpi in grp:
            w = np.isnan(grpi['magnr'])
            w2 = grpi[w].index
            dflc.loc[w2, 'magnr'] = impute_magnr[idx]
            dflc.loc[w2, 'sigmagnr'] = impute_sigmagnr[idx]

        # fix weird isdiffpos'es:
        w_1 = dflc['isdiffpos'] == '1'
        dflc.loc[w_1, 'isdiffpos'] = 't'

        dflc['sign'] = 2 * (dflc['isdiffpos'] == 't') - 1
        
        # Eric Bellm 20190722: Convert to DC magnitudes (see p.102 of the Explanatory Supplement)
        dflc['dc_flux'] = 10**(-0.4*dflc['magnr']) + dflc['sign'] * 10**(-0.4*dflc['magpsf'])
        w_dc_flux_good = dflc['dc_flux'] > 0
        dflc.loc[w_dc_flux_good, 'dc_mag'] = -2.5 * np.log10(dflc.loc[w_dc_flux_good, 'dc_flux'])
        dflc.loc[w_dc_flux_good, 'dc_sigmag'] = np.sqrt(
            (10**(-0.4*dflc['magnr'])* dflc['sigmagnr']) **2. + 
            (10**(-0.4*dflc['magpsf']) * dflc['sigmapsf'])**2.) / dflc.loc[w_dc_flux_good, 'dc_flux']
        
        dflc['dc_flux_ulim'] = 10**(-0.4*dflc['magnr']) + 10**(-0.4*dflc['diffmaglim'])
        dflc['dc_flux_llim'] = 10**(-0.4*dflc['magnr']) - 10**(-0.4*dflc['diffmaglim'])
        
        w_dc_flux_ulim_good = dflc['dc_flux_ulim'] > 0
        w_dc_flux_llim_good = dflc['dc_flux_llim'] > 0
        
        dflc.loc[w_dc_flux_ulim_good, 'dc_mag_ulim'] = -2.5 * np.log10(10**(-0.4*dflc.loc[w_dc_flux_ulim_good, 'magnr']) + 
                                                                       10**(-0.4*dflc.loc[w_dc_flux_ulim_good, 'diffmaglim']))
        dflc.loc[w_dc_flux_llim_good, 'dc_mag_llim'] = -2.5 * np.log10(10**(-0.4*dflc.loc[w_dc_flux_llim_good, 'magnr']) - 
                                                                       10**(-0.4*dflc.loc[w_dc_flux_llim_good, 'diffmaglim']))

        # if some of the above produces NaNs for some reason, try fixing it sloppy way:
        for fid in (1, 2, 3):
            if fid in dflc.fid.values:
                ref_flux = None
                w = (dflc.fid == fid) & ~dflc.magpsf.isnull() & (dflc.distnr <= match_radius_arcsec)
                if np.sum(w):
                    ref_mag = np.float64(dflc.loc[w].iloc[0]['magnr'])
                    ref_flux = np.float64(10 ** (0.4 * (27 - ref_mag)))
                    # print(fid, ref_mag, ref_flux)

                wnodet_old = (dflc.fid == fid) & dflc.magpsf.isnull() &                              dflc.dc_mag_ulim.isnull() & (dflc.diffmaglim > 0)

                if np.sum(wnodet_old) and (ref_flux is not None):
                    # if we have a non-detection that means that there's no flux +/- 5 sigma from
                    # the ref flux (unless it's a bad subtraction)
                    dflc.loc[wnodet_old, 'difference_fluxlim'] = 10 ** (0.4 * (27 - dflc.loc[wnodet_old, 'diffmaglim']))
                    dflc.loc[wnodet_old, 'dc_flux_ulim'] = ref_flux + dflc.loc[wnodet_old, 'difference_fluxlim']
                    dflc.loc[wnodet_old, 'dc_flux_llim'] = ref_flux - dflc.loc[wnodet_old, 'difference_fluxlim']

                    # mask bad values:
                    w_u_good = (dflc.fid == fid) & dflc.magpsf.isnull() &                                dflc.dc_mag_ulim.isnull() & (dflc.diffmaglim > 0) & (dflc.dc_flux_ulim > 0)
                    w_l_good = (dflc.fid == fid) & dflc.magpsf.isnull() &                                dflc.dc_mag_ulim.isnull() & (dflc.diffmaglim > 0) & (dflc.dc_flux_llim > 0)

                    dflc.loc[w_u_good, 'dc_mag_ulim'] = 27 - 2.5 * np.log10(dflc.loc[w_u_good, 'dc_flux_ulim'])
                    dflc.loc[w_l_good, 'dc_mag_llim'] = 27 - 2.5 * np.log10(dflc.loc[w_l_good, 'dc_flux_llim'])

        # corrections done, now proceed with assembly
        for fid in (1, 2, 3):
            # print(fid)
            # get detections in this filter:
            w = (dflc.fid == fid) & ~dflc.magpsf.isnull()
            lc_dets = pd.concat([dflc.loc[w, 'jd'], dflc.loc[w, 'dt'], 
                                 dflc.loc[w, 'datetime'], dflc.loc[w, 'days_ago'],
                                 dflc.loc[w, 'mjd'], dflc.loc[w, 'dc_mag'], dflc.loc[w, 'dc_sigmag']],
                                axis=1, ignore_index=True, sort=False) if np.sum(w) else None
            if lc_dets is not None:
                lc_dets.columns = ['jd', 'dt',  'datetime', 'days_ago', 'mjd', 'mag', 'magerr']

            wnodet = (dflc.fid == fid) & dflc.magpsf.isnull()
            # print(wnodet)

            lc_non_dets = pd.concat([dflc.loc[wnodet, 'jd'], dflc.loc[wnodet, 'dt'], 
                                     dflc.loc[wnodet, 'datetime'], dflc.loc[wnodet, 'days_ago'],
                                     dflc.loc[wnodet, 'mjd'], dflc.loc[wnodet, 'dc_mag_llim'],
                                     dflc.loc[wnodet, 'dc_mag_ulim']],
                                    axis=1, ignore_index=True, sort=False) if np.sum(wnodet) else None
            if lc_non_dets is not None:
                lc_non_dets.columns = ['jd', 'dt', 'datetime', 'days_ago', 'mjd', 'mag_llim', 'mag_ulim']

            if lc_dets is None and lc_non_dets is None:
                continue

            lc_joint = None

            if lc_dets is not None:
                # print(lc_dets)
                # print(lc_dets.to_dict('records'))
                lc_joint = lc_dets
            if lc_non_dets is not None:
                # print(lc_non_dets.to_dict('records'))
                lc_joint = lc_non_dets if lc_joint is None else pd.concat([lc_joint, lc_non_dets],
                                                                          axis=0, ignore_index=True, sort=False)

            # sort by date and fill NaNs with zeros
            lc_joint.sort_values(by=['mjd'], inplace=True)
            # print(lc_joint)
            lc_joint = lc_joint.fillna(0)

            # single or multiple alert packets used?
            lc_id = f"{objectId}_composite_{datetime.datetime.utcnow().strftime('%Y%m%d_%H%M%S')}"                 if composite else f"{objectId}_{int(dflc.loc[0, 'candid'])}"
            # print(lc_id)

            lc_save = {"telescope": "PO:1.2m",
                       "instrument": "ZTF",
                       "filter": fid,
                       "source": "alert_stream",
                       "comment": "corrected for flux in reference image",
                       "id": lc_id,
                       "lc_type": "temporal",
                       "data": lc_joint.to_dict('records')
                       }
            lc.append(lc_save)

    else:
        # print('Not a star!')
        # not a star (transient): up to three individual lcs
        lc = []

        for fid in (1, 2, 3):
            # print(fid)
            # get detections in this filter:
            w = (dflc.fid == fid) & ~dflc.magpsf.isnull()
            lc_dets = pd.concat([dflc.loc[w, 'jd'], dflc.loc[w, 'dt'], 
                                 dflc.loc[w, 'datetime'], dflc.loc[w, 'days_ago'],
                                 dflc.loc[w, 'mjd'], dflc.loc[w, 'magpsf'], dflc.loc[w, 'sigmapsf']],
                                axis=1, ignore_index=True, sort=False) if np.sum(w) else None
            if lc_dets is not None:
                lc_dets.columns = ['jd', 'dt', 'datetime', 'days_ago', 'mjd', 'mag', 'magerr']

            wnodet = (dflc.fid == fid) & dflc.magpsf.isnull()

            lc_non_dets = pd.concat([dflc.loc[wnodet, 'jd'], dflc.loc[wnodet, 'dt'], 
                                     dflc.loc[wnodet, 'datetime'], dflc.loc[wnodet, 'days_ago'],
                                     dflc.loc[wnodet, 'mjd'], dflc.loc[wnodet, 'diffmaglim']],
                                    axis=1, ignore_index=True, sort=False) if np.sum(wnodet) else None
            if lc_non_dets is not None:
                lc_non_dets.columns = ['jd', 'dt', 'datetime', 'days_ago', 'mjd', 'mag_ulim']

            if lc_dets is None and lc_non_dets is None:
                continue

            lc_joint = None

            if lc_dets is not None:
                lc_joint = lc_dets
            if lc_non_dets is not None:
                lc_joint = lc_non_dets if lc_joint is None else pd.concat([lc_joint, lc_non_dets],
                                                                          axis=0, ignore_index=True, sort=False)

            # sort by date and fill NaNs with zeros
            lc_joint.sort_values(by=['mjd'], inplace=True)
            lc_joint = lc_joint.fillna(0)

            # single or multiple alert packets used?
            lc_id = f"{objectId}_composite_{datetime.datetime.utcnow().strftime('%Y%m%d_%H%M%S')}"                 if composite else f"{objectId}_{int(dflc.loc[0, 'candid'])}"

            lc_save = {"telescope": "PO:1.2m",
                       "instrument": "ZTF",
                       "filter": fid,
                       "source": "alert_stream",
                       "comment": "no corrections applied. using raw magpsf, sigmapsf, and diffmaglim",
                       "id": lc_id,
                       "lc_type": "temporal",
                       "data": lc_joint.to_dict('records')
                       }
            lc.append(lc_save)

    return dflc, lc


# In[64]:


def stack_lc(tbl, days_stack=1., snt_det=3, snt_ul=5):
    """Given a dataframe with a maxlike light curve,
    stack the flux """

    if 'jdobs' in list(tbl.colnames):
        key_jd = 'jdobs'
    elif 'jd' in list(tbl.colnames):
        key_jd = 'jd'
    else:
        print("What is the column for the JD??")
        pdb.set_trace()
    t_out = Table([[],[],[],[],[],[],[],[],[]],
                  names=(key_jd, 'flux', 'flux_unc', 'zp', 'ezp',
                         'mag', 'mag_unc', 'limmag', 'filter'),
                  dtype=('double', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'S'))
    # Bin separately by filter
    filters = list(set(tbl['filter']))
    for f in filters:

        t = tbl[tbl['filter'] == f]

        bins = np.arange(int(np.max(t[key_jd]) - np.min(t[key_jd]))+2)
        dt0 = np.min(t[key_jd]) - int(np.min(t[key_jd]))
        if dt0 <= 0.4:
            start = int(np.min(t[key_jd])) - 0.6
        else:
            start = int(np.min(t[key_jd])) + 0.4
        bins = bins + start
        for b in bins:
            temp = t[(t[key_jd] > b) & (t[key_jd] < b+1)]
            if len(temp) == 0:
                continue
            new_jd = np.mean(np.array(temp[key_jd]))

            if len(set(temp['zp'])) == 1:
                zp = temp['zp'][0]
                flux = np.array(temp['Flux_maxlike'])
                flux_unc = np.array(temp['Flux_maxlike_unc'])
                flux[np.isnan(flux)] = 0
                # Use weights only if there are only detections
                if np.min(flux/flux_unc) >= snt_det:
                    weights = (flux/flux_unc)**2
                    new_flux = np.sum(np.array(temp['Flux_maxlike'])*weights)/np.sum(weights)
                else:
                    new_flux = np.mean(np.array(temp['Flux_maxlike']))
                new_flux_unc = np.sqrt(np.sum(np.array(temp['Flux_maxlike_unc'])**2))/len(temp)
            else:
                zp = temp['zp'][0]
                flux1 = np.array(temp['Flux_maxlike'])
                flux1_unc = np.array(temp['Flux_maxlike_unc'])
                zp1 = np.array(temp['zp'])
                flux = 10**((2.5*np.log10(flux1) - zp1 + zp) / 2.5)
                flux_unc = 10**((2.5*np.log10(flux1_unc) - zp1 + zp) / 2.5)
                flux[np.isnan(flux)] = 0
                # Use weights only if there are only detections
                if np.min(flux/flux_unc) >= snt_det:
                    weights = (flux/flux_unc)**2
                    new_flux = np.sum(flux*weights)/np.sum(weights)
                else:
                    new_flux = np.mean(flux)
                new_flux_unc = np.sqrt(np.sum(flux_unc**2))/len(temp)
            if new_flux/new_flux_unc > snt_det:
                mag_stack = -2.5*np.log10(new_flux) + zp
                mag_unc_stack = np.abs(-2.5*np.log10(new_flux-new_flux_unc) + 2.5*np.log10(new_flux))
                maglim_stack = 99.
            else:
                mag_stack = 99.
                mag_unc_stack = 99.
                maglim_stack = -2.5 * np.log10(snt_ul * new_flux_unc) + zp
            ezp = np.sum(temp['ezp']**2)/len(temp)
            t_out.add_row([new_jd, new_flux, new_flux_unc, zp, ezp, mag_stack,
                           mag_unc_stack, maglim_stack, f])

    return t_out


# In[ ]:





# In[70]:


basename = '/Users/tahumada/Downloads/GRB200826/'
jd_trigger = 2459087.6868055556
files = glob.glob(basename+'*')
names = [s[s.find('ZTF'):12+s.find('ZTF')] for s in files]


# In[67]:





# In[68]:



program1=True

# path1 = basename+'fitsfiles/'
# files1 = os.listdir(path1)
for file in files:
#     name = file[11:23]
    name = file[file.find('ZTF'):12+file.find('ZTF')]
    hdu = fits.open(file)

    data = Table()
    data.add_column(hdu[1].data['jdobs'],name= 'jd')
    data.add_column(hdu[1].data['mag'],name= 'magpsf')
    data.add_column(hdu[1].data['mag_unc'],name= 'mag_unc')
    data.add_column(hdu[1].data['zp'],name= 'zp')
    data.add_column(hdu[1].data['ezp'],name= 'ezp')
    data.add_column(hdu[1].data['Flux_maxlike'],name= 'Flux_maxlike')
    data.add_column(hdu[1].data['Flux_maxlike_unc'],name= 'Flux_maxlike_unc')
    
    ul = hdu[1].data['limmag']
    ul[hdu[1].data['mag']<98]=-10
    data.add_column(ul,name= 'ul')
    data.add_column(hdu[1].data['filter'],name= 'filter')

    nondet = data['ul'] > 0.
    times = Time(data['jd'], format='jd')
    t0 = Time(2459087.6868055556, format='jd')
    days = TimeDelta(times - t0, scale='tt')
    upperlim = data[nondet]
    dayslim = days[nondet]
    data = data[nondet^True] # remove the upper limits from the data
    days = days[nondet^True]
    gband = data['filter'] == 'g'
    rband = data['filter'] == 'r'
    iband = data['filter'] == 'i'
    glim = upperlim['filter'] == 'g'
    rlim = upperlim['filter'] == 'r'
    ilim = upperlim['filter'] == 'i'
    cond_g = days[gband].value < 30
    cond_r = days[rband].value < 30
    cond_i = days[iband].value < 30
    condlim_g = dayslim[glim].value < 30
    condlim_r = dayslim[rlim].value < 30
    condlim_i = dayslim[ilim].value < 30

    plt.figure()
    plt.errorbar(days[gband].value[cond_g], data[gband]['magpsf'][cond_g], yerr=data[gband]['mag_unc'][cond_g], fmt='o', color='g', label='g-band detection')
    plt.errorbar(days[rband].value[cond_r], data[rband]['magpsf'][cond_r], yerr=data[rband]['mag_unc'][cond_r], fmt='o', color='r', label='r-band detection')
    plt.errorbar(days[iband].value[cond_i], data[iband]['magpsf'][cond_i], yerr=data[iband]['mag_unc'][cond_i], fmt='o', color='goldenrod', label='i-band detection')
    plt.scatter(dayslim[glim].value[condlim_g], upperlim[glim]['ul'][condlim_g], marker='v',alpha=0.3,s=100, color='g', label='g-band upperlim')
    plt.scatter(dayslim[rlim].value[condlim_r], upperlim[rlim]['ul'][condlim_r], marker='v',alpha=0.3,s=100, color='r', label='r-band upperlim')
    plt.scatter(dayslim[ilim].value[condlim_i], upperlim[ilim]['ul'][condlim_i], marker='v',alpha=0.3,s=100, color='goldenrod', label='i-band upperlim')
    
    if program1:
        k_dflc,k_lc = k_query(name)

        colors = {1: 'g', 2: 'r', 3: 'goldenrod'}
        filter_names = {1: "ZTF g", 2: "ZTF r", 3: "ZTF i"}
        t_format='datetime'
        k_time = 0
        mag_tot = np.array([])
        time_tot = np.array([])
        for lc in k_lc:
            fid = lc['filter']
            df_lc = pd.DataFrame.from_records(lc['data'])
    #         display(df_lc)
            # mags:
            if 'mag' in df_lc:
                w = df_lc['mag'] > 1e-6
                if np.sum(w) > 0:
                    t = df_lc.loc[w, t_format]
                    mag = df_lc.loc[w, 'mag']
                    mag_error = df_lc.loc[w, 'magerr']
                    cond = TimeDelta(Time(t) - t0, scale='tt').value < 30
                    plt.errorbar(TimeDelta(Time(t) - t0, scale='tt').value[cond], mag[cond], yerr=mag_error[cond],fmt='o',c=colors[fid], label=f'{filter_names[fid]}')
                    k_time = max(np.amax(TimeDelta(Time(t) - t0, scale='tt').value[cond]),k_time)
                    time_tot = np.hstack((time_tot,TimeDelta(Time(t) - t0, scale='tt').value)) 
                    mag_tot = np.hstack((mag_tot,mag.values)) 
            # upper limits:
            if 'mag_ulim' in df_lc:
                w = df_lc['mag_ulim'] > 1e-6
                if np.sum(w) > 0:
                    t = df_lc.loc[w, t_format]
                    mag_ulim = df_lc.loc[w, 'mag_ulim']
                    plt.scatter(TimeDelta(Time(t) - t0, scale='tt').value, mag_ulim, marker='v', s=100, alpha=0.3,c=colors[fid], label=f'{filter_names[fid]} upper limit')
                    time_tot = np.hstack((time_tot,TimeDelta(Time(t) - t0, scale='tt').value))
                    mag_tot = np.hstack((mag_tot,np.ones(len(mag_ulim))*99))
    ymin_fp = 0
    ymax_fp = 100
    tmax_fp = 0
    mag_tot = np.asarray(mag_tot).flatten()
    time_tot = np.asarray(time_tot).flatten()
    try:
        ymin_fp = np.amax(data['magpsf']) + 0.5
        ymax_fp = np.amin(data['magpsf']) - 0.5
        tmax_fp = np.amax(np.array(days.value[days.value<30])) + 3.0
    except ValueError:
        print(name, ' has no detections using forced photometry')
    if len(mag_tot[(mag_tot<98) * (time_tot>0)])>0:
        ymin = max(ymin_fp,np.amax(mag_tot[(mag_tot<98) * (time_tot>0)])) + 0.5
        ymax = min(ymax_fp,np.amin(mag_tot[(mag_tot<98) * (time_tot>0)])) - 0.5
    else:
        ymin = ymin_fp + 0.5
        ymax = ymax_fp - 0.5    
    tmax = max(tmax_fp,k_time) + 3.0
    plt.xlabel('Days after the merger')
    plt.ylabel('AB magnitude')
    plt.ylim(ymin, ymax)
    plt.xlim(-10., tmax)
    plt.title(name,fontsize=14,fontweight='bold')
    plt.grid()
    if not os.path.exists('plots'):
        os.system('mkdir plots')
        
    # stacking
    stck = stack_lc(data)
    
    nondet = stck['limmag'] != 99.0
    times = Time(stck['jd'], format='jd')
    t0 = Time(jd_trigger, format='jd')
    days = TimeDelta(times - t0, scale='tt')
    upperlim = stck[nondet]
    dayslim = days[nondet]
    data = stck[nondet^True] # remove the upper limits from the data
    days = days[nondet^True]
    gband = stck['filter'] == 'g'
    rband = stck['filter'] == 'r'
    iband = stck['filter'] == 'i'
    glim = upperlim['filter'] == 'g'
    rlim = upperlim['filter'] == 'r'
    ilim = upperlim['filter'] == 'i'
    cond_g = days[gband].value < 30
    cond_r = days[rband].value < 30
    cond_i = days[iband].value < 30
    condlim_g = dayslim[glim].value < 30
    condlim_r = dayslim[rlim].value < 30
    condlim_i = dayslim[ilim].value < 30
    if (days<0).any(): print(name, 'has a predetection from stacked forced photometry', days)
    plt.errorbar(days[gband].value[cond_g], stck[gband]['mag'][cond_g], yerr=stck[gband]['mag_unc'][cond_g], fmt='x', color='g', label='g-band detection')
    plt.errorbar(days[rband].value[cond_r], stck[rband]['mag'][cond_r], yerr=stck[rband]['mag_unc'][cond_r], fmt='x', color='r', label='r-band detection')
    plt.errorbar(days[iband].value[cond_i], stck[iband]['mag'][cond_i], yerr=stck[iband]['mag_unc'][cond_i], fmt='x', color='goldenrod', label='i-band detection')

    plt.savefig('plots/stacked_'+name+'_lc.pdf')
    t_tot = np.hstack((time_tot,days.value))
#     np.amax(t_tot[t_tot<0]),t_tot
    if len(t_tot[t_tot<0])> 0: print(name,'last detection (hrs)',np.round(np.amax(t_tot[t_tot<0]),3)*-24)
    else: print(name,'first detection!')
    # plot_force_lcs.py
    # # Mostrando plot_force_lcs.py.


# In[ ]:





# In[ ]:




