import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import io
import astropy.io.fits as fits
from penquins import Kowalski
import astropy.stats


def plot_LC_cutout(lc_cand,mu_cand,name_cand,t_gw,t_i,t_f,m_min=0,m_max=0,loc='best',output_dir='./'):
    
    ####################################
    # Cutouts from arvo packages
    
    username = 'kped'
    password = 'queryitEDdy!'

    ko = Kowalski(protocol='https', host='kowalski.caltech.edu', port=443, verbose=True,
                 username=username, password=password)
    q = {"query_type": "general_search", "query": "db['ZTF_alerts'].find({'objectId': {'$eq': '"+name_cand+"'}})" }
    r = ko.query(query=q,timeout=30)

    cm = plt.cm.cubehelix
    stamp_ext = 'fits'
    candidate = r['result_data']['query_result'][0]
    # getting cutouts
    sdir = output_dir
    f=open(output_dir+candidate['objectId']+'-sci.fits','wb').write(io.BytesIO(candidate['cutoutScience']['stampData']).getvalue()) 
    f=open(output_dir+candidate['objectId']+'-ref.fits','wb').write(io.BytesIO(candidate['cutoutTemplate']['stampData']).getvalue()) 
    f=open(output_dir+candidate['objectId']+'-dif.fits','wb').write(io.BytesIO(candidate['cutoutDifference']['stampData']).getvalue())

    

    fig = plt.figure(figsize=(7, 9))
    
    # plotting cutouts
    # transient
    ax1 = plt.subplot(222)
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    itype='-sci.'
    img_data = fits.getdata(output_dir+candidate['objectId']+itype+stamp_ext)
    # getting scale values
    imgd_scaled = np.log10(img_data)
    vmax = np.sort(imgd_scaled[28:36,28:36].flatten())[-3]
    npixel = (len(img_data)+1)**2
    imgd_flat = img_data.flatten()
    v_onesig = np.log10(np.nanmedian(img_data) - astropy.stats.mad_std(imgd_flat[np.isnan(imgd_flat)==False])*1.4826)
    vmin= max(v_onesig, np.nanmin(imgd_scaled))
    
    ax1.imshow(imgd_scaled, cmap=cm, vmax=1.015*vmax, vmin=1.005*vmin)
    
    
    # reference
    ax2 = plt.subplot(221)
    ax2.set_xticks([])
    ax2.set_yticks([])
    
    itype='-ref.'
    img_data = fits.getdata(output_dir+candidate['objectId']+itype+stamp_ext)
    imgd_scaled = np.log10(img_data)
    ax2.imshow(imgd_scaled, cmap=cm, vmax=vmax, vmin=vmin)
    
    # plotting location of transient in ref image
    ax2.axhline(y=31,xmin=0.5-0.1,xmax=0.5-0.05, c='white',linewidth=3.75,alpha=0.7)
    ax2.axhline(y=31,xmin=0.5+0.05,xmax=0.5+0.1, c='white',linewidth=3.75,alpha=0.7)
    ax2.axvline(x=31,ymin=0.5-0.1,ymax=0.5-0.05, c='white',linewidth=3.75,alpha=0.7)
    ax2.axvline(x=31,ymin=0.5+0.05,ymax=0.5+0.1, c='white',linewidth=3.75,alpha=0.7)


    ####################################
    # Reading photometric data
    
    time=np.asarray(lc_cand).T[0].astype('float64')
# uncomment to check the loaded time, change this value as needed
#     last_relevant_obs = -4 
#     print('The mean time between observations: ',np.mean((time[1:last_relevant_obs]-time[:last_relevant_obs-1])*24))

    g_band = []
    r_band = []
    i_band = []
    z_band = []


    for obs in lc_cand:
        if obs[1] =='g':
            g_band.append(obs)
        elif obs[1] == 'r':
            r_band.append(obs)
        elif obs[1] == 'i':
            i_band.append(obs)
        elif obs[1] == 'z':
            z_band.append(obs)

    g_band=np.asarray(g_band)
    r_band=np.asarray(r_band)
    i_band=np.asarray(i_band)
    z_band=np.asarray(z_band)

    # plotting LC with two y-axis 
    ax3 = plt.subplot(212)
    ax4 = ax3.twinx()
    print(m_min,m_max)
    if m_min == 0 and m_max == 0:
        lc_temp = np.asarray(lc_cand)
        cond = lc_temp.T[3].astype('float64')==99.0
        m_max = np.amax(lc_temp.T[3].astype('float64')[~cond])+1.6
        m_min = np.amin(lc_temp.T[3].astype('float64')[~cond])-0.2        
        print(m_min,m_max)
    for filt in [g_band,r_band,i_band,z_band]:
        if len(filt)>0:
            if filt[0][1] =='i':
                c='orange'
            elif filt[0][1] == 'z':
                c='k'
            elif filt[0][1] == 'g':
                c='green'
            else:
                c=filt[0][1]

            t_jd =filt.T[0].astype('float64')
            Mag = filt.T[2].astype('float64') 
            mag = filt.T[3].astype('float64')
            e_mag = filt.T[4].astype('float64') 
            lim_mag = filt.T[5].astype('float64') 

            cond = filt.T[3].astype('float64')==99.0
            
            ax3.errorbar(t_jd[~cond]-t_gw,mag[~cond],e_mag[~cond],marker='o',color=c,label=filt[0][1]+'-band')
            ax4.errorbar(t_jd[~cond]-t_gw,Mag[~cond],e_mag[~cond],marker='o',color=c,alpha=0)

            if cond.any():
                ax3.scatter(t_jd[cond]-t_gw,lim_mag[cond],color=c,marker='^',alpha=0.5,s=100,label=filt[0][1]+'-band upper limit')

    ax3.legend(loc=loc,fontsize=12)
    
    ax3.set_xlim(t_i,t_f)
    ax3.set_ylim(m_max,m_min)
    ax4.set_ylim(m_max-mu_cand,m_min-mu_cand)

    ax3.set_xlabel('Days after the merger',fontsize=14)
    ax3.set_ylabel('AB magnitude',fontsize=14)
    ax4.set_ylabel('Absolute magnitude',fontsize=14)

    plt.subplots_adjust(wspace=.1, hspace=0)
    plt.savefig(output_dir+candidate['objectId']+'.pdf',bbox_inches='tight')
    plt.show()


# for 180425z

t_gw = 2458598.846134 #mjd

names=['ZTF19aasckkq','ZTF19aarykkb','ZTF19aasckwd','ZTF19aarzaod']

mu = [16.37+20.44,15.93+19.11,19.15+20.13,14.92+20.47]
t_i = [-1,-0.1,-1,0]
t_f = [4,1.15,1.15,1.15]
m_max = [21,20.6,20.5,22]
m_min = [18.95,17.65,19.8,19.45]
loc=['upper left','lower left','lower left','upper right']


lc=[
# ckkq
[[2458597.9629,'g',99.0,99.0,99.0,20.08,'P48+ZTF'],
[2458598.8976,'g',-16.37,20.44,0.25,20.48,'ZTF'],
[2458598.8993,'r',-16.27,20.55,0.35,20.24,'ZTF'],
[2458598.9722,'r',99.0,99.0,99.0,20.31,'ZTF'],
[2458599.8869,'r',-16.31,20.51,0.20,20.68,'ZTF'],
[2458600.6384,'g',-16.15,20.66,0.18,21.44,'LT'],
[2458600.6395,'r',-16.44,20.37,0.17,21.81,'LT'],
[2458600.6406,'i',-16.53,20.28,0.23,21.73,'LT'],
[2458600.6416,'z',-16.72,20.09,0.38,21.25,'LT'],
[2458600.9291,'r',-16.52,20.29,0.10,20.50,'KPED'],
[2458601.6058,'g',-16.11,20.70,0.09,22.80,'LT'],
[2458601.6087,'r',-16.39,20.42,0.10,21.85,'LT'],
[2458601.6116,'i',-16.60,20.21,0.07,22.36,'LT'],
[2458601.6145,'z',-16.66,20.15,0.15,20.66,'LT'],
[2458601.7544,'r',-16.61,20.20,0.09,21.39,'KPED'],
[2458602.4739,'g',-16.08,20.73,0.24,21.28,'KPED'],
[2458602.4746,'r',-16.71,20.10,0.12,21.18,'KPED'],
[2458605.8970,'r',-16.75,20.06,0.18,20.57,'ZTF'],
[2458608.8767,'r',-17.13,19.68,0.19,20.28,'ZTF'],
[2458608.8994,'g',-16.56,20.25,0.20,20.63,'ZTF'],
[2458612.9206,'r',-16.99,19.83,0.13,20.74,'ZTF']],

#ckkb
[[2458592.9585,'g',99.0,99.0,99.0,18.74,'P48+ZTF'],
[2458598.9353,'g',-15.93,19.11,0.11,19.99,'P48'],
[2458598.9353,'g',-15.93,19.11,0.11,19.99,'P48'],
[2458598.9358,'g',-15.93,19.11,0.11,20.02,'P48'],
[2458598.9523,'g',-15.92,19.12,0.14,19.98,'P48'],
[2458598.9835,'r',-16.41,18.63,0.10,20.04,'P48'],
[2458599.1658,'g',-16.20,18.84,0.04,20.00,'LOT'],
[2458599.1681,'r',-16.81,18.23,0.02,20.00,'LOT'],
[2458599.1705,'i',-17.12,17.92,0.02,20.00,'LOT'],
[2458599.25462297,'r',18.0485-(15.93+19.11),18.0485, 0.1,20.85,'GIT'],
[2458599.41775125,'g',18.64-(15.93+19.11),18.64, 0.27,20.85,'GIT'],
[2458599.44477927,'i',17.786-(15.93+19.11),17.786, 0.09,20.85,'GIT'],
[2458599.5445,'g',-15.95,19.09,0.12,21.82,'LT'],
[2458599.5456,'r',-16.46,18.58,0.03,21.23,'LT'],
[2458599.5466,'i',-16.64,18.40,0.02,21.05,'LT'],
[2458599.5477,'z',-16.60,18.44,0.05,20.30,'LT'],
[2458599.8487,'g',-15.96,19.08,0.05,21.03,'P48'],
[2458599.8487,'g',-15.96,19.08,0.05,21.03,'P48'],
[2458599.9140,'r',-16.53,18.51,0.04,20.70,'P48'],
[2458599.9140,'r',-16.53,18.51,0.04,20.70,'P48'],
[2458606.8320,'r',-16.65,18.38,0.07,20.12,'P48'],
[2458606.9384,'g',-16.00,19.04,0.09,20.50,'P48'],
[2458617.8364,'r',-16.69,18.35,0.08,19.92,'P48'],
[2458617.8958,'g',-15.92,19.12,0.13,19.89,'P48']],

#ckwd
[[2458597.9629,'g',99.0,99.0,99.0,20.01,'P48+ZTF'],
[2458598.8976,'g',-19.15,20.13,0.23,20.26,'P48+ZTF'],
[2458598.8993,'r',99.0,99.0,99.0,20.34,'ztf'],
[2458598.9722,'r',99.0,99.0,99.0,20.44,'P48+ZTF'],
[2458599.7597,'g',-19.16,20.12,0.12,21.05,'P48+ZTF'],
[2458599.7597,'g',-19.16,20.12,0.12,21.05,'P48+ZTF'],
[2458599.8869,'r',-19.12,20.16,0.13,20.62,'P48+ZTF'],
[2458600.9344,'r',99.0,99.0,99.0,29.84,'KPED+KPED']],


# zaod
[[2458598.9353,'g',99.0,99.0,99.0,20.01,'P48+ZTF'],
[2458598.9358,'g',-14.92,20.47,0.34,20.02,'P48+ZTF'],
[2458598.9523,'g',99.0,99.0,99.0,20.02,'P48+ZTF'],
[2458598.9835,'r',-15.28,20.11,0.18,20.10,'P48+ZTF'],
[2458598.9835,'r',-15.28,20.11,0.18,20.10,'P48+ZTF'],
[2458599.2374,'i',-15.46,19.93,0.13,21.00,'LOT+LOT'],
[2458599.2419,'g',-13.94,21.45,0.47,22.00,'LOT+LOT'],
[2458599.2456,'r',-15.13,20.26,0.25,21.00,'LOT+LOT'],
[2458599.5756,'r',-15.24,20.15,0.05,22.01,'LT+IOO'],
[2458599.5820,'i',-15.84,19.56,0.03,22.14,'LT+IOO'],
[2458599.6512,'g',99.0,99.0,99.0,20.32,'LT+IOO'],
[2458599.6523,'r',-15.20,20.19,0.07,21.20,'LT+IOO'],
[2458599.6533,'i',-15.82,19.57,0.04,21.69,'LT+IOO'],
[2458599.9140,'r',-15.29,20.10,0.14,20.76,'P48+ZTF']]
]

# plotting

for i in range(len(names)):
    plot_LC_cutout(lc[i],mu[i],names[i],t_gw,t_i[i],t_f[i],output_dir='../output/GW190425/')#,m_min=m_min[i],m_max=m_max[i],loc=loc[i]'upper left')


